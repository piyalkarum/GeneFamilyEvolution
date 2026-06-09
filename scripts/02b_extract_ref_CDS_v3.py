#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


TRANSCRIPT_FEATURES = {
    "mRNA", "mrna", "transcript", "RNA",
    "pseudogenic_transcript", "lnc_RNA", "ncRNA", "tRNA", "rRNA"
}


def parse_args():
    ap = argparse.ArgumentParser(
        description="Run Helixer on each gene in family FASTAs and extract predicted CDS."
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--in-fasta", help="Single gene-family FASTA file")
    src.add_argument("--in-dir", help="Directory of gene-family FASTA files")

    ap.add_argument("--out-dir", required=True, help="Output directory for CDS family FASTAs")
    ap.add_argument("--suffix", default=".fasta", help="Input FASTA suffix when using --in-dir")
    ap.add_argument("--lineage", default="land_plant", help="Helixer lineage")
    ap.add_argument("--species", required=True, help="Helixer species name, e.g. Arabidopsis_thaliana")
    ap.add_argument("--helixer-bin", default="helixer", help="Helixer executable name/path")
    ap.add_argument("--keep-temp", action="store_true", help="Keep per-record temp files for debugging")
    return ap.parse_args()


def parse_gff_attributes(attr_str):
    attrs = {}
    for item in attr_str.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs


def run_helixer(input_fa: Path, out_gff: Path, helixer_bin: str, lineage: str, species: str) -> None:
    cmd = [
        helixer_bin,
        "--lineage", lineage,
        "--fasta-path", str(input_fa),
        "--species", species,
        "--gff-output-path", str(out_gff),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            f"Helixer failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"STDOUT:\n{res.stdout}\n"
            f"STDERR:\n{res.stderr}"
        )
    if not out_gff.exists():
        raise RuntimeError(f"Helixer finished but GFF was not created: {out_gff}")


def parse_helixer_gff(gff_file: Path):
    transcripts = {}
    cds_by_transcript = defaultdict(list)

    with open(gff_file, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            if line.startswith(">"):
                break

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = parts
            if not start.isdigit() or not end.isdigit():
                continue

            start = int(start)
            end = int(end)
            attrd = parse_gff_attributes(attrs)

            if feature in TRANSCRIPT_FEATURES:
                tx_id = attrd.get("ID") or attrd.get("transcript_id") or attrd.get("Name")
                if not tx_id:
                    continue
                transcripts[tx_id] = {
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                }

            elif feature == "CDS":
                parent = attrd.get("Parent") or attrd.get("transcript_id")
                if not parent:
                    continue
                for tx_id in parent.split(","):
                    tx_id = tx_id.strip()
                    if not tx_id:
                        continue
                    cds_by_transcript[tx_id].append({
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "phase": phase,
                    })

    return transcripts, cds_by_transcript


def choose_best_transcript(transcripts, cds_by_transcript):
    candidates = []
    for tx_id, parts in cds_by_transcript.items():
        if not parts:
            continue
        cds_len = sum(abs(p["end"] - p["start"]) + 1 for p in parts)
        candidates.append((tx_id, cds_len))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[1], reverse=True)
    return candidates[0][0]


def extract_spliced_cds_from_record(record_seq, cds_parts):
    """
    Extract CDS from a single-record locus sequence using Helixer-predicted coordinates.
    CDS pieces are concatenated in transcript order.
    """
    if not cds_parts:
        return None, None, "no_cds_parts"

    strands = {p["strand"] for p in cds_parts}
    seqids = {p["seqid"] for p in cds_parts}

    if len(strands) != 1:
        return None, None, "mixed_strands"
    if len(seqids) != 1:
        return None, None, "mixed_seqids"

    strand = list(strands)[0]

    if strand == "+":
        ordered = sorted(cds_parts, key=lambda x: x["start"])
    else:
        ordered = sorted(cds_parts, key=lambda x: x["start"], reverse=True)

    pieces = []
    seq_len = len(record_seq)

    for p in ordered:
        start = min(p["start"], p["end"])
        end = max(p["start"], p["end"])

        if start < 1 or end > seq_len:
            return None, strand, f"cds_out_of_bounds:{start}-{end}/{seq_len}"

        piece = str(record_seq[start - 1:end]).upper()
        pieces.append(piece)

    cds = "".join(pieces)
    if strand == "-":
        cds = str(Seq(cds).reverse_complement())

    if not cds:
        return None, strand, "empty_cds"

    return cds, strand, "ok"


def qc_cds(seq):
    seq = str(seq).upper()
    cds_len = len(seq)
    starts_with_atg = seq.startswith("ATG")
    len_mod3 = cds_len % 3

    usable = cds_len - len_mod3
    trimmed = seq[:usable]
    protein = str(Seq(trimmed).translate(to_stop=False)) if usable > 0 else ""

    if protein.endswith("*"):
        internal_stops = protein[:-1].count("*")
        has_terminal_stop = True
    else:
        internal_stops = protein.count("*")
        has_terminal_stop = False

    return {
        "cds_len": cds_len,
        "len_mod3": len_mod3,
        "starts_with_ATG": starts_with_atg,
        "internal_stop_count": internal_stops,
        "has_terminal_stop": has_terminal_stop,
        "protein_len": len(protein),
        "protein_preview": protein[:60],
    }


def process_family_fasta(family_fa: Path, out_dir: Path, helixer_bin: str, lineage: str, species: str, keep_temp: bool):
    out_records = []
    qc_rows = []
    fail_rows = []

    records = list(SeqIO.parse(str(family_fa), "fasta"))

    for rec in tqdm(records, desc=family_fa.name, unit="gene", leave=False):
        header = rec.description if rec.description else rec.id
        gene_id = rec.description.strip().lstrip(">").split()[0] if rec.description else rec.id

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            rec_fa = td_path / "input.fa"
            out_gff = td_path / "helixer.gff3"

            SeqIO.write([rec], str(rec_fa), "fasta")

            try:
                run_helixer(rec_fa, out_gff, helixer_bin, lineage, species)
                transcripts, cds_by_transcript = parse_helixer_gff(out_gff)
                tx_id = choose_best_transcript(transcripts, cds_by_transcript)

                if tx_id is None:
                    fail_rows.append({
                        "family_fasta": family_fa.name,
                        "gene_id": gene_id,
                        "header_raw": header,
                        "status": "no_transcript_with_cds",
                    })
                    continue

                cds_seq, strand, status = extract_spliced_cds_from_record(
                    rec.seq,
                    cds_by_transcript[tx_id]
                )

                if cds_seq is None:
                    fail_rows.append({
                        "family_fasta": family_fa.name,
                        "gene_id": gene_id,
                        "header_raw": header,
                        "transcript_id": tx_id,
                        "status": status,
                    })
                    continue

                qc = qc_cds(cds_seq)
                qc_rows.append({
                    "family_fasta": family_fa.name,
                    "gene_id": gene_id,
                    "header_raw": header,
                    "transcript_id": tx_id,
                    "strand": strand,
                    "status": "ok",
                    **qc,
                })

                out_records.append(
                    SeqRecord(
                        Seq(cds_seq),
                        id=gene_id,
                        description=f"transcript={tx_id} strand={strand}"
                    )
                )

                if keep_temp:
                    dbg_dir = out_dir / "_helixer_debug" / family_fa.stem / gene_id
                    dbg_dir.mkdir(parents=True, exist_ok=True)
                    SeqIO.write([rec], str(dbg_dir / "input.fa"), "fasta")
                    if out_gff.exists():
                        dbg_dir.joinpath("helixer.gff3").write_text(out_gff.read_text())

            except Exception as e:
                fail_rows.append({
                    "family_fasta": family_fa.name,
                    "gene_id": gene_id,
                    "header_raw": header,
                    "status": f"helixer_failed: {str(e)[:300]}",
                })

    out_fa = out_dir / family_fa.name
    if out_records:
        SeqIO.write(out_records, str(out_fa), "fasta")

    return {
        "family_fasta": family_fa.name,
        "n_input": len(records),
        "cds_written": len(out_records),
        "failed": len(records) - len(out_records),
    }, qc_rows, fail_rows


def main():
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.in_fasta:
        fasta_files = [Path(args.in_fasta)]
    else:
        fasta_files = sorted(Path(args.in_dir).glob(f"*{args.suffix}"))

    if not fasta_files:
        raise SystemExit("No input FASTA files found.")

    summary_rows = []
    all_qc = []
    all_fail = []

    for fam_fa in tqdm(fasta_files, desc="Families", unit="family"):
        summary, qc_rows, fail_rows = process_family_fasta(
            family_fa=fam_fa,
            out_dir=out_dir,
            helixer_bin=args.helixer_bin,
            lineage=args.lineage,
            species=args.species,
            keep_temp=args.keep_temp,
        )
        summary_rows.append(summary)
        all_qc.extend(qc_rows)
        all_fail.extend(fail_rows)

    pd.DataFrame(summary_rows).to_csv(out_dir / "cds_extraction_summary.tsv", sep="\t", index=False)
    pd.DataFrame(all_qc).to_csv(out_dir / "cds_qc.tsv", sep="\t", index=False)

    if all_fail:
        pd.DataFrame(all_fail).to_csv(out_dir / "failed_records.tsv", sep="\t", index=False)

    print(f"[OK] Wrote CDS family FASTAs to: {out_dir}")
    print(f"[OK] Wrote summary: {out_dir / 'cds_extraction_summary.tsv'}")
    print(f"[OK] Wrote QC table: {out_dir / 'cds_qc.tsv'}")
    if all_fail:
        print(f"[WARN] Failed records: {out_dir / 'failed_records.tsv'}")


if __name__ == "__main__":
    main()