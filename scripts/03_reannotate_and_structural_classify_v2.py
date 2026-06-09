#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


# -----------------------------
# basic FASTA helpers
# -----------------------------

def load_fasta_as_dict(path):
    return {rec.id: rec for rec in SeqIO.parse(str(path), "fasta")}


def translate_cds_safe(seq):
    seq = str(seq).upper().replace("U", "T")
    usable = len(seq) - (len(seq) % 3)
    seq = seq[:usable]
    if usable == 0:
        return ""
    return str(Seq(seq).translate(to_stop=False))


def count_internal_stops(protein):
    if not protein:
        return np.nan
    if protein.endswith("*"):
        protein_core = protein[:-1]
    else:
        protein_core = protein
    return protein_core.count("*")


def has_terminal_stop(protein):
    return protein.endswith("*") if protein else False


def starts_with_methionine(protein):
    return protein.startswith("M") if protein else False


# -----------------------------
# exonerate running/parsing
# -----------------------------

def run_exonerate(query_protein_fa, target_dna_fa):
    cmd = [
        "exonerate",
        "--model", "protein2genome",
        "--bestn", "1",
        "--showalignment", "no",
        "--showvulgar", "yes",
        "--showtargetgff", "yes",
        str(query_protein_fa),
        str(target_dna_fa),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Exonerate failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
    return res.stdout


def parse_exonerate_vulgar(text):
    """
    Parse first vulgar line.
    Example structure:
    vulgar: query qstart qend qstrand target tstart tend tstrand score ops...
    """
    vulgar_lines = [x for x in text.splitlines() if x.startswith("vulgar:")]
    if not vulgar_lines:
        return None

    parts = vulgar_lines[0].strip().split()
    if len(parts) < 10:
        return None

    out = {
        "query_id": parts[1],
        "qstart": int(parts[2]),
        "qend": int(parts[3]),
        "qstrand": parts[4],
        "target_id": parts[5],
        "tstart": int(parts[6]),
        "tend": int(parts[7]),
        "tstrand": parts[8],
        "score": float(parts[9]),
        "raw_ops": parts[10:]
    }

    ops = parts[10:]
    op_counts = {}
    for i in range(0, len(ops), 3):
        if i + 2 < len(ops):
            op = ops[i]
            qlen = int(ops[i + 1])
            tlen = int(ops[i + 2])
            op_counts[op] = op_counts.get(op, 0) + 1
    out["op_counts"] = op_counts
    out["frameshift_ops"] = op_counts.get("F", 0)
    return out


def parse_exonerate_target_gff(text):
    """
    Extract CDS coordinates from exonerate target GFF block.
    Return list of tuples: (seqid, start, end, strand)
    """
    cds_rows = []

    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue

        seqid, source, feature, start, end, score, strand, phase, attrs = parts
        feature_l = feature.lower()

        if feature_l == "cds":
            cds_rows.append((seqid, int(start), int(end), strand))

    return cds_rows


def extract_spliced_cds(target_record, cds_coords):
    """
    cds_coords: list of (seqid, start, end, strand), 1-based inclusive
    """
    if not cds_coords:
        return None, None, None

    strands = set([x[3] for x in cds_coords])
    if len(strands) != 1:
        return None, None, None

    strand = list(strands)[0]

    # sort in transcript order
    if strand == "+":
        cds_coords_sorted = sorted(cds_coords, key=lambda x: x[1])
    else:
        cds_coords_sorted = sorted(cds_coords, key=lambda x: x[1], reverse=True)

    pieces = []
    for _, start, end, _ in cds_coords_sorted:
        s = min(start, end) - 1
        e = max(start, end)
        pieces.append(str(target_record.seq[s:e]))

    cds = "".join(pieces).upper()
    if strand == "-":
        cds = str(Seq(cds).reverse_complement())

    exon_count = len(cds_coords_sorted)
    return cds, strand, exon_count


# -----------------------------
# alignment stats
# -----------------------------

def protein_global_stats(ref_prot, qry_prot):
    if not ref_prot or not qry_prot:
        return {
            "prot_matches": np.nan,
            "prot_mismatches": np.nan,
            "prot_gap_columns": np.nan,
            "prot_percent_identity": np.nan,
            "prot_ref_coverage": np.nan,
            "prot_qry_coverage": np.nan,
            "prot_aligned_columns": np.nan,
        }

    aln = pairwise2.align.globalms(ref_prot, qry_prot, 2, -1, -5, -1, one_alignment_only=True)
    if not aln:
        return {
            "prot_matches": np.nan,
            "prot_mismatches": np.nan,
            "prot_gap_columns": np.nan,
            "prot_percent_identity": np.nan,
            "prot_ref_coverage": np.nan,
            "prot_qry_coverage": np.nan,
            "prot_aligned_columns": np.nan,
        }

    a = aln[0].seqA
    b = aln[0].seqB

    matches = 0
    mismatches = 0
    gap_cols = 0
    ref_non_gap = 0
    qry_non_gap = 0

    for x, y in zip(a, b):
        if x != "-":
            ref_non_gap += 1
        if y != "-":
            qry_non_gap += 1

        if x == "-" or y == "-":
            gap_cols += 1
        elif x == y:
            matches += 1
        else:
            mismatches += 1

    denom = matches + mismatches
    pid = (matches / denom) * 100 if denom > 0 else np.nan

    return {
        "prot_matches": matches,
        "prot_mismatches": mismatches,
        "prot_gap_columns": gap_cols,
        "prot_percent_identity": round(pid, 3) if not pd.isna(pid) else np.nan,
        "prot_ref_coverage": round(ref_non_gap / len(ref_prot) * 100, 3) if len(ref_prot) > 0 else np.nan,
        "prot_qry_coverage": round(qry_non_gap / len(qry_prot) * 100, 3) if len(qry_prot) > 0 else np.nan,
        "prot_aligned_columns": len(a),
    }


def count_nontriplet_gap_runs(ref_cds, qry_cds):
    """
    Align ref CDS vs predicted CDS and count contiguous gap runs with length % 3 != 0.
    This is a 'frameshift-like indel' heuristic.
    """
    if not ref_cds or not qry_cds:
        return np.nan

    aln = pairwise2.align.globalms(ref_cds, qry_cds, 2, -1, -5, -1, one_alignment_only=True)
    if not aln:
        return np.nan

    a = aln[0].seqA
    b = aln[0].seqB

    nontriplet_runs = 0
    i = 0
    while i < len(a):
        if a[i] == "-" or b[i] == "-":
            run_len = 0
            while i < len(a) and (a[i] == "-" or b[i] == "-"):
                run_len += 1
                i += 1
            if run_len % 3 != 0:
                nontriplet_runs += 1
        else:
            i += 1

    return nontriplet_runs


# -----------------------------
# classification
# -----------------------------

def structural_classify(
    predicted_cds_len,
    ref_cds_len,
    predicted_prot,
    ref_prot,
    internal_stop_count,
    nontriplet_gap_runs,
    prot_ref_coverage,
    exonerate_score,
    exon_count
):
    """
    Returns P, U0, or A
    U0 = structurally intact candidate, to be split into U/F in the next script
    """
    if pd.isna(predicted_cds_len) or predicted_cds_len == 0:
        return "A", "no_predicted_cds"

    if pd.isna(exonerate_score):
        return "A", "no_exonerate_model"

    if not predicted_prot:
        return "A", "empty_predicted_protein"

    ref_prot_len = len(ref_prot) if ref_prot else np.nan
    pred_prot_len = len(predicted_prot)
    prot_len_ratio = pred_prot_len / ref_prot_len if ref_prot_len and ref_prot_len > 0 else np.nan

    if internal_stop_count > 0:
        return "P", "internal_stop"

    if not pd.isna(nontriplet_gap_runs) and nontriplet_gap_runs > 0:
        return "P", "frameshift_like_indel"

    if not pd.isna(prot_len_ratio) and prot_len_ratio < 0.75:
        return "P", "severely_truncated"

    if not pd.isna(prot_ref_coverage) and prot_ref_coverage < 75:
        return "A", "low_ref_coverage"

    if exon_count is None or exon_count == 0:
        return "A", "no_exons"

    return "U0", "structurally_intact"


# -----------------------------
# main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Re-annotate copies with exonerate and classify P/U0/A.")
    ap.add_argument("--master", required=True, help="copy_locus_master.tsv")
    ap.add_argument("--step2", required=True, help="copy_reference_alignment_metrics.tsv")
    ap.add_argument("--out-dir", required=True, help="output directory")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(args.master, sep="\t")
    step2 = pd.read_csv(args.step2, sep="\t")

    key_cols = ["family_id", "assembly_short", "gene", "copy", "copy_header"]
    merged = master.merge(
        step2[[c for c in step2.columns if c in set(key_cols + ["status", "strand", "alignment_score"])]],
        on=key_cols,
        how="left",
        suffixes=("", "_step2")
    )

    work = merged[(merged["match_status"] == "matched")].copy()

    copy_cache = {}
    ref_cache = {}

    cds_records = []
    prot_records = []
    rows = []

    missing_copy_fasta_count = 0
    missing_ref_fasta_count = 0

    for _, row in tqdm(work.iterrows(), total=len(work), desc="Processing copies", unit="copy"):
        family_id = row["family_id"]
        gene = row["gene"]
        copy_header = row["copy_header"]
        copy_fasta = row["copy_fasta"]
        ref_fasta = row["ref_fasta"]

        # Safe FASTA loading
        if copy_fasta not in copy_cache:
            try:
                if pd.isna(copy_fasta) or not os.path.exists(copy_fasta):
                    copy_cache[copy_fasta] = None
                    print(f"[WARN] Missing copy FASTA: {copy_fasta}")
                else:
                    copy_cache[copy_fasta] = load_fasta_as_dict(copy_fasta)
            except Exception as e:
                print(f"[WARN] Failed loading copy FASTA: {copy_fasta} | {str(e)}")
                copy_cache[copy_fasta] = None

        if ref_fasta not in ref_cache:
            try:
                if pd.isna(ref_fasta) or not os.path.exists(ref_fasta):
                    ref_cache[ref_fasta] = None
                    print(f"[WARN] Missing reference FASTA: {ref_fasta}")
                else:
                    ref_cache[ref_fasta] = load_fasta_as_dict(ref_fasta)
            except Exception as e:
                print(f"[WARN] Failed loading reference FASTA: {ref_fasta} | {str(e)}")
                ref_cache[ref_fasta] = None

        base = {
            "family_id": family_id,
            "assembly": row.get("assembly", np.nan),
            "assembly_short": row.get("assembly_short", np.nan),
            "scaffold": row.get("scaffold", np.nan),
            "gene": gene,
            "copy": row.get("copy", np.nan),
            "copy_header": copy_header,
            "copy_fasta": row.get("copy_fasta", np.nan),
            "ref_fasta": row.get("ref_fasta", np.nan),
        }

        copy_dict = copy_cache.get(copy_fasta)
        ref_dict = ref_cache.get(ref_fasta)

        if copy_dict is None:
            missing_copy_fasta_count += 1
            rows.append({**base, "structural_status": "A", "structural_reason": "missing_copy_fasta"})
            continue

        if ref_dict is None:
            missing_ref_fasta_count += 1
            rows.append({**base, "structural_status": "A", "structural_reason": "missing_reference_fasta"})
            continue

        if copy_header not in copy_dict:
            rows.append({**base, "structural_status": "A", "structural_reason": "copy_missing_in_fasta"})
            continue

        if gene not in ref_dict:
            rows.append({**base, "structural_status": "A", "structural_reason": "reference_gene_missing"})
            continue

        target_rec = copy_dict[copy_header]
        ref_rec = ref_dict[gene]

        ref_cds = str(ref_rec.seq).upper()
        ref_prot = translate_cds_safe(ref_cds)

        with tempfile.TemporaryDirectory() as td:
            td = Path(td)

            query_prot_fa = td / "query_protein.fa"
            target_dna_fa = td / "target_dna.fa"

            SeqIO.write([SeqRecord(Seq(ref_prot), id=gene, description="")], str(query_prot_fa), "fasta")
            SeqIO.write([SeqRecord(target_rec.seq, id=copy_header, description="")], str(target_dna_fa), "fasta")

            try:
                exo_text = run_exonerate(query_prot_fa, target_dna_fa)
            except Exception as e:
                rows.append({**base, "structural_status": "A", "structural_reason": f"exonerate_failed: {str(e)[:180]}"})
                continue

        vulgar = parse_exonerate_vulgar(exo_text)
        cds_coords = parse_exonerate_target_gff(exo_text)

        predicted_cds, exo_strand, exon_count = extract_spliced_cds(target_rec, cds_coords)

        if predicted_cds is None:
            rows.append({
                **base,
                "structural_status": "A",
                "structural_reason": "no_cds_recovered",
                "predicted_cds_len": np.nan,
                "predicted_protein_len": np.nan,
                "ref_cds_len": len(ref_cds),
                "ref_protein_len": len(ref_prot),
            })
            continue

        predicted_prot = translate_cds_safe(predicted_cds)
        internal_stop_count = count_internal_stops(predicted_prot)
        terminal_stop = has_terminal_stop(predicted_prot)
        start_m = starts_with_methionine(predicted_prot)
        nontriplet_gap_runs = count_nontriplet_gap_runs(ref_cds, predicted_cds)

        prot_stats = protein_global_stats(ref_prot, predicted_prot)

        structural_status, structural_reason = structural_classify(
            predicted_cds_len=len(predicted_cds),
            ref_cds_len=len(ref_cds),
            predicted_prot=predicted_prot,
            ref_prot=ref_prot,
            internal_stop_count=internal_stop_count,
            nontriplet_gap_runs=nontriplet_gap_runs,
            prot_ref_coverage=prot_stats["prot_ref_coverage"],
            exonerate_score=vulgar["score"] if vulgar else np.nan,
            exon_count=exon_count,
        )

        rows.append({
            **base,
            "exo_score": vulgar["score"] if vulgar else np.nan,
            "exo_frameshift_ops": vulgar["frameshift_ops"] if vulgar else np.nan,
            "exo_strand": exo_strand,
            "exon_count": exon_count,
            "predicted_cds_len": len(predicted_cds),
            "predicted_protein_len": len(predicted_prot),
            "ref_cds_len": len(ref_cds),
            "ref_protein_len": len(ref_prot),
            "protein_len_ratio": round(len(predicted_prot) / len(ref_prot), 4) if len(ref_prot) > 0 else np.nan,
            "internal_stop_count": internal_stop_count,
            "has_terminal_stop": terminal_stop,
            "starts_with_M": start_m,
            "frameshift_like_gap_runs": nontriplet_gap_runs,
            **prot_stats,
            "structural_status": structural_status,
            "structural_reason": structural_reason,
        })

        cds_records.append(
            SeqRecord(Seq(predicted_cds), id=copy_header, description=f"family={family_id} gene={gene} structural={structural_status}")
        )
        prot_records.append(
            SeqRecord(Seq(predicted_prot), id=copy_header, description=f"family={family_id} gene={gene} structural={structural_status}")
        )

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_dir / "step3_structural_classification.tsv", sep="\t", index=False)

    if cds_records:
        SeqIO.write(cds_records, str(out_dir / "step3_predicted_cds.fasta"), "fasta")
    if prot_records:
        SeqIO.write(prot_records, str(out_dir / "step3_predicted_proteins.fasta"), "fasta")

    print(f"[OK] Wrote: {out_dir / 'step3_structural_classification.tsv'}")
    print(f"[OK] Wrote: {out_dir / 'step3_predicted_cds.fasta'}")
    print(f"[OK] Wrote: {out_dir / 'step3_predicted_proteins.fasta'}")
    print(f"[INFO] Missing copy FASTA rows: {missing_copy_fasta_count}")
    print(f"[INFO] Missing reference FASTA rows: {missing_ref_fasta_count}")


if __name__ == "__main__":
    main()