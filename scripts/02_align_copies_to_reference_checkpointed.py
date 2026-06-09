#!/usr/bin/env python3

import argparse
import warnings
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import BiopythonDeprecationWarning

warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


VALID_MATCH_STATUSES = {"strict_match", "relaxed_match_no_copy", "matched"}


def load_fasta_as_dict(fasta_path):
    return {rec.id: rec for rec in SeqIO.parse(str(fasta_path), "fasta")}



def compute_alignment_stats(aln_ref, aln_qry):
    """
    Given aligned strings from pairwise2, compute:
    - ref start/end
    - qry start/end
    - matches / mismatches / gaps
    - aligned lengths
    """
    ref_pos = 0
    qry_pos = 0

    ref_start = None
    ref_end = None
    qry_start = None
    qry_end = None

    matches = 0
    mismatches = 0
    gap_cols = 0
    ref_aln_bases = 0
    qry_aln_bases = 0

    for a, b in zip(aln_ref, aln_qry):
        ref_has = a != "-"
        qry_has = b != "-"

        if ref_has:
            ref_pos += 1
        if qry_has:
            qry_pos += 1

        if ref_has or qry_has:
            if ref_start is None and ref_has:
                ref_start = ref_pos
            if qry_start is None and qry_has:
                qry_start = qry_pos

        if ref_has and qry_has:
            ref_end = ref_pos
            qry_end = qry_pos
            ref_aln_bases += 1
            qry_aln_bases += 1
            if a.upper() == b.upper():
                matches += 1
            else:
                mismatches += 1
        elif ref_has or qry_has:
            gap_cols += 1
            if ref_has:
                ref_end = ref_pos
                ref_aln_bases += 1
            if qry_has:
                qry_end = qry_pos
                qry_aln_bases += 1

    aligned_columns = matches + mismatches + gap_cols
    pid = (matches / (matches + mismatches)) * 100 if (matches + mismatches) > 0 else 0.0

    return {
        "ref_aln_start_1based": ref_start,
        "ref_aln_end_1based": ref_end,
        "qry_aln_start_1based": qry_start,
        "qry_aln_end_1based": qry_end,
        "matches": matches,
        "mismatches": mismatches,
        "gap_columns": gap_cols,
        "ref_aln_bases": ref_aln_bases,
        "qry_aln_bases": qry_aln_bases,
        "aligned_columns": aligned_columns,
        "percent_identity_no_gaps": round(pid, 3),
    }



def best_local_alignment(ref_seq, qry_seq):
    """
    Align both forward and reverse complement and keep the better score.
    """
    ref_seq = str(ref_seq).upper()
    qry_fwd = str(qry_seq).upper()
    qry_rev = str(Seq(qry_fwd).reverse_complement())

    alignments_fwd = pairwise2.align.localms(
        ref_seq, qry_fwd, 2, -1, -5, -1, one_alignment_only=True
    )
    alignments_rev = pairwise2.align.localms(
        ref_seq, qry_rev, 2, -1, -5, -1, one_alignment_only=True
    )

    best = None
    strand = None
    qry_used = None

    if alignments_fwd and alignments_rev:
        if alignments_fwd[0].score >= alignments_rev[0].score:
            best = alignments_fwd[0]
            strand = "+"
            qry_used = qry_fwd
        else:
            best = alignments_rev[0]
            strand = "-"
            qry_used = qry_rev
    elif alignments_fwd:
        best = alignments_fwd[0]
        strand = "+"
        qry_used = qry_fwd
    elif alignments_rev:
        best = alignments_rev[0]
        strand = "-"
        qry_used = qry_rev
    else:
        return None

    aln_ref = best.seqA
    aln_qry = best.seqB

    stats = compute_alignment_stats(aln_ref, aln_qry)
    stats["strand"] = strand
    stats["alignment_score"] = best.score
    stats["ref_len"] = len(ref_seq)
    stats["qry_len"] = len(qry_used)
    stats["ref_coverage"] = round((stats["ref_aln_bases"] / len(ref_seq)) * 100, 3) if len(ref_seq) > 0 else 0.0
    stats["qry_coverage"] = round((stats["qry_aln_bases"] / len(qry_used)) * 100, 3) if len(qry_used) > 0 else 0.0
    stats["aligned_ref_string"] = aln_ref
    stats["aligned_qry_string"] = aln_qry
    stats["oriented_qry_seq"] = qry_used

    return stats



def pick_row_value(row, candidates):
    for c in candidates:
        if c in row.index and pd.notna(row[c]):
            return row[c]
    return None



def append_tsv_rows(rows, out_tsv, include_alignment_strings=False, wrote_header=False):
    if not rows:
        return wrote_header

    df = pd.DataFrame(rows)
    if not include_alignment_strings:
        for c in ["aligned_ref_string", "aligned_qry_string"]:
            if c in df.columns:
                df.drop(columns=c, inplace=True)

    mode = "a" if wrote_header else "w"
    df.to_csv(out_tsv, sep="\t", index=False, mode=mode, header=not wrote_header)
    return True



def append_alignment_records(records_by_gene, combined_fa=None, gene_dir=None):
    total = 0

    if combined_fa is not None:
        with open(combined_fa, "a") as handle:
            for gene, records in records_by_gene.items():
                if records:
                    SeqIO.write(records, handle, "fasta")
                    total += len(records)

    if gene_dir is not None:
        gene_dir.mkdir(parents=True, exist_ok=True)
        for gene, records in records_by_gene.items():
            if not records:
                continue
            gene_fa = gene_dir / f"{gene}_oriented.fasta"
            with open(gene_fa, "a") as handle:
                SeqIO.write(records, handle, "fasta")

    return total



def main():
    ap = argparse.ArgumentParser(
        description="Reference-guided alignment of extracted loci to reference genes with incremental TSV and FASTA writing."
    )
    ap.add_argument("--master", required=True, help="copy_locus_master.tsv from script 1")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument(
        "--write-oriented-fasta",
        action="store_true",
        help="Write reoriented copy sequences as a combined FASTA and per-gene FASTAs"
    )
    ap.add_argument(
        "--write-alignment-strings",
        action="store_true",
        help="Write aligned strings into output TSV"
    )
    ap.add_argument(
        "--batch-size",
        type=int,
        default=200,
        help="Flush TSV and FASTA output every N processed rows"
    )
    ap.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="Print progress every N processed rows"
    )
    args = ap.parse_args()

    master = pd.read_csv(args.master, sep="\t")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    required = {"family_id", "gene", "copy_header", "copy_fasta", "ref_fasta", "match_status"}
    missing = required - set(master.columns)
    if missing:
        raise ValueError(f"Master table missing required columns: {missing}")

    work = master[master["match_status"].isin(VALID_MATCH_STATUSES)].copy()
    work = work.dropna(subset=["copy_header", "copy_fasta", "ref_fasta", "gene"])

    if len(work) == 0:
        raise ValueError("No rows with successful match_status values were found in the master table.")

    total_rows = len(work)
    out_tsv = out_dir / "copy_reference_alignment_metrics.tsv"
    combined_fa = out_dir / "copy_sequences_best_orientation.fasta"
    gene_dir = out_dir / "oriented_fastas_by_gene"
    progress_file = out_dir / "run_progress.tsv"
    summary_file = out_dir / "alignment_status_summary.tsv"

    # fresh output files for a new run
    if out_tsv.exists():
        out_tsv.unlink()
    if progress_file.exists():
        progress_file.unlink()
    if summary_file.exists():
        summary_file.unlink()
    if combined_fa.exists() and args.write_oriented_fasta:
        combined_fa.unlink()
    if gene_dir.exists() and args.write_oriented_fasta:
        for p in gene_dir.glob("*.fasta"):
            p.unlink()

    copy_cache = {}
    ref_cache = {}

    batch_rows = []
    batch_records_by_gene = defaultdict(list)
    wrote_header = False
    processed = 0
    oriented_written = 0
    status_counts = defaultdict(int)

    with open(progress_file, "w") as ph:
        ph.write("processed\tstatus_ok\tstatus_other\toriented_written\n")

    print(f"[INFO] Rows to process: {total_rows}")
    print(f"[INFO] Writing metrics to: {out_tsv}")
    if args.write_oriented_fasta:
        print(f"[INFO] Writing combined oriented FASTA to: {combined_fa}")
        print(f"[INFO] Writing per-gene oriented FASTAs to: {gene_dir}")

    for _, row in work.iterrows():
        copy_fasta = row["copy_fasta"]
        ref_fasta = row["ref_fasta"]
        copy_header = row["copy_header"]
        gene = row["gene"]

        assembly_val = pick_row_value(row, ["assembly", "assembly_full", "assembly_short"])
        assembly_short_val = pick_row_value(row, ["assembly_short", "assembly"])
        scaffold_val = pick_row_value(row, ["scaffold", "scaffold_full", "scaffold_key", "scaffold_raw"])
        copy_val = pick_row_value(row, ["copy", "copy_tt", "copy_extract_idx", "copy_header_idx"])

        if copy_fasta not in copy_cache:
            copy_cache[copy_fasta] = load_fasta_as_dict(copy_fasta)
        if ref_fasta not in ref_cache:
            ref_cache[ref_fasta] = load_fasta_as_dict(ref_fasta)

        copy_dict = copy_cache[copy_fasta]
        ref_dict = ref_cache[ref_fasta]

        if copy_header not in copy_dict:
            batch_rows.append({
                "family_id": row["family_id"],
                "assembly": assembly_val,
                "assembly_short": assembly_short_val,
                "gene": gene,
                "copy": copy_val,
                "copy_header": copy_header,
                "match_status": row["match_status"],
                "status": "copy_header_not_found_in_fasta"
            })
            status_counts["copy_header_not_found_in_fasta"] += 1
            processed += 1
        elif gene not in ref_dict:
            batch_rows.append({
                "family_id": row["family_id"],
                "assembly": assembly_val,
                "assembly_short": assembly_short_val,
                "gene": gene,
                "copy": copy_val,
                "copy_header": copy_header,
                "match_status": row["match_status"],
                "status": "reference_gene_not_found_in_ref_fasta"
            })
            status_counts["reference_gene_not_found_in_ref_fasta"] += 1
            processed += 1
        else:
            copy_rec = copy_dict[copy_header]
            ref_rec = ref_dict[gene]

            aln = best_local_alignment(ref_rec.seq, copy_rec.seq)
            if aln is None:
                batch_rows.append({
                    "family_id": row["family_id"],
                    "assembly": assembly_val,
                    "assembly_short": assembly_short_val,
                    "gene": gene,
                    "copy": copy_val,
                    "copy_header": copy_header,
                    "match_status": row["match_status"],
                    "status": "alignment_failed"
                })
                status_counts["alignment_failed"] += 1
                processed += 1
            else:
                out_row = {
                    "family_id": row["family_id"],
                    "assembly": assembly_val,
                    "assembly_short": assembly_short_val,
                    "scaffold": scaffold_val,
                    "gene": gene,
                    "copy": copy_val,
                    "copy_header": copy_header,
                    "match_status": row["match_status"],
                    "status": "ok",
                    "strand": aln["strand"],
                    "alignment_score": aln["alignment_score"],
                    "ref_len": aln["ref_len"],
                    "qry_len": aln["qry_len"],
                    "ref_aln_start_1based": aln["ref_aln_start_1based"],
                    "ref_aln_end_1based": aln["ref_aln_end_1based"],
                    "qry_aln_start_1based": aln["qry_aln_start_1based"],
                    "qry_aln_end_1based": aln["qry_aln_end_1based"],
                    "matches": aln["matches"],
                    "mismatches": aln["mismatches"],
                    "gap_columns": aln["gap_columns"],
                    "ref_aln_bases": aln["ref_aln_bases"],
                    "qry_aln_bases": aln["qry_aln_bases"],
                    "aligned_columns": aln["aligned_columns"],
                    "percent_identity_no_gaps": aln["percent_identity_no_gaps"],
                    "ref_coverage": aln["ref_coverage"],
                    "qry_coverage": aln["qry_coverage"],
                }

                if args.write_alignment_strings:
                    out_row["aligned_ref_string"] = aln["aligned_ref_string"]
                    out_row["aligned_qry_string"] = aln["aligned_qry_string"]

                batch_rows.append(out_row)
                status_counts["ok"] += 1
                processed += 1

                if args.write_oriented_fasta:
                    batch_records_by_gene[gene].append(
                        SeqRecord(
                            Seq(aln["oriented_qry_seq"]),
                            id=copy_header,
                            description=f"strand={aln['strand']} gene={gene} family={row['family_id']}"
                        )
                    )

        if processed % args.batch_size == 0 or processed == total_rows:
            wrote_header = append_tsv_rows(
                batch_rows,
                out_tsv,
                include_alignment_strings=args.write_alignment_strings,
                wrote_header=wrote_header,
            )
            batch_rows = []

            if args.write_oriented_fasta:
                oriented_written += append_alignment_records(
                    batch_records_by_gene,
                    combined_fa=combined_fa,
                    gene_dir=gene_dir,
                )
                batch_records_by_gene = defaultdict(list)

            with open(progress_file, "a") as ph:
                non_ok = sum(v for k, v in status_counts.items() if k != "ok")
                ph.write(f"{processed}\t{status_counts['ok']}\t{non_ok}\t{oriented_written}\n")

        if processed % args.progress_every == 0 or processed == total_rows:
            non_ok = sum(v for k, v in status_counts.items() if k != "ok")
            print(
                f"[PROGRESS] {processed}/{total_rows} processed | ok={status_counts['ok']} | "
                f"other={non_ok} | oriented_written={oriented_written}"
            )

    summary_rows = [{"status": k, "n": status_counts[k]} for k in sorted(status_counts.keys())]
    pd.DataFrame(summary_rows).to_csv(summary_file, sep="\t", index=False)

    print(f"[OK] Alignment metrics written to: {out_tsv}")
    if args.write_oriented_fasta:
        print(f"[OK] Combined oriented FASTA written to: {combined_fa}")
        print(f"[OK] Per-gene oriented FASTAs written to: {gene_dir}")
    print(f"[OK] Progress log written to: {progress_file}")
    print(f"[OK] Status summary written to: {summary_file}")
    print(f"[OK] Total rows processed: {processed}")


if __name__ == "__main__":
    main()
