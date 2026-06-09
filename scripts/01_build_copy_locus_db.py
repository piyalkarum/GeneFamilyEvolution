#!/usr/bin/env python3

import re
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO


def infer_family_id(path: Path) -> str:
    return path.name[:9]


def detect_delimiter(file_path: Path) -> str:
    with open(file_path, "r") as f:
        first_line = f.readline()
    if "\t" in first_line:
        return "\t"
    elif "," in first_line:
        return ","
    else:
        return r"\s+"


def extract_assembly_short_from_tt(x: str) -> str:
    """
    Handles:
      GCA_020911765.2_ASM2091176v2_genomic -> ASM2091176v2
      GCA_025388435.1_AtA1-141.20181003_genomic -> AtA1-141.20181003
      GCA_946402385.1_AUZE-A-5.genomic -> AUZE-A-5.
    """
    s = str(x).strip()

    m = re.match(r"^GCA_[^_]+_(.+)$", s)
    if m:
        s = m.group(1)

    if s.endswith("_genomic"):
        s = s[:-len("_genomic")]
    elif s.endswith("genomic"):
        s = s[:-len("genomic")]

    return s


def load_tt_table(tt_file: Path) -> pd.DataFrame:
    sep = detect_delimiter(tt_file)
    if sep == r"\s+":
        df = pd.read_csv(tt_file, sep=sep, engine="python")
    else:
        df = pd.read_csv(tt_file, sep=sep)

    expected_cols = {"assembly", "gene", "scaffold", "start", "end", "length", "max_length", "copy"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"{tt_file} is missing required columns: {missing}")

    df = df.copy()
    df["family_id"] = infer_family_id(tt_file)
    df["tt_file"] = str(tt_file)
    df["assembly_short"] = df["assembly"].apply(extract_assembly_short_from_tt)

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["length"] = pd.to_numeric(df["length"], errors="coerce")
    df["max_length"] = pd.to_numeric(df["max_length"], errors="coerce")
    df["copy"] = pd.to_numeric(df["copy"], errors="coerce").astype("Int64")

    # preserve original order, because the R extraction script used row order
    df["_tt_row_order"] = range(1, len(df) + 1)

    # this is the rank actually used by the R extraction script for naming _copy_n
    group_cols = ["family_id", "assembly_short", "gene", "scaffold"]
    df["copy_rank"] = df.groupby(group_cols).cumcount() + 1

    return df


def parse_reference_fasta(ref_fasta: Path) -> pd.DataFrame:
    rows = []
    family_id = infer_family_id(ref_fasta)

    for rec in SeqIO.parse(str(ref_fasta), "fasta"):
        gene = rec.id.strip()
        rows.append({
            "family_id": family_id,
            "gene": gene,
            "ref_header": rec.id,
            "ref_seq_len": len(rec.seq),
            "ref_fasta": str(ref_fasta),
        })

    return pd.DataFrame(rows)


def parse_extracted_header(header: str):
    """
    Supports both:
      ASM2091176v2__AT2G28700_CP087127.2_copy_1
      AUZE-A-5._AT5G48670_OX291433.1_copy_1

    Returns:
      assembly_short, gene, scaffold, copy_rank
    """
    header = header.strip()

    if "_copy_" not in header:
        return None

    try:
        left, copy_str = header.rsplit("_copy_", 1)
    except ValueError:
        return None

    try:
        copy_rank = int(copy_str)
    except ValueError:
        return None

    if "__" in left:
        assembly_short, rest = left.split("__", 1)
    else:
        if "_" not in left:
            return None
        assembly_short, rest = left.split("_", 1)

    parts = rest.split("_")
    if len(parts) < 2:
        return None

    gene = parts[0]
    scaffold = "_".join(parts[1:])

    return assembly_short, gene, scaffold, copy_rank


def parse_extracted_fasta(extracted_fasta: Path) -> pd.DataFrame:
    rows = []
    family_id = infer_family_id(extracted_fasta)

    for rec in SeqIO.parse(str(extracted_fasta), "fasta"):
        parsed = parse_extracted_header(rec.id)
        if parsed is None:
            rows.append({
                "family_id": family_id,
                "copy_header": rec.id,
                "assembly_short": None,
                "gene": None,
                "scaffold": None,
                "copy_rank": None,
                "copy_seq_len": len(rec.seq),
                "copy_fasta": str(extracted_fasta),
                "parse_ok": False,
            })
            continue

        assembly_short, gene, scaffold, copy_rank = parsed
        rows.append({
            "family_id": family_id,
            "copy_header": rec.id,
            "assembly_short": assembly_short,
            "gene": gene,
            "scaffold": scaffold,
            "copy_rank": copy_rank,
            "copy_seq_len": len(rec.seq),
            "copy_fasta": str(extracted_fasta),
            "parse_ok": True,
        })

    return pd.DataFrame(rows)


def map_files_by_family(file_list):
    d = {}
    for f in file_list:
        fam = infer_family_id(f)
        d.setdefault(fam, []).append(f)
    return d


def choose_single_file_per_family(files_by_family, label):
    out = {}
    for fam, files in files_by_family.items():
        if len(files) == 1:
            out[fam] = files[0]
        else:
            raise ValueError(
                f"Family {fam} has multiple {label} files:\n" +
                "\n".join(str(x) for x in files)
            )
    return out


def main():
    ap = argparse.ArgumentParser(description="Build master per-copy locus database.")
    ap.add_argument("--tt-dir", required=True, help="Directory with family-wise TT tables")
    ap.add_argument("--ref-dir", required=True, help="Directory with family reference FASTA files")
    ap.add_argument("--copy-dir", required=True, help="Directory with extracted-copy FASTA files")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    args = ap.parse_args()

    tt_dir = Path(args.tt_dir)
    ref_dir = Path(args.ref_dir)
    copy_dir = Path(args.copy_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tt_files = sorted([p for p in tt_dir.iterdir() if p.is_file() and not p.name.startswith(".")])
    ref_files = sorted([p for p in ref_dir.iterdir() if p.is_file() and not p.name.startswith(".") and p.suffix.lower() in [".fa", ".fasta", ".fna"]])
    copy_files = sorted([p for p in copy_dir.iterdir() if p.is_file() and not p.name.startswith(".") and p.suffix.lower() in [".fa", ".fasta", ".fna"]])

    tt_by_family = choose_single_file_per_family(map_files_by_family(tt_files), "TT table")
    ref_by_family = choose_single_file_per_family(map_files_by_family(ref_files), "reference FASTA")
    copy_by_family = choose_single_file_per_family(map_files_by_family(copy_files), "copy FASTA")

    all_fams = sorted(set(tt_by_family) & set(ref_by_family) & set(copy_by_family))
    if not all_fams:
        raise ValueError("No overlapping family IDs found across TT / reference / copy files.")

    all_master = []
    all_unmatched_tt = []
    all_unmatched_copy = []
    all_bad_headers = []

    for fam in all_fams:
        tt_df = load_tt_table(tt_by_family[fam])
        ref_df = parse_reference_fasta(ref_by_family[fam])
        copy_df = parse_extracted_fasta(copy_by_family[fam])

        bad = copy_df[copy_df["parse_ok"] == False].copy()
        if len(bad) > 0:
            all_bad_headers.append(bad)

        copy_good = copy_df[copy_df["parse_ok"] == True].copy()

        # Match on copy_rank, not TT copy column
        merged = tt_df.merge(
            copy_good,
            on=["family_id", "assembly_short", "gene", "scaffold", "copy_rank"],
            how="left",
            suffixes=("", "_copy")
        )

        merged = merged.merge(
            ref_df[["family_id", "gene", "ref_header", "ref_seq_len", "ref_fasta"]],
            on=["family_id", "gene"],
            how="left"
        )

        merged["match_status"] = "matched"
        merged.loc[merged["copy_header"].isna(), "match_status"] = "tt_without_copy"

        unmatched_tt = merged[merged["match_status"] == "tt_without_copy"].copy()
        if len(unmatched_tt) > 0:
            all_unmatched_tt.append(unmatched_tt)

        matched_keys = merged[["family_id", "assembly_short", "gene", "scaffold", "copy_rank"]].drop_duplicates()
        copy_only = copy_good.merge(
            matched_keys,
            on=["family_id", "assembly_short", "gene", "scaffold", "copy_rank"],
            how="left",
            indicator=True
        )
        copy_only = copy_only[copy_only["_merge"] == "left_only"].drop(columns=["_merge"])
        if len(copy_only) > 0:
            all_unmatched_copy.append(copy_only)

        all_master.append(merged)

    master_df = pd.concat(all_master, ignore_index=True)

    preferred_cols = [
        "family_id", "assembly", "assembly_short", "gene", "scaffold",
        "start", "end", "length", "max_length", "copy", "copy_rank",
        "copy_header", "copy_seq_len", "copy_fasta",
        "ref_header", "ref_seq_len", "ref_fasta",
        "tt_file", "match_status"
    ]
    cols = [c for c in preferred_cols if c in master_df.columns] + [c for c in master_df.columns if c not in preferred_cols]
    master_df = master_df[cols]

    master_out = out_dir / "copy_locus_master.tsv"
    master_df.to_csv(master_out, sep="\t", index=False)

    if all_unmatched_tt:
        pd.concat(all_unmatched_tt, ignore_index=True).to_csv(out_dir / "unmatched_tt_rows.tsv", sep="\t", index=False)

    if all_unmatched_copy:
        pd.concat(all_unmatched_copy, ignore_index=True).to_csv(out_dir / "unmatched_copy_headers.tsv", sep="\t", index=False)

    if all_bad_headers:
        pd.concat(all_bad_headers, ignore_index=True).to_csv(out_dir / "bad_copy_headers.tsv", sep="\t", index=False)

    print(f"[OK] Master table written to: {master_out}")
    print(f"[OK] Families processed: {len(all_fams)}")
    print(f"[OK] Total TT rows: {len(master_df)}")


if __name__ == "__main__":
    main()
