# coding:utf-8
"""
Build final AlloPipe cohort tables after Nextflow pair fan-out.
"""
import argparse
import re
import shutil
from pathlib import Path

import pandas as pd
from sklearn import linear_model as lm

from table_operations import create_AAMS_df, get_ref_ratio_pair


def pair_sort_key(pair):
    match = re.search(r"\d+", str(pair))
    return int(match.group()) if match else 0


def pair_from_name(path):
    name = Path(path).name
    match = re.match(r"(P\d+)_", name)
    return match.group(1) if match else None


def first_for_pair(paths, pair, marker):
    candidates = [
        path for path in paths
        if marker in Path(path).name and (pair is None or Path(path).name.startswith(f"{pair}_"))
    ]
    if not candidates and len(paths) == 1:
        return paths[0]
    if not candidates:
        raise FileNotFoundError(f"No {marker} table found for pair {pair}")
    return sorted(candidates)[0]


def normalize_ams(ams_pkls, donor_tables, recipient_tables):
    rows = []
    for ams_pkl in sorted(ams_pkls, key=lambda p: pair_sort_key(pair_from_name(p))):
        ams_df = pd.read_pickle(ams_pkl)
        pair = ams_df["pair"].iloc[0]
        donor_table = first_for_pair(donor_tables, pair, "_D0_")
        recipient_table = first_for_pair(recipient_tables, pair, "_R0_")
        donor_df = pd.read_csv(donor_table, sep="\t", dtype={"CHROM": str})
        recipient_df = pd.read_csv(recipient_table, sep="\t", dtype={"CHROM": str})
        common_ref, total_ref, ref_ratio = get_ref_ratio_pair(donor_df, recipient_df)
        ams_index = ams_df.columns.get_loc("ams")
        ams_df.insert(ams_index + 1, "common_ref", common_ref)
        ams_df.insert(ams_index + 2, "total_ref", total_ref)
        ams_df.insert(ams_index + 3, "ref_ratio", ref_ratio)
        rows.append(ams_df)

    ams_df = pd.concat(rows).reset_index(drop=True)
    ams_df = ams_df.sort_values("pair", key=lambda col: col.map(pair_sort_key))

    if ams_df["ref_ratio"].notna().all():
        x = ams_df["ref_ratio"].values.reshape(-1, 1)
        y = ams_df["ams"].values.reshape(-1, 1)
        reg = lm.LinearRegression().fit(x, y)
        ref_ratio_mean = ams_df["ref_ratio"].mean()
        ams_df["ams_norm"] = (
            float(reg.coef_) * (ref_ratio_mean - ams_df["ref_ratio"]) + ams_df["ams"]
        ).astype(int)
        ams_df["ref_ratio"] = round(ams_df["ref_ratio"], 3)
    else:
        ams_df["ams_norm"] = "NA"
        ams_df["ref_ratio"] = "NA"

    ams_df.insert(ams_df.columns.get_loc("ref_ratio") + 1, "ams_norm", ams_df.pop("ams_norm"))
    return ams_df.rename({"ams": "ams_giab"}, axis=1)


def write_normalized_ams_tables(ams_pkls, donor_tables, recipient_tables):
    by_ams_dir = {}
    for ams_pkl in ams_pkls:
        by_ams_dir.setdefault(Path(ams_pkl).parent, []).append(ams_pkl)

    for ams_dir, grouped_pkls in by_ams_dir.items():
        ams_df = normalize_ams(grouped_pkls, donor_tables, recipient_tables)
        ams_df.to_csv(ams_dir / "AMS_df.tsv", sep="\t", index=False)


def write_aams_tables(aams_pkls):
    by_aams_dir = {}
    for aams_pkl in aams_pkls:
        by_aams_dir.setdefault(Path(aams_pkl).parent, []).append(aams_pkl)

    for aams_dir in by_aams_dir:
        create_AAMS_df(aams_dir)


def find_run_path(run_dir, run_name):
    run_dir = Path(run_dir)
    candidates = [
        run_dir / run_name,
        run_dir / "runs" / run_name,
    ]
    if run_dir.name == run_name:
        candidates.append(run_dir)

    for candidate in candidates:
        if candidate.exists():
            return candidate

    matches = [
        path for path in run_dir.rglob(run_name)
        if path.is_dir() and path.parent.name == "runs"
    ]
    if matches:
        return sorted(matches)[0]

    raise FileNotFoundError(f"No such run directory for {run_name} under {run_dir}")


def collect_run_files(run_dirs, run_name):
    ams_pkls = []
    aams_pkls = []
    donor_tables = []
    recipient_tables = []

    for run_dir in run_dirs:
        run_path = find_run_path(run_dir, run_name)

        ams_pkls.extend(str(path) for path in run_path.glob("AMS/**/*.pkl"))
        aams_pkls.extend(str(path) for path in run_path.glob("AAMS/**/*AAMS_df*.pkl"))
        donor_tables.extend(str(path) for path in run_path.glob("run_tables/*_D0_*.tsv"))
        recipient_tables.extend(str(path) for path in run_path.glob("run_tables/*_R0_*.tsv"))

    missing = []
    if not ams_pkls:
        missing.append("AMS pickle files")
    if not donor_tables:
        missing.append("donor tables")
    if not recipient_tables:
        missing.append("recipient tables")
    if missing:
        raise FileNotFoundError("Missing expected output files: " + ", ".join(missing))

    return ams_pkls, aams_pkls, donor_tables, recipient_tables


def merge_run_dirs(run_dirs, run_name, output_dir):
    final_run_dir = Path(output_dir) / "runs" / run_name
    final_run_dir.mkdir(parents=True, exist_ok=True)

    for run_dir in run_dirs:
        source_run_dir = find_run_path(run_dir, run_name)
        shutil.copytree(source_run_dir, final_run_dir, dirs_exist_ok=True)

    return final_run_dir


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--run-dir", nargs="+", required=True)
    args = parser.parse_args()

    run_dir = merge_run_dirs(args.run_dir, args.run_name, args.output_dir)
    ams_pkls, aams_pkls, donor_tables, recipient_tables = collect_run_files([run_dir], args.run_name)
    write_normalized_ams_tables(ams_pkls, donor_tables, recipient_tables)
    write_aams_tables(aams_pkls)


if __name__ == "__main__":
    main()
