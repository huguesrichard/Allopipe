"""Test-only replacements for VEP and the two expensive Python scoring stages.

Nextflow modules, sample extraction and cohort normalization remain real.
"""
import argparse
import base64
import gzip
import json
import os
from pathlib import Path
import sys

import pandas as pd


def sample_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if not line.startswith("##")]
    header = lines[0].split("\t")
    assert len(header) == 10, "A scoring task must receive a single-sample VCF"
    sample = header[-1]
    rows = []
    for line in lines[1:]:
        cols = line.split("\t")
        gt_index = cols[8].split(":").index("GT")
        rows.append([cols[0], int(cols[1]), cols[9].split(":")[gt_index]])
    return sample, pd.DataFrame(rows, columns=["CHROM", "POS", "GT"])


def run_count(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", required=True)
    parser.add_argument("--pair", nargs="?", default="", const="")
    parser.add_argument("-o", required=True)
    parser.add_argument("--frameshift", action="store_true")
    for name in ["donor", "recipient", "orientation", "imputation"]:
        parser.add_argument(name)
    args = parser.parse_args(argv)
    run = Path(args.o) / "runs" / args.n
    tables, logs, scores = run / "run_tables", run / "logs", run / "AMS" / "test"
    for folder in [tables, logs, scores]:
        folder.mkdir(parents=True, exist_ok=True)
    donor, donor_df = sample_table(args.donor)
    recipient, recipient_df = sample_table(args.recipient)
    prefix = args.pair + "_" if args.pair else ""
    donor_df.to_csv(tables / f"{prefix}D0_{donor}.tsv", sep="\t", index=False)
    recipient_df.to_csv(tables / f"{prefix}R0_{recipient}.tsv", sep="\t", index=False)
    record = {
        "pair": args.pair, "donor": donor, "recipient": recipient,
        "orientation": args.orientation, "imputation": args.imputation,
        "frameshift": args.frameshift,
        "output_dir": os.environ["ALLOPIPE_PUBLISHED_OUTPUT_DIR"],
        "nextflow_command": base64.b64decode(os.environ["ALLOPIPE_NEXTFLOW_COMMAND_BASE64"]).decode(),
    }
    (logs / f"{prefix}count.json").write_text(json.dumps(record), encoding="utf-8")
    score = int(args.pair[1:]) if args.pair else 1
    pd.DataFrame({"pair": [args.pair or "-"], "ams": [score]}).to_pickle(scores / f"{prefix}AMS.pkl")


def run_affinity(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", required=True)
    parser.add_argument("-p", nargs="?", default="", const="")
    parser.add_argument("-d", required=True)
    parser.add_argument("-a", required=True)
    parser.add_argument("-o", required=True)
    args = parser.parse_args(argv)
    run = Path(args.o) / "runs" / args.n
    prefix = args.p + "_" if args.p else ""
    # The affinity task must receive the count results for its own pair.
    record = json.loads((run / "logs" / f"{prefix}count.json").read_text())
    assert record["pair"] == args.p
    record["hla"] = args.a
    record["ensembl_exists"] = Path(args.d).is_dir()
    assert record["ensembl_exists"]
    (run / "logs" / f"{prefix}affinity.json").write_text(json.dumps(record), encoding="utf-8")
    scores = run / "AAMS" / "test"
    scores.mkdir(parents=True, exist_ok=True)
    score = int(args.p[1:]) if args.p else 1
    pd.DataFrame({"pair": [args.p or "-"], "aams": [score * 10]}).to_pickle(
        scores / f"{prefix}AAMS_df.pkl"
    )


def run_vep(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True)
    args, _ = parser.parse_known_args(argv)
    opener = gzip.open if args.i.endswith(".gz") else open
    with opener(args.i, "rt", encoding="utf-8") as source:
        with gzip.open(args.o, "wt", encoding="utf-8") as output:
            output.write(source.read())


if __name__ == "__main__":
    target, *argv = sys.argv[1:]
    if Path(target).name == "ams_pipeline.py":
        run_count(argv)
    elif Path(target).name == "aams_pipeline.py":
        run_affinity(argv)
    elif target == "vep":
        run_vep(argv)
    else:
        # In particular, execute the actual normalize_cohort.py unchanged.
        executable = os.environ["ALLOPIPE_TEST_PYTHON"]
        os.execv(executable, [executable, target, *argv])
