"""Validate the arguments passed to the scoring scripts before running them."""

import base64
import json
from pathlib import Path
import shlex
import sys

from tools import arguments_handling, netmhc_arguments


def validate(request):
    try:
        count_opts = shlex.split(request["count_opts"])
        affinity_opts = shlex.split(request["affinity_opts"])
    except ValueError as error:
        raise ValueError(f"Invalid quoting in scoring options: {error}") from error

    required_executables = set()
    for pair in request["pairs"]:
        label = pair["pair_id"] or "pair"
        try:
            count = ["-n", request["run_name"], "donor.vcf", "recipient.vcf",
                     request["orientation"], request["imputation"]]
            if request["frameshift"]:
                count.append("--frameshift")
            count += ["--pair", pair["pair_id"], *count_opts, "--mode", request["mode"], "-o", "./"]
            arguments_handling.arguments(count, check_generated_inputs=False)
            affinity = ["-n", request["run_name"], "-p", pair["pair_id"],
                        "-d", request["ensembl_path"], "-a", pair["hla_typing"],
                        *affinity_opts, "-o", "./"]
            affinity_args = netmhc_arguments.netmhc_arguments(affinity, check_generated_inputs=False)
            if not affinity_args.dry_run:
                required_executables.add("netMHCpan" if affinity_args.class_type == 1 else "netMHCIIpan")
            if affinity_args.cleavage:
                required_executables.add("netchop")
        except (SystemExit, ValueError) as error:
            print(f"ERROR: argument validation failed for {label}.", file=sys.stderr, flush=True)
            if isinstance(error, ValueError):
                print(f"ERROR: {error}", file=sys.stderr, flush=True)
            raise
    return {"count_opts": shlex.join(count_opts), "affinity_opts": shlex.join(affinity_opts),
            "required_executables": sorted(required_executables)}


if __name__ == "__main__":
    try:
        options = validate(json.loads(base64.b64decode(sys.argv[1])))
        Path("validated-options.json").write_text(json.dumps(options), encoding="utf-8")
    except (ValueError, KeyError) as error:
        sys.exit(f"ERROR: argument validation failed: {error}")
