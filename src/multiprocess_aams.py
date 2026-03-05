# coding:utf-8
"""
Script to launch multiple aams pipeline processes
The AAMS pipeline estimates the mismatch between a donor and a recipient after running NetMHCpan
command line help : python multiprocess_aams.py [-h]
"""
import os
import glob
import concurrent.futures
import time
import shlex
from tools import netmhc_arguments
from pathlib import Path


# one argument function to be used by ProcessPoolExecutor.map multiprocessing function
def launch_aams_pipeline(command_line):
    """
    Returns a message saying the couple was successfully compared, after having executed the command line
    Parameters :
                    command_line (str): command line of the ams pipeline
    Returns :
                    (str): message saying the comparison was successfully done
    """
    # execute command line
    os.system(command_line)
    parts = shlex.split(command_line)
    pair = parts[parts.index("-p") + 1] if "-p" in parts else "-"
    return f"Done estimating AAMS for pair {pair}"


def main():
    """
    Script to launch multiple aams pipeline processes
    The AAMS pipeline estimates the mismatch between a donor and a recipient after running NetMHCpan
    command line help : python multiprocess_aams.py [-h]
    """
    print("Pipeline starting...")
    args = netmhc_arguments.netmhc_arguments()
    q = shlex.quote
    run_tables_dir = os.path.join(args.output_dir, "runs", args.run_name, "run_tables")

    # get list of mismatches
    MISMATCHES = [
        file
        for file in glob.glob(os.path.join(run_tables_dir, "*.tsv"))
        if "_mismatches_" in file
    ]

    # get list of transcripts
    TRANSCRIPTS = [
        file
        for file in glob.glob(os.path.join(run_tables_dir, "*.tsv"))
        if "_transcripts_" in file
    ]

    # sort lists
    MISMATCHES.sort()
    TRANSCRIPTS.sort()

    # build a list of all commands
    commands = [
        f"python aams_pipeline.py"
        f" -M {q(MISMATCH)} -T {q(TRANSCRIPT)}"
        f" -n {q(args.run_name)} -d {q(str(args.ensembl_path))}"
        f" --output_dir {q(args.output_dir)}"
        f" --el_rank {args.el_rank} -p {Path(MISMATCH).name.split('_')[0]}"
        f" -a {args.hla_typing}{' --cleavage' if args.cleavage else ''}"
        f"{' --dry_run' if args.dry_run else ''}"
        for (MISMATCH, TRANSCRIPT) in zip(MISMATCHES,TRANSCRIPTS)
    ]

    print("Number of pairs:",len(commands))
    print("Processing transcripts...")

    # start mesuring time
    start = time.time()
    # multiprocessing
    with concurrent.futures.ProcessPoolExecutor(max_workers = args.workers) as executor:
        # use all commands in list on the function
        results = executor.map(launch_aams_pipeline, commands)
    # stop time count
    end = time.time()
    print("Total elapsed time:", int(end - start), "seconds")

if __name__ == "__main__":
    main()
