# coding:utf-8
"""
Script to launch multiple aams pipeline processes
The AAMS pipeline estimates the mismatch between a donor and a recipient after running NetMHCpan
command line help : python3 multiprocess_aams.py [-h]
"""
import os
import sys
import glob
import concurrent.futures
import time
import re
from tools import parsing_functions, table_operations
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
    # get donor and recipient names
    donor, recipient = command_line.split(" ")[2:4]
    return f"Done estimating AAMS for the couple {donor} {recipient}"


def main():
    """
    Script to launch multiple aams pipeline processes
    The AAMS pipeline estimates the mismatch between a donor and a recipient after running NetMHCpan
    command line help : python3 multiprocess_aams.py [-h]
    """
    args = netmhc_arguments.netmhc_arguments()

    # get list of mismatches
    MISMATCHES = [
        file
        for file in glob.glob(f"../output/runs/{args.run_name}/run_tables/*.tsv")
        if "_mismatches_" in file
    ]

    # get list of transcripts
    TRANSCRIPTS = [
        file
        for file in glob.glob(f"../output/runs/{args.run_name}/run_tables/*.tsv")
        if "_transcripts_" in file
    ]

    # sort lists
    MISMATCHES.sort()
    TRANSCRIPTS.sort()

    print(len(TRANSCRIPTS))
    print(len(MISMATCHES))

    # build a list of all commands
    commands = [
        f"python3 aams_pipeline.py -M {MISMATCH}"
        f" -T {args.transcripts} -E {args.ensembl_transcripts} -P {args.peptides}"
        f" -R {args.refseq} -n {args.run_name} -l {args.length}"
        f" --el_rank {args.el_rank} -p {Path(MISMATCH).name.split('_')[0]}"
        f" -a {args.hla_typing}"
        for (MISMATCH, TRANSCRIPT) in zip(MISMATCHES,TRANSCRIPTS)
    ]

    print("commands : ",len(commands))
    print(commands)

    # start mesuring time
    start = time.time()
    # multiprocessing
    COUNT = 0
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # use all commands in list on the function
        results = executor.map(launch_aams_pipeline, commands)
        for res in results:
            COUNT += 1
            print(COUNT)
            print(res)
    # stop time count
    end = time.time()
    print(end - start)

if __name__ == "__main__":
    main()
