# coding:utf-8
"""
Script to launch multiple ams pipeline processes
"""
import os
import glob
import concurrent.futures
import time
import re
from tools import parsing_functions, table_operations
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


# get list of mismatches
# get list of transcripts

### HSC genoidentical with gnomad ###
MISMATCHES = [
    file
    for file in glob.glob("../output/runs/gnomad_run/run_tables/*.tsv")
    if "mismatches" in file
]

TRANSCRIPTS = [
    file
    for file in glob.glob("../output/runs/gnomad_run/run_tables/*.tsv")
    if "transcripts" in file
]

# sort lists
MISMATCHES = sorted(MISMATCHES,
    key=lambda x: int(''.join(filter(str.isdigit, str(Path(x).name).split("_")[0]))))
TRANSCRIPTS = sorted(TRANSCRIPTS, 
    key=lambda x: int(''.join(filter(str.isdigit, str(Path(x).name).split("_")[0]))))

print(len(TRANSCRIPTS))
print(len(MISMATCHES))

# parameters
ENSEMBL_TR = "../output/Ensembl/Homo_sapiens.GRCh38.cdna.all.103.fa"
ENSEMBL_PEP = "../output/Ensembl/Homo_sapiens.GRCh38.pep.all.103.fa"
REFSEQ_TR = "../output/Ensembl/Homo_sapiens.GRCh38.103.refseq.tsv"
RUN_NAME = "gnomad_run"
PEP_LENGTH = 9
EL_RANK = 2
HLA_FILE = ("../output/indiv_vcf/joint_genotyping/"
    "hard-filtered/geno-full-typages-HLA.xlsx")
SHEET = "Typage HLA final 16-03-22"
CLASS_I = parsing_functions.read_hla_file(HLA_FILE, SHEET)[0]
GENO_CLASS_I = table_operations.process_geno_hla(CLASS_I)
print(len(GENO_CLASS_I))
# build a list of all commands

### HSC genoidentical normal consequences ###
commands = [
    f"python3 aams_pipeline.py -M {MISMATCH}"
    f" -T {TRANSCRIPT} -E {ENSEMBL_TR} -P {ENSEMBL_PEP}"
    f" -R {REFSEQ_TR} -n {RUN_NAME} -l {PEP_LENGTH}"
    f" --el_rank {EL_RANK} -p {Path(MISMATCH).name.split('_')[0]}"
    f" -a {HLA_PAIR}"
    for (MISMATCH, TRANSCRIPT, HLA_PAIR) in zip(MISMATCHES,TRANSCRIPTS,GENO_CLASS_I)
]

print("commands : ",len(commands))
# print(commands)
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