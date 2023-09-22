# coding:utf-8
"""
Script to launch multiple ams pipeline processes
"""
import os
import glob
import concurrent.futures
import time
import re

# one argument function to be used by ProcessPoolExecutor.map multiprocessing function
def launch_ams_pipeline(command_line):
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
    return f"Done comparing the couple {donor} {recipient}"


# get list of donors
# get list of recipients

### FULL MATCH KIDNEY ###
# donors = [file for file in glob.glob("../output/indiv_vcf/Full_Match/giab_filtered/**/*.vcf.gz") if "Twist_D_" in file and "P21" not in file and "P24" not in file]
# recipients = [file for file in glob.glob("../output/indiv_vcf/Full_Match/giab_filtered/**/*.vcf.gz") if "Twist_R_" in file and "P21" not in file and "P24" not in file]

### HSC genoidentical ###
# donors = [
#     file
#     for file in glob.glob("../output/indiv_vcf/joint_genotyping/hard-filtered/**/*.vcf")
#     if "D0-" in file and "bed." in file
# ]

# recipients = [
#     file
#     for file in glob.glob("../output/indiv_vcf/joint_genotyping/hard-filtered/**/*.vcf")
#     if "R0-" in file and "bed." in file
# ]

### HSC genoidentical with gnomad ###
donors = [
    file
    for file in glob.glob("../output/indiv_vcf/joint_genotyping/hard-filtered/**/*.vcf")
    if "D0-" in file and "bed" in file and "gnomad" in file
]

recipients = [
    file
    for file in glob.glob("../output/indiv_vcf/joint_genotyping/hard-filtered/**/*.vcf")
    if "R0-" in file and "bed" in file and "gnomad" in file
]

### podocytes ###
# donors = [f for f in glob.glob("../output/indiv_vcf/podocytes/**/*.vcf.gz") if not f.split("SRR")[1].split("_")[0][-1].isdigit() and "giab" in f]
# recipients = [f for f in glob.glob("../output/indiv_vcf/podocytes/**/*.vcf.gz") if f.split("SRR")[1].split("_")[0][-1].isdigit() and "giab" in f]

### KIDNEY ###
# donors = [f for f in glob.glob("../output/indiv_vcf/cohort_GQPDOMB/**/*.vcf.gz") if ("D0-" in f or "donor" in f) and ("filtered" in f)]
# recipients = [f for f in glob.glob("../output/indiv_vcf/cohort_GQPDOMB/**/*.vcf.gz") if ("R0-" in f or "recipient" in f) and ("filtered" in f)]

print("donors : ",len(donors))
print("recipients : ",len(recipients))
couples = [
    (f1, f2)
    for f1 in donors
    for f2 in recipients
    if re.compile(r"R\d+|/P\d+").search(f1).group(0)
    == re.compile(r"R\d+|/P\d+").search(f2).group(0)
]
print("couples : ",len(couples))
# parameters
MIN_DP = 20
MAX_DP_CSH = 400
MAX_DP_KID = 500
MIN_AD = 5
HOM_THR = 0.8
GENO_QUAL = 20
ORIENT_KID = "dr"
ORIENT_CSH = "rd"
WORST_DONOR = "../output/GQPDOMB-VEP.txt"
WORST_RECIPIENT = "../output/GQPDOMB-VEP.txt"
RUN_NAME = "gnomad_run"
BASE_LENGTH = 3
# build a list of all commands

### HSC genoidentical normal consequences ###
commands = [
    f"python3 ams_pipeline.py {donor}"
    f" {recipient} {ORIENT_CSH} --min_dp {MIN_DP}"
    f" --max_dp {MAX_DP_CSH} --min_ad {MIN_AD}"
    f" -t {HOM_THR} --min_gq {GENO_QUAL} -l {BASE_LENGTH}"
    f" -f -n {RUN_NAME} -p {'/'.join(donor.split('/')[0:-1]).rsplit('/', maxsplit=1)[-1]}"
    for (donor, recipient) in couples
]

### KIDNEY normal consequences ###
# commands = [f"python3 ams_pipeline.py {donor}"
# f" {recipient} --min_dp {MIN_DP} --max_dp {MAX_DP_KID}"
# f" --min_ad {MIN_AD} -t {HOM_THR} --min_gq {GENO_QUAL} {ORIENT_KID}"
# f" -f -n {RUN_NAME} -p {'/'.join(donor.split('/')[0:-1]).rsplit('/', maxsplit=1)[-1]}"
# for (donor,recipient) in couples]


### KIDNEY --most-severe VEP ###
# commands = [f"python3 comparison_old.py {donor}"
# f" {recipient} {MIN_DP} {MAX_DP_KID}"
# f" {MIN_AD} {HOM_THR} {GENO_QUAL} {ORIENT_KID}"
# f" -f -n {RUN_NAME} -p {'/'.join(donor.split('/')[0:-1]).rsplit('/', maxsplit=1)[-1]}"
# f" -wc -wcd {WORST_DONOR} -wcr {WORST_RECIPIENT}" for (donor,recipient) in couples]

print("commands : ",len(commands))
print(commands)
# start mesuring time
start = time.time()
# multiprocessing
COUNT = 0
with concurrent.futures.ProcessPoolExecutor() as executor:
    # use all commands in list on the function
    results = executor.map(launch_ams_pipeline, commands)
    for res in results:
        COUNT += 1
        print(COUNT)
        print(res)
# stop time count
end = time.time()
print(end - start)