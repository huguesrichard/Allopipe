# coding:utf-8
"""
Script to launch multiple AMS pipeline processes
The AMS pipeline estimates the mismatch between a donor and a recipient based on VCF information
command line help : python3 multiprocess_ams.py [-h]
"""
import os
import concurrent.futures
import time
import sys
from tools import arguments_handling, multivcf_extract, table_operations

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
    
def launch_multivcf_extract(command_line):
    """
    Returns a message saying the couple was successfully extracted, after having executed the command line
    Parameters :
                    command_line (str): command line of the multivcf extraction
    Returns :
                    (str): message saying the vcf were successfully extracted
    """
    # execute command line
    os.system(command_line)
    # get donor and recipient names
    donor, recipient = command_line.split(" ")[-2:]
    return f"Done extracting the couple {donor} {recipient}"    


def read_pairs_info(donors_recipients_table, couples):
    """
    Reads a list of donors and recipients from a table
    Parameters:
                    donors_recipients_table (str): path of table with path of all donors and recipients files
                    couples (list): empty list
    Returns:
                    couples (list): list of couples of donors and recipients
    """
    with open(donors_recipients_table, "r") as f:
        donors, recipients = [], []
        # skips header
        next(f)
        for line in f:
            donor = line.split(",")[0]
            donors.append(donor)
            recipient = line.split(",")[1].rstrip(os.linesep)
            recipients.append(recipient)
            couples.append((donor,recipient))
    return couples


def main():
    """
    Script to launch multiple AMS pipeline processes
    The AMS pipeline estimates the mismatch between a donor and a recipient based on VCF information
    command line help : python3 multiprocess_ams.py [-h]
    """
    args = arguments_handling.arguments(sys.argv[1])
    if args.full:
        full = " --full"
    else:
        full = ""
    couples = []
    read_pairs_info(args.file_pairs, couples)
    
    # extract donor and recipient columns from multi VCF
    path_couples = []
    path_out_couples = multivcf_extract.create_dependencies()
    commands_multivcf = [
    f"python3 tools/multivcf_extract.py {args.multi_vcf} {donor} {recipient}"
    for (donor, recipient) in couples
    ]
    
    for (donor, recipient) in couples:
        path_donor = path_out_couples + "/{}.vcf.gz".format(donor)
        path_recipient = path_out_couples + "/{}.vcf.gz".format(recipient)
        path_couples.append((path_donor, path_recipient))
    print(commands_multivcf)
    
    # start mesuring time
    start = time.time()
    # multiprocessing
    COUNT = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers = args.workers) as executor:
        # use all commands in list on the function
        results = executor.map(launch_multivcf_extract, commands_multivcf)
        for res in results:
            print(res)
            COUNT += 1
            print(COUNT)
            print(res)
    # stop time count
    end = time.time()
    print("Elapsed time:", end - start, "seconds")
    
    # leading zeros for pairs ID : 01 instead of 1 (nb > 10), 001 instead of 1 (nb > 100)
    leading_zeros_number = len(str(len(path_couples)))
    commands = [
    f"python3 ams_pipeline.py {path_donor} {path_recipient} {args.orientation} "
    f"--min_dp {args.min_dp} --max_dp {args.max_dp} --min_ad {args.min_ad} "
    f"--homozygosity_thr {args.homozygosity_thr} --gnomad_af {args.gnomad_af} "
    f"--min_gq {args.min_gq} --base_length {args.base_length} "
    f"--run_name {args.run_name} --pair P{pair_number:0{leading_zeros_number}d}"
    f"{full} --norm_score"
#   Doesn't take into account -wc
#    f"-wc --wc_donor {args.wc_donor} --wc_recipient {args.wc_recipient}"
    for pair_number, (path_donor, path_recipient) in enumerate(path_couples, 1)
    ]
    
    print("couples : ",len(couples))
    print("commands: ",len(commands))
    print(commands)
    
    # start mesuring time
    start = time.time()
    # multiprocessing
    COUNT = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers = args.workers) as executor:
        # use all commands in list on the function
        results = executor.map(launch_ams_pipeline, commands)
        for res in results:
            print(res)
            COUNT += 1
            print(COUNT)
            print(res)
    # stop time count
    end = time.time()
    print("Elapsed time:", end - start, "seconds")

if __name__ == "__main__":
    main()
