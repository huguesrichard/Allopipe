# coding : utf-8
import glob
import os
import pandas as pd
import concurrent.futures
import tools.parsing_functions as parsing
from pathlib import Path


def prepare_netmhc_commands(x):
    x = "netMHCpan -BA -f {} -inptype 0 -l 9 -xls -xlsfile {} -a " + x
    return x


def get_commands_list(directory, hla_class_df):
    """
    Input : directory name, class I or class II HLA dataframe
    Output : list of netMHCpan commands
    """
    # add shape of netmhc cmd line to rows in dataframe, adding HLAs
    hla_class_df = hla_class_df.apply(lambda x: prepare_netmhc_commands(x))
    # print(hla_class_df)
    # dataframe to dictionary
    hla_dict = hla_class_df.to_dict(orient="list")
    # get format {str:str} instead of {str:list}
    hla_dict = {
        key: hla_dict[key][0]
        for key in hla_dict
    }
    # get all peptides files designed for netMHCpan
    fa_files = [
        file
        for file in glob.glob(os.path.join(directory, "*netmhc_fasta.fa"))
        for key in hla_dict
        if key.split("-")[0] in file
    ]
    # remove duplicates
    fa_files = list(set(fa_files))
    done = []
    # loop through keys of dictionary
    for key in hla_dict:
        # loop through fasta file paths
        for fasta in fa_files:
            # check if key in fasta file path
            if key.split("-")[0] in fasta:
                # check if fasta path not already in done list
                if fasta not in done:
                    # complete the format fields of the value associated to the key
                    hla_dict[key] = (
                        hla_dict[key].format(
                            fasta,
                            "../output/runs/VIP2/run_tables/netMHCpan_out/"
                            + key
                            + ".out",
                        )
                        + " > /dev/null"
                    )
                    # append to done list
                    done.append(key)

    # keys or not necessary, return only the values (command line for each key)
    return list(hla_dict.values())


# one argument function taking the command line to give to the multiprocessing script
def launch_netmhc(command_line):
    # create the netMHCpan_out directory if it was not created previously
    Path("../output/runs/VIP2/run_tables/netMHCpan_out").mkdir(
        parents=True, exist_ok=True
    )
    os.system(command_line)
    return "End of run"


hla_file = (
    "../output/indiv_vcf/joint_genotyping/hard-filtered/geno-full-typages-HLA.xlsx"
)
sheet = "Typage HLA final 16-03-22"
directory = "../output/runs/VIP2/run_tables"
cI, cII, pivoted = parsing.read_hla_file(hla_file, sheet)
print(cI)
commands = get_commands_list(directory, cI)
print(commands)
print(len(commands))


count = 0
# create multiprocess
with concurrent.futures.ProcessPoolExecutor() as executor:
    # map all commands in list to the launch_netmhc function
    results = executor.map(launch_netmhc, commands)
    for res in results:
        count += 1
        print(count)
        print(res)
