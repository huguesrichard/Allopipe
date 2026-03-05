# coding:utf-8
"""
This script splits a multi-sample VCF into single VCFs then groupped by pair
"""
import sys
import gzip
import os
from pathlib import Path

# str path -> content of multi-individual vcf : [str,...], str, [str,...]
def get_infos_mvcf(multi_vcf_path):
    if multi_vcf_path.split(".")[-1] == "vcf":
        headers = []
        lines = []
        with open(multi_vcf_path, "r") as file:
            for line in file:
                # get the headers
                if line.startswith("##"):
                    headers.append(line)
                # get column names and content
                else:
                    lines.append(line[:-1])
    else:
        headers = []
        lines = []
        with gzip.open(multi_vcf_path, "rb") as file:
            for line in file:
                # get the headers
                if line.startswith(b"##"):
                    headers.append(line.decode("utf-8"))
                # get column names and content
                else:
                    lines.append(line[:-1].decode("utf-8"))
    return (headers, lines[0], lines[1:])

# str directory_name -> create str path and make directory
def create_dependencies(run_name, output_dir):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    path = os.path.join(output_dir, "runs", run_name, "vcf_indiv")
    Path(path).mkdir(parents=True, exist_ok=True)
    return path

# checks GT field status -> True: isn't empty; False: is empty
def isgenotype(row):
    # get GT field
    if row.split(":")[0] == "./.":
        return False
    else:
        return True
    
# content of multi vcf : [str,...], str, [str,...] -> 1 vcf file per indiv
def write_vcf(headers, names, infos, path, donor_colname, recipient_colname):
    # split names of columns in a list
    ls_names = names.split("\t")
    # get column indices of individuals
    indiv_dic = {donor_colname : ls_names.index(donor_colname), recipient_colname : ls_names.index(recipient_colname)}
    # loop over individuals
    for k,v in indiv_dic.items():
        file_name = path + "/{}.vcf.gz".format(ls_names[v])
        with gzip.open(file_name, "wb") as file:
            # write the headers
            for header in headers:
                file.write(header.encode())
            # write the column names
            file.write(("\t".join(ls_names[0:9]) + "\t" + ls_names[v] + "\n").encode())
            for j in range(len(infos)):
                line_info = infos[j].split("\t")
                # select rows where there is information
                if isgenotype(line_info[v]):
                    # write the row of information
                    file.write((
                        "\t".join(line_info[0:9])
                        + "\t"
                        + line_info[v]
                        + "\n"
                    ).encode())
    

def main(multi_vcf, donor_colname, recipient_colname, run_name, output_dir):
    headers, names, infos = get_infos_mvcf(multi_vcf)
    path = create_dependencies(run_name, output_dir)
    print(f"Extracting {donor_colname} and {recipient_colname}...")
    write_vcf(headers, names, infos, path, donor_colname, recipient_colname)
    path_donor = path + "/{}.vcf.gz".format(donor_colname)
    path_recipient = path + "/{}.vcf.gz".format(recipient_colname)
    return(path_donor, path_recipient)
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
