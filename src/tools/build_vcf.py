# coding : utf-8
# command line : python3 build_vcf.py <file_name> <directory_name>

import sys
import os
import glob
import re
import argparse
import gzip
from arguments_handling import check_file
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
def create_dependencies(directory_name):
    Path("../../output").mkdir(parents=True, exist_ok=True)
    Path("../../output/indiv_vcf").mkdir(parents=True, exist_ok=True)
    path = "../../output/indiv_vcf/{}".format(directory_name)
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


# content of multi vcf : [str,...], str, [str,...] -> 1 vcf file per indiv
def write_vcf(headers, names, infos, path):
    # split names of columns in a list
    ls_names = names.split("\t")
    # get column indices of individuals
    print(ls_names)
    indiv_columns = [i for i, e in enumerate(ls_names) if "-" in e or "_" in e or "B0" in e]
    print(indiv_columns)
    sys.exit()
    # loop over individuals
    for i in indiv_columns:
        file_name = path + "/{}.vcf.gz".format(ls_names[i])
        with gzip.open(file_name, "wb") as file:
            # write the headers
            for header in headers:
                file.write(header.encode())
            # write the column names
            file.write(("\t".join(ls_names[: indiv_columns[0]]) + "\t" + ls_names[i] + "\n").encode())
            for j in range(len(infos)):
                line_info = infos[j].split("\t")
                # select rows where there is information
                # change this for new vcfs
                if "not-typed" not in line_info[i]:
                    # write the row of information
                    file.write((
                        "\t".join(line_info[: indiv_columns[0]])
                        + "\t"
                        + line_info[i]
                        + "\n"
                    ).encode())


# compress vcf files in given path
def to_gz(path):
    # walk through directories
    for subdir, directory, files in os.walk(path):
        for fname in files:
            fpath = subdir + os.sep + fname
            # find files with .vcf extension
            if fpath.endswith(".vcf"):
                # gzip file without keeping the .vcf file
                os.system("gzip {}".format(fpath))


# move gzipped files in pair directories
def move_files(path):
    # find .vcf and .vcf.gz files
    files = glob.glob(os.path.join(path, "*.vcf.gz")) + glob.glob(
        os.path.join(path, "*.vcf")
    )
    for file in files:
        # find pair name
        pair_dir = re.compile(r"[P|R]\d+").search(file).group()
        # create path
        pair_path = os.path.join(path, pair_dir)
        # create directory of pair if it doesn't exist
        Path(pair_path).mkdir(parents=True, exist_ok=True)
        # move file to directory
        os.system("mv {} {}".format(file, pair_path))


def arguments():
    parser = argparse.ArgumentParser(
        prog="build_vcf.py",
        usage="python3 %(prog)s [options] multi_vcf directory_name",
        description="Split a multi-sample VCF into single VCFs then groupped by pair",
    )
    parser.add_argument(
        "multi_vcf",
        help="multi-sample VCF file, accepted formats are vcf and vcf.gz",
        type=lambda x: check_file(parser, x),
    )
    parser.add_argument(
        "directory_name", help="subdirectory name in the output/indiv_vcf directory"
    )
    return parser.parse_args()


def main(args):
    args = arguments()
    headers, names, infos = get_infos_mvcf(args.multi_vcf)
    path = create_dependencies(args.directory_name)
    write_vcf(headers, names, infos, path)
    # to_gz(path)
    # move_files(path)


if __name__ == "__main__":
    main(sys.argv[1:])
