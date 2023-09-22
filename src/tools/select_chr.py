# coding:utf-8
import sys
import pandas as pd


def select_chr_vcf(vcf_file, chr_nb):
    lines = []
    with open(vcf_file, "r") as f:
        count = 0
        for line in f:
            if "##" in line:
                lines.append(line)
            # get row number and break the loop if column names are found
            if "#CHROM" in line:
                header_index = count
                break
            # increase count otherwise
            else:
                count += 1
    df = pd.read_csv(vcf_file, header=header_index, dtype="str", sep="\t")
    df = df[df["#CHROM"] == chr_nb]
    print(vcf_file)
    write_file = vcf_file.split(".vcf")[0] + "-chr{}".format(chr_nb) + ".vcf"
    print(write_file)
    w = open(write_file, "w")
    for line in lines:
        w.write(line)
    w.close()
    df.to_csv(write_file, mode="a", index=False, sep="\t")


def main(args):
    vcf_file, chr_nb = args[0], args[1]
    select_chr_vcf(vcf_file, chr_nb)


if __name__ == "__main__":
    main(sys.argv[1:])
