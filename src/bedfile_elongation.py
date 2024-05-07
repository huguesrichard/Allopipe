# coding:utf-8
import re
import os
import glob
import pandas as pd
import numpy as np
from pathlib import Path

# warning
# must change the filtering so that the original bedfile contains only exons instead of keeping only the intervals containing the sum of AMS positions
################################################################################
######################## Bedfile with genes of interest ########################
################################################################################

# function to manually take ensembl coords format Chromosome x:xxx,xxx,xxx-xxx,xxx,xxx
def convert_ensembl_coords(coords):
    """
    In : Ensembl coords
    Out : Coords with bed format
    """
    coords = re.split(" |:|-", coords)
    return (
        coords[1]
        + "\t"
        + coords[-2].replace(",", "")
        + "\t"
        + coords[-1].replace(",", "")
        + "\n"
    )


# writes bed format coordinates in file
def add_genes_of_interest(ls_genes, file_name):
    with open(file_name, "w") as f:
        for i in range(len(ls_genes)):
            if i == len(ls_genes) - 1:
                f.write(ls_genes[i][:-1])
            else:
                f.write(ls_genes[i])


# mhags_file was previously obtained on biomart from the file provided
def bedfile_genes(hla_file, mhags_file):
    """
    In : File names of files containing genes of interest
    Out : bedfile of merged regions of all genes of interest
    """
    hla = pd.read_excel(hla_file, sheet_name="HLA")
    hla = hla.dropna(subset=["Ensembl"])
    # kir = pd.read_excel(file_name,sheet_name="KIR")
    mhags = pd.read_csv(mhags_file, sep="\t")
    hla = hla[["Chromosome", "Position"]].copy()
    hla["Chromosome"] = "chr" + hla["Chromosome"].astype(int).astype(str)
    hla[["start", "end"]] = hla["Position"].str.split("-", expand=True)
    hla["Group"] = "HLA"
    #     kir = kir[["Chromosome","Position"]].copy()
    kir = pd.DataFrame(
        [["chr19", "54700000-54900000", "KIR"]],
        columns=["Chromosome", "Position", "Group"],
    )
    kir[["start", "end"]] = kir["Position"].str.split("-", expand=True)
    mhags = mhags[["Chromosome", "start", "end"]].copy()
    mhags["Group"] = "mHAgs"
    genes_coordinates = pd.concat([hla, kir, mhags])
    directory = Path(hla_file).parent
    genes_coordinates[["Chromosome", "start", "end"]].to_csv(
        os.path.join(directory, "genes_of_interest.bed"),
        sep="\t",
        index=False,
        header=False,
    )
    os.system(
        "sort -V -k1,1 -k2,2n {} > {}".format(
            os.path.join(directory, "genes_of_interest.bed"),
            os.path.join(directory, "goi_sorted.bed"),
        )
    )
    os.system(
        "bedtools merge -i {} > {}".format(
            os.path.join(directory, "goi_sorted.bed"),
            os.path.join(directory, "goi_merged.bed"),
        )
    )


#     return()


#########################################################################
######################## bedtools subtract noisy ########################
#########################################################################

# bed subtract cmd line
def subtract_ams_regions(file, bedfile):
    os.system(
        "bedtools subtract -header -a {} -b {} > ../{}_ams_prime.vcf".format(
            file, bedfile, file.split(".")[-2]
        )
    )


# get the giab filtered vcf and filter them using the noisy regions bedfile
def filter_all_vcf(path, bedfile, minDP, maxDP, minAD, homozygosity_threshold, gq):
    files = [
        file
        for file in glob.glob(
            os.path.join(
                path, "**/*.vcf".format(minDP, maxDP, minAD, homozygosity_threshold, gq)
            )
        )
        if "bed" in file
    ]
    for file in files:
        print(file)
        subtract_ams_regions(file, bedfile)


####################################################################
######################## Bedfile elongation ########################
####################################################################


def merged_ams_to_mockbed(merged_ams):
    """
    In : mismatch table of a pair
    Out : Bed-like file containing chr pos pos for each mismatch
    """
    df = pd.read_csv(merged_ams, sep="\t")
    pair = Path(merged_ams).name.split("_")[0]
    mockbed = df[["CHROM", "POS", "POS"]].copy()
    mockbed["CHROM"] = "chr" + mockbed["CHROM"].astype(str)
    path_mockbed = os.path.join(
        Path(merged_ams).parent, "{}_ams_mockbed.bed".format(pair)
    )
    mockbed.to_csv(path_mockbed, sep="\t", index=False, header=False)
    return path_mockbed


def intersect_mockbed(path_mockbed, flank_bed, out_mockbed):
    path_out_mockbed = os.path.join(Path(path_mockbed).parent, out_mockbed)
    os.system(
        "bedtools intersect -a {} -b {} > {}".format(
            path_mockbed, flank_bed, path_out_mockbed
        )
    )
    return path_out_mockbed


def estimate_ams(path_out_mockbed, merged_ams, flank_bed, genome_file):
    mock = pd.read_csv(path_out_mockbed, sep="\t", names=["CHROM", "POS", "dupl_POS"])
    mock["CHROM"] = mock["CHROM"].str.replace("chr", "")
    merged = pd.read_csv(merged_ams, sep="\t")
    merged["CHROM"] = merged["CHROM"].astype(str)
    ams_df = pd.merge(merged, mock, how="inner", on=["CHROM", "POS"])
    ams = ams_df["mismatch"].sum()
    pair = str(Path(path_out_mockbed).parent.absolute()).split("/")[-1]
    df = pd.DataFrame([[pair, ams]], columns=["pair", "AMS"])
    flank_size = re.search("_(\d+)_", flank_bed).group(1)
    quantile = re.search("\d\.\d+", flank_bed).group(0)
    df["flanks"] = flank_size
    df["quantile_filter"] = quantile
    bedcov_file = os.path.join(
        Path(flank_bed).parent, "bedcov_{}_{}.tsv".format(flank_size, quantile)
    )
    os.system(
        "bedtools genomecov -i {} -g {} > {}".format(
            flank_bed, genome_file, bedcov_file
        )
    )
    bedcov = pd.read_csv(
        bedcov_file,
        sep="\t",
        names=["CHROM", "DP", "nb_bases_corresp_DP", "chrom_size", "cov_fraction"],
    )
    bedcov = bedcov[(bedcov["CHROM"] == "genome") & (bedcov["DP"] == 1)]
    df["genome_cov"] = float(bedcov["cov_fraction"])
    return (df, ams)


# merged_ams_to_mockbed("../../output/indiv_vcf/joint_genotyping/hard-filtered/R603/bed_merged_df_20_400_5_0.8_20_TRANSCRIPTS_GENES_STOPLOST.tsv")
# path_mockbed = intersect_mockbed("../../output/indiv_vcf/joint_genotyping/hard-filtered/R603/ams_mockbed.bed","../../output/indiv_vcf/joint_genotyping/hard-filtered/bedfiles/flank_5000/filtered_0.25_giab_without_noisy_regions_5000_merged.bed","ams_merged_mockbed.bed")
# df,ams = estimate_ams(path_mockbed,"../../output/indiv_vcf/joint_genotyping/hard-filtered/R603/bed_merged_df_20_400_5_0.8_20_TRANSCRIPTS_GENES_STOPLOST.tsv","../../output/indiv_vcf/joint_genotyping/hard-filtered/bedfiles/flank_5000/filtered_0.25_giab_without_noisy_regions_5000_merged.bed")
# print(ams)


# flank_bedfiles = [file for file in glob.glob("../../output/indiv_vcf/joint_genotyping/hard-filtered/bedfiles/**/*") if "filtered" in file]
# ls_merged_ams = glob.glob("../../output/indiv_vcf/joint_genotyping/hard-filtered/**/bed_merged_df_20_400_5_0.8_20_TRANSCRIPTS_GENES_STOPLOST.tsv")

# for merged in ls_merged_ams:
#     lst_ams_df = []
#     for flank in flank_bedfiles:
#         path_mockbed = merged_ams_to_mockbed(merged)
#         path_out_mockbed = intersect_mockbed(path_mockbed,flank,"ams_merged_mockbed.bed")
#         df,ams = estimate_ams(path_out_mockbed,merged,flank,"../../output/indiv_vcf/joint_genotyping/hard-filtered/bedfiles/chrom_sizes.txt")
#         print(ams)
#         lst_ams_df.append(df)
#     pair_flanks_df = pd.concat(lst_ams_df)
#     pair_flanks_df.to_csv(os.path.join(Path(path_mockbed).parent,"flanks_with_cov_dataframe.csv"),index=False)

#################################################################
######################## bedtools flanks ########################
#################################################################


# make bigger bedfile intervals
def sloping_regions(bedfile, genome_file, slop_size):
    sloped_bed = os.path.join(
        Path(bedfile).parent, Path(bedfile).stem
    ) + "_{}_slops_ams.bed".format(slop_size)
    merged_bed = os.path.join(
        Path(bedfile).parent, Path(bedfile).stem
    ) + "_{}_merged_ams.bed".format(slop_size)
    os.system(
        "bedtools slop -i {} -g {} -b {} > {}".format(
            bedfile, genome_file, slop_size, sloped_bed
        )
    )
    os.system(
        "sort -V -k1,1 -k2,2n {} > {}_sorted.bed".format(
            sloped_bed, os.path.join(Path(sloped_bed).parent, Path(sloped_bed).stem)
        )
    )
    os.system(
        "bedtools merge -i {}_sorted.bed > {}".format(
            os.path.join(Path(sloped_bed).parent, Path(sloped_bed).stem), merged_bed
        )
    )
    return merged_bed


# def ams_concentration(merged_bedfile,merged_ams):
#     bed = pd.read_csv(merged_bedfile,sep="\t",names=["CHROM","start","end"])
#     ams = pd.read_csv(merged_ams,sep="\t")


def filter_bedfile_intervals(bedfile, threshold):
    bed = pd.read_csv(bedfile, sep="\t", names=["CHROM", "chromStart", "chromEnd"])
    bed[["chromStart", "chromEnd"]] = bed[["chromStart", "chromEnd"]].astype(int)
    bed["size"] = bed["chromEnd"] - bed["chromStart"]
    bed = bed[bed["size"] >= bed["size"].quantile(threshold)]
    bed = bed.drop("size", axis=1)
    bed.to_csv(
        "filtered_{}_{}".format(threshold, bedfile), index=False, header=False, sep="\t"
    )
    return bed


########################################################################
######################## Bedfile elongation v.2 ########################
########################################################################


def make_general_ams_bed(directory, ams_bed):
    # get all the mockbeds containing chrom pos pos for AMS
    files = glob.glob(os.path.join(directory, "**/ams_merged_mockbed.bed"))
    save_dir = os.path.join(directory, "bedfiles")
    ls_df = []
    # read all mockbeds and add to list
    for file in files:
        ls_df.append(pd.read_csv(file, sep="\t", names=["CHROM", "start", "end"]))
    # concat list to get a full dataframe of the AMS positions
    df = pd.concat(ls_df)
    # remove duplicates
    df = df.drop_duplicates()
    # sort by chromosome
    df["CHROM"] = df["CHROM"].str.replace("chr", "")
    df["sorter"] = np.select(
        [df["CHROM"].str.isdigit(), df["CHROM"] == "X", df["CHROM"] == "Y"],
        [df["CHROM"], "23", "24"],
    )
    df["sorter"] = df["sorter"].astype(int)
    df = df.sort_values(by=["CHROM", "start"])
    df = df.drop("sorter", axis=1)
    df["CHROM"] = "chr" + df["CHROM"]
    # save to bed
    df.to_csv(os.path.join(save_dir, ams_bed), sep="\t", index=False, header=False)
    return save_dir


def intersect_with_giab(base_bedfile, ams_all_bed, intersection_bed):
    out_bed_path = os.path.join(Path(base_bedfile).parent, intersection_bed)
    # intersect general mockbed file with base bedfile and returns
    os.system(
        "bedtools intersect -wa -a {} -b {} > {}".format(
            base_bedfile, ams_all_bed, out_bed_path
        )
    )
    df = pd.read_csv(out_bed_path, sep="\t", names=["CHROM", "start", "end"])
    df = df.drop_duplicates()
    df.to_csv(out_bed_path, sep="\t", index=False, header=False)
    return out_bed_path


def cat_genes_of_interest(sloped_merged_bed, goi_bed, full_bed_name):
    out_path = os.path.join(Path(sloped_merged_bed).parent, full_bed_name)
    os.system(
        "cat {} {} | sort -V -k1,1 -k2,2n | bedtools merge -i stdin > {}".format(
            sloped_merged_bed, goi_bed, out_path
        )
    )
    return out_path


def get_ams_ratio(final_bed, ams_pair_mockbed, merged_ams_pair, mock_name):
    mock_path = os.path.join(Path(ams_pair_mockbed).parent, mock_name)
    os.system(
        "bedtools intersect -a {} -b {} > {}".format(
            ams_pair_mockbed, final_bed, mock_path
        )
    )
    mock = pd.read_csv(mock_path, sep="\t", names=["CHROM", "POS", "mockPOS"])
    mock["CHROM"] = mock["CHROM"].str.replace("chr", "")
    ams = pd.read_csv(merged_ams_pair, sep="\t")
    ams["CHROM"] = ams["CHROM"].astype(str)
    new_ams = pd.merge(mock, ams, how="inner", on=["CHROM", "POS"])
    mismatch = new_ams["mismatch"].sum()
    return mismatch


def get_genomecov(final_bed, genome_path, out_name):
    out_path = os.path.join(Path(final_bed).parent, out_name)
    os.system(
        "bedtools genomecov -i {} -g {} > {}".format(final_bed, genome_path, out_path)
    )


############################################################################
######################## Bedfile elongation v.twist ########################
############################################################################


def estimate_cov(bedfile_with_goi, ams_mockbed, merged_ams, mock_name, genome_size):
    mm = get_ams_ratio(bedfile_with_goi, ams_mockbed, merged_ams, mock_name)
    bed = pd.read_csv(bedfile_with_goi, sep="\t", names=["CHROM", "start", "end"])
    ams = pd.read_csv(merged_ams, sep="\t")["mismatch"].sum()
    pair = Path(ams_mockbed).name.split("_")[0]
    flank_size = Path(bedfile_with_goi).name.split("_")[4]
    bed["size"] = bed["end"] - bed["start"]
    coverage = bed["size"].sum() / genome_size
    med = bed["size"].median()
    mean_int = bed["size"].mean()
    return pd.DataFrame(
        [[pair, flank_size, ams, mm, mm / ams, med, mean_int, coverage]],
        columns=[
            "pair",
            "flank_size",
            "AMS",
            "AMS_remaining",
            "AMS_ratio",
            "median_interval",
            "mean_interval",
            "genome_cov",
        ],
    )


# bedfiles = [file for file in glob.glob("../../output/indiv_vcf/joint_genotyping/hard-filtered/bedfiles/Twist/*") if "goi" in file]
# mocks = [file for file in glob.glob("../../output/runs/AMS_twist_giab/run_tables/*") if "mockbed" in file and "pairs" not in file]
# mergeds = [file for file in glob.glob("../../output/runs/AMS_twist_giab/run_tables/*") if "merged" in file and "pairs" not in file]
# sorted_mocks = sorted(mocks,key=lambda x:int(x.split("/")[-1].split("_")[0].split("R")[1]))
# sorted_merged = sorted(mergeds,key=lambda x:int(x.split("/")[-1].split("_")[0].split("R")[1]))
# lst = []
# for bed in bedfiles:
#     for i in range(len(sorted_mocks)):
#         lst.append(estimate_cov(bed,sorted_mocks[i],sorted_merged[i],Path(sorted_mocks[i]).name.split(".")[0]+"_f.bed",3.3*10**9))


def main():
    directory = "../output/indiv_vcf/joint_genotyping/hard-filtered/"

    save_dir = make_general_ams_bed(directory, "general_ams.bed")

    base_bedfile = os.path.join(save_dir, "hg38_giab_haloplex.bed")
    ams_all_bed = os.path.join(save_dir, "general_ams.bed")
    intersection_bed = "selected_ams_giab_halo_hg38.bed"

    out_bed_path = intersect_with_giab(base_bedfile, ams_all_bed, intersection_bed)

    genome_path = os.path.join(save_dir, "chrom_sizes.txt")

    lst_final = []
    for flank in [500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]:
        sloped_merged_bed = sloping_regions(out_bed_path, genome_path, flank)
        out_name = "giab_halo_hg38_{}_bed_final.bed".format(flank)
        out_path = cat_genes_of_interest(
            sloped_merged_bed, os.path.join(save_dir, "goi_merged.bed"), out_name
        )
        lst_final.append(out_path)

    pairs = [p.path for p in os.scandir(directory) if re.search(r"\/R\d+", p.path)]
    lst_remain = []
    for pair in pairs:
        for final_out in lst_final:
            # print(final_out)
            ams_pair_mockbed = os.path.join(pair, "ams_merged_mockbed.bed")
            merged_ams_pair = os.path.join(
                pair, "bed_merged_df_20_400_5_0.8_20_TRANSCRIPTS_GENES_STOPLOST.tsv"
            )
            slop_size = final_out.split("_")[-3]
            mock_name = "remaining_ams_giab_halo_hg38_{}.bed".format(slop_size)
            mismatch = get_ams_ratio(
                final_out, ams_pair_mockbed, merged_ams_pair, mock_name
            )
            ams_old = pd.read_csv(merged_ams_pair, sep="\t")["mismatch"].sum()
            rpair = Path(pair).stem
            lst_remain.append(
                pd.DataFrame(
                    [[rpair, slop_size, mismatch, ams_old, mismatch / ams_old]],
                    columns=["pair", "flanks_size", "AMS_remain", "AMS", "ratio"],
                )
            )
    conc = pd.concat(lst_remain)
    conc.to_csv(
        os.path.join(save_dir, "conc_all_ams_giab_halo_hg38.tsv"), sep="\t", index=False
    )


if __name__ == "__main__":
    main()
