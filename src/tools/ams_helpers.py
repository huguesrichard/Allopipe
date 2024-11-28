# coding:utf-8
"""
This file contains all the helpers required to run the ams pipeline
"""
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
from tools import parsing_functions


def create_run_directory(run_name):
    """
    Returns the paths created at the beginning of the run
                    Parameters :
                                    run_name (str): name of the run
                    Returns :
                                    run_path (str): path of the run directory
                                    run_tables (str): path of tables subdir
                                    run_plots (str): path of the plots subdir
                                    run_ams (str): path of the AMS subdir
    """
    # create output directory if it does not exist
    Path("../output").mkdir(parents=True, exist_ok=True)
    # create runs directory, subdict of runs
    Path("../output/runs").mkdir(parents=True, exist_ok=True)
    # create path with run_name
    run_path = os.path.join("../output/runs", run_name)
    # create other needed dependencies for the output
    run_tables = os.path.join(run_path, "run_tables")
    run_plots = os.path.join(run_path, "plots")
    # run_beds = os.path.join(run_path,"bedfiles")
    run_ams = os.path.join(run_path, "AMS")
    # create subdirectories inside of the run_name directory
    Path(run_path).mkdir(parents=True, exist_ok=True)
    Path(run_tables).mkdir(parents=True, exist_ok=True)
    Path(run_plots).mkdir(parents=True, exist_ok=True)
    # Path(run_beds).mkdir(parents=True,exist_ok=True)
    Path(run_ams).mkdir(parents=True, exist_ok=True)
    return run_path, run_tables, run_plots, run_ams


def handle_overwrite(args):
    """
    Returns a boolean after checking if the file already exists
    Parameters:
            run_name (str): string containing the name of the run
    Returns:
            overwrite (bool): True if the file exists, False otherwise
    """
    run_path = f"../output/runs/{args.run_name}"
    if os.path.isdir(run_path) and os.listdir(run_path)!=[]:
        if Path(os.path.join(run_path),f"AMS/{args.run_name}_AMS_{args.min_dp}"
            f"_{args.max_dp}_{args.min_ad}_{args.homozygosity_thr}_{args.min_gq}"
            f"_{args.orientation}_{args.base_length}/"
            f"AMS_{args.run_name}_{args.pair}_{args.min_dp}_{args.max_dp}"
            f"_{args.min_ad}_{args.homozygosity_thr}_{args.min_gq}"
            f"_{args.orientation}_{args.base_length}").is_file():
            overwrite = True
            return overwrite
    return False


#############
# giab v1.3 #
#############
# allow genotypes defined in multi sample vcf to be taken into account (other than 0/0 0/1 1/1)
def explode_gt_ad(df_indiv):
    """
    Returns a subset of a dataframe in exploded format on GT and AD fields
                    Parameters :
                                    df_indiv (pd.DataFrame): dataframe of VCF information with processed sample column
                    Returns :
                                    subset (pd.DataFrame): dataframe containing CHROM, POS, REF, ALT information and GT, AD in exploded format
    """
    # get a subset of the dataframe
    subset = df_indiv[["CHROM", "POS", "REF", "ALT", "GT", "AD"]].copy()
    indiv_name = list(df_indiv.columns)[-1]
    # split Allelic Depth field (format usually as follows for 0/0 : 221,0 for 0/1: 138,116 for 1/1: 0,30)
    subset["AD"] = subset["AD"].apply(lambda x: x.split(","))
    # explode dataframe -> duplicate lines except for the AD field, the previously contained list is exploded
    subset = subset.explode("AD")
    # check and remove lines with "." as AD
    subset = subset[subset["AD"] != "."].copy()
    # change type of AD field from object to int
    subset["AD"] = subset["AD"].astype(int)
    ## operations to keep only the right allelic depths
    subset["ind"] = list(subset.index)
    # index counting the number of values in the AD field, ex : 36,2,5 -> 0 1 2
    subset["ind"] = subset.groupby("ind").cumcount()
    subset = subset.reset_index()
    # proceed exactly the same for GT
    subset["GT"] = subset["GT"].apply(lambda x: x.split("/"))
    subset = subset.explode("GT", ignore_index=True)
    subset["GT"] = subset["GT"].astype(int)
    # keep only the right GT fields, ex : GT=0/2, AD=36,2,15 -> keep 0 = 0 and 2 = 2 (36 and 15 will be kept)
    subset = subset[subset["GT"] == subset["ind"]]
    # remove duplicate rows
    subset = subset.drop_duplicates()
    return subset


# helper function to parse the FORMAT field of the VCF based on indices (useful if FORMAT field is not consistent)
def parse_format(line):
    """
    Returns the indices of the FORMAT field in a line of VCF
                    Parameters :
                                    line (list): list containing the FORMAT fields
                    Returns :
                                    format_indices (list): list containing the indices of interest
    """
    if "FT" in line:
        format_indices = [
            line.index("GT"),
            line.index("GQ"),
            line.index("AD"),
            line.index("DP"),
            line.index("FT"),
        ]
    else:
        format_indices = [
            line.index("GT"),
            line.index("GQ"),
            line.index("AD"),
            line.index("DP"),
        ]
    return format_indices


def apply_format(line, format_indices):
    """
    Returns a concatenated string of all relevant fields of the FORMAT field
                    Parameters :
                                    line (list): list containing the FORMAT field
                                    format_indices (list): list containing the indices of interest
                    Returns :
                                    format_chain (str): concatenated string of FORMAT field
    """
    format_chain = ""
    for i, format_index in enumerate(format_indices):
        if i == len(format_indices) - 1:
            format_chain += line[format_index]
        else:
            format_chain += line[format_index] + "_"
    return format_chain


# pandas.DataFrame -> pandas.DataFrame
def get_read_counts(df_indiv, indiv_file, min_ad, min_gq, base_length):
    """
    Returns the dataframe of the individual and the associated subset, after filtration on the given parameters
                    Parameters :
                                    df_indiv (pd.DataFrame): dataframe of the individual
                                    indiv_file (str): the indiv file name
                                    min_ad (int): minimal Allelic Depth
                                    min_gq (float): minimal Genotype Quality
                                    base_length (int): maximal base length for REF AND ALT
                    Returns:
                            df_indiv (pd.DataFrame): filtered dataframe of the individual
                            subset (pd.DataFrame): filtered subset associated to the dataframe
    """
    # rename the #CHROM column to CHROM
    df_indiv = df_indiv.rename(columns={"#CHROM": "CHROM"})
    # switch type of POS column to int
    df_indiv["POS"] = df_indiv["POS"].astype(int)
    df_indiv = df_indiv.drop_duplicates()
    # remove chr in front of the chromosome number
    df_indiv["CHROM"] = df_indiv["CHROM"].str.replace("chr", "")
    ###############
    # legacy v1.2 #
    ###############
    # filter CHROM column
    df_indiv = df_indiv[
        (df_indiv["CHROM"].str.isdigit())
        | (df_indiv["CHROM"] == "X")
        | (df_indiv["CHROM"] == "Y")
    ]
    # get the name of the indiv
    # indiv_path = "/".join(indiv_file.split("/")[:-1])
    # indiv_pair = indiv_file.split("/")[-2]
    indiv_name = indiv_file.split("/")[-1].split(".")[0]
    if indiv_name not in list(df_indiv.columns):
        indiv_name = list(df_indiv.columns)[-1]
    df_indiv = df_indiv[df_indiv[indiv_name] != "."].reset_index(drop=True).copy()
    df_indiv = df_indiv
    # split and expand the indiv field using FORMAT
    df_indiv["FORMAT"] = df_indiv["FORMAT"].str.split(":")
    df_indiv["select"] = df_indiv[indiv_name].str.split(":")
    df_indiv["split_indices"] = df_indiv["FORMAT"].apply(lambda x: parse_format(x))
    df_indiv["select"] = df_indiv.apply(
        lambda x: apply_format(x["select"], x["split_indices"]), axis=1
    )
    if len(df_indiv["split_indices"][0]) == 5:
        df_indiv[["GT", "GQ", "AD", "DP", "FT"]] = df_indiv["select"].str.split(
            "_", expand=True
        )
        # Filter FT field
        ###############
        # legacy v1.2 #
        ###############
        df_indiv = df_indiv[df_indiv["FT"] == "PASS"]
    else:
        df_indiv[["GT", "GQ", "AD", "DP"]] = df_indiv["select"].str.split(
            "_", expand=True
        )
    df_indiv = df_indiv[~df_indiv["GT"].str.contains(r"\.")]
    # some positions are annotated 0 or 1 and GT field looks like it should be 0/0 and 1/1
    df_indiv["GT"] = df_indiv["GT"].replace("0", "0/0")
    df_indiv["GT"] = df_indiv["GT"].replace("1", "1/1")
    df_indiv = df_indiv.dropna(subset=["GT"])
    df_indiv = df_indiv.drop(["select", "split_indices"], axis=1)
    #############
    # giab v1.3 #
    # #############
    # filter genotype quality
    df_indiv = df_indiv[df_indiv["GQ"].astype(float) > min_gq]
    # filter REF calls
    ###############
    # legacy v1.2 #
    ###############
    # filter REF calls with higher length than 3, value should be selected in cmd line instead of written manually here
    df_indiv = df_indiv[df_indiv["REF"].str.len() <= base_length]
    #### ADD CSV EXCLUDED
    ###############
    # legacy v1.2 #
    ###############
    df_indiv["phased"] = np.where(
        (df_indiv["GT"].str.contains("|")), df_indiv["GT"], np.nan
    )
    # some GT are phased, remove this info in GT field
    df_indiv["GT"] = df_indiv["GT"].str.replace("|", "/", regex=True)
    # when not phased, 1/0 is equivalent to 0/1
    df_indiv["GT"] = df_indiv["GT"].str.replace("1/0", "0/1", regex=True)
    # remove DP field if present
    if "DP" in df_indiv.columns:
        df_indiv = df_indiv.drop("DP", axis=1)
    # ###################################
    # explode df_indiv to get only the right AD values for the given GT
    subset = explode_gt_ad(df_indiv)
    ###############
    # legacy v1.1 #
    ###############
    # filter AD
    subset["AD"] = subset["AD"].astype(int)
    subset = subset[subset["AD"] >= min_ad]
    # filter ALT
    subset["ALT_length"] = np.where(
        subset["GT"] == 0,
        0,
        subset.apply(lambda x: len(x["ALT"].split(",")[x["GT"] - 1]), axis=1),
    )
    #############################################
    sub = subset.groupby(["CHROM", "POS"], as_index=False)["AD"].sum()
    sub = sub.rename({"AD": "sumAD"}, axis=1)
    ###############
    # legacy v1.2 #
    ###############
    ###################################
    subset = subset[subset["ALT_length"] <= base_length]
    ###################################
    # put the DP column for the exploded subset
    subset = pd.merge(subset, sub, how="outer", on=["CHROM", "POS"])
    subset = subset[
        (subset[["CHROM", "POS"]].duplicated(keep=False))
        | (subset["AD"] == subset["sumAD"])
    ]
    subset = subset.drop("sumAD", axis=1)
    sub = subset.groupby(["CHROM", "POS"], as_index=False)["AD"].sum()
    sub = sub.rename({"AD": "DP"}, axis=1)
    subset = pd.merge(subset, sub, how="outer", on=["CHROM", "POS"])
    df_indiv = pd.merge(df_indiv, sub, how="outer", on=["CHROM", "POS"])
    df_indiv = df_indiv.dropna(subset=["DP"])
    return df_indiv, subset


# pandas.DataFrame,int -> pandas.DataFrame
def filter_on_depth(df_indiv, subset, min_dp, max_dp):
    """
    Returns the dataframe of the individual and the associated subset both filtered on DP
                    Parameters :
                                    df_indiv (pd.DataFrame): dataframe of individual
                                    subset (pd.DataFrame): subset associated to the dataframe
                                    min_dp (int): minimal Depth
                                    max_dp (int): maximal Depth
                    Returns :
                                    df_indiv (pd.DataFrame): dataframe filtered on Depth
                                    subset (pd.DataFrame): subset of dataframe filtered on Depth
    """
    ###############
    # legacy v1.1 #
    ###############
    # filter df_indiv and subset on DP
    df_indiv = df_indiv[df_indiv["DP"].between(min_dp, max_dp)]
    subset = subset[subset["DP"].between(min_dp, max_dp)]
    return df_indiv, subset

def convert(df_indiv, subset, homozygosity_threshold):
    """
    Returns the dataframe of the individual, after filtration replacing the GT codes of the heterozygous positions that were above the threshold (symmetric) => converting ambiguous heterozygous positions
                    Parameters :
                                    df_indiv (pd.DataFrame): dataframe of the individual
                                    subset (pd.DataFrame): subset associated to the dataframe
                                    homozygosity_threshold (float): threshold between 0 and 1, ratios above the threshold (symmetric) but called as heterozygous will be converted to homozygous
                    Returns:
                                    df_indiv (pd.DataFrame): dataframe of the individual
    """
    ###############
    # legacy v1.1 #
    ###############
    # create TYPE field
    df_indiv["TYPE"] = np.where(
        df_indiv["GT"].str.split("/").str[0] == df_indiv["GT"].str.split("/").str[1],
        "homozygous",
        "heterozygous",
    )
    subset = pd.merge(
        subset,
        df_indiv[["CHROM", "POS", "DP", "TYPE"]],
        how="outer",
        on=["CHROM", "POS", "DP"],
    )
    # compute allelic ratio
    subset["ratio"] = subset["AD"].astype(int) / subset["DP"]
    # keep only rows where the ratio is higher than (1 - the threshold)
    subset = subset[subset["ratio"] >= 1 - homozygosity_threshold]    
    # change the TYPE column of the rows with ratio >= threshold to homozygous
    subset.loc[
        ~(subset[["CHROM", "POS"]].duplicated(keep=False))
        & (subset["ratio"].between(homozygosity_threshold, 1, inclusive="left")),
        "TYPE",
    ] = "homozygous"
    subset.loc[
        (subset["ratio"] == 1) & (subset["TYPE"] == "heterozygous"), "TYPE"
    ] = "homozygous"
    type_subset = subset.groupby(["CHROM", "POS"], as_index=False)["TYPE"].first()
    # include above changes in individual df
    df_indiv = pd.merge(df_indiv, type_subset, on=["CHROM", "POS"], how="outer")
    df_indiv = df_indiv.drop("TYPE_x", axis=1)
    df_indiv = df_indiv.rename({"TYPE_y": "TYPE"}, axis=1)
    return df_indiv


# pandas.DataFrame,int -> pandas.DataFrame
def get_aa_indiv(df_indiv):
    """
    Returns the dataframe of the individual, with aminoacid information
                    Parameters :
                                    df_indiv (pd.DataFrame): dataframe of the individual
                    Returns:
                                    df_indiv (pd.DataFrame): dataframe of the individual
    """
    # get aa ref and aa alt from VEP information filtered on read counts
    df_indiv["aa_ref_indiv"] = np.where(
        df_indiv["GT"].str.contains("0"), df_indiv["aa_REF"], ""
    )
    df_indiv["aa_alt_indiv"] = np.where(
        (
            ~(df_indiv["GT"].str.contains("0"))
            | (
                df_indiv["GT"].str.split("/").str[0]
                != df_indiv["GT"].str.split("/").str[1]
            )
        ),
        df_indiv["aa_ALT"],
        "",
    )
    df_indiv[["aa_ref_indiv", "aa_alt_indiv"]] = df_indiv[
        ["aa_ref_indiv", "aa_alt_indiv"]
    ].replace("", np.nan)
    df_indiv["aa_indiv"] = df_indiv[["aa_ref_indiv", "aa_alt_indiv"]].apply(
        lambda x: ",".join(x.dropna().values.tolist()), axis=1
    )
    return df_indiv

def filter_on_gnomad_af(df_indiv,min_af):
    """
    Returns the filtered df on gnomad AF information for the individual
    Parameters : 
            df_indiv (pd.DataFrame): dataframe of the individual
            min_af (float): value between 0 and 1 for the minimal accepted AF in the gnomad field of VEP
    Returns :
            df_indiv (pd.DataFrame): filtered dataframe of the individual
    """
    df_indiv = df_indiv[df_indiv["gnomADe_AF"] > min_af]
    return df_indiv

# pandas.DataFrame, str -> pandas.DataFrame
def clean_df(df_indiv, vcf_path_indiv):
    """
    Returns the dataframe of the individual after removing unused columns
                    Parameters :
                                    df_indiv (pd.DataFrame): dataframe of the individual
                                    vcf_path_indiv (str): path of the individual's VCF
                    Returns:
                                    df_indiv (pd.DataFrame): dataframe of the individual
    """
    indiv_name = vcf_path_indiv.split("/")[-1].split(".")[0]
    # remove unwanted known columns
    ls_drop = ["INFO"]
    if indiv_name in df_indiv.columns:
        ls_drop.append(indiv_name)
    df_indiv = df_indiv.drop(ls_drop, axis=1)
    # remove columns containing only zeros
    df_indiv = df_indiv.loc[:, (df_indiv != 0).any(axis=0)]
    # remove rows where there is no aa in both indiv columns
    df_indiv = df_indiv.dropna(subset=["aa_alt_indiv", "aa_ref_indiv"], how="all")
    # remove non rsid SNP
    # df_indiv = df_indiv[df_indiv["ID"].str.contains("rs",na=False)]
    return df_indiv


# build DataFrame from vcf file
# str,int,int -> pandas.DataFrame
def prepare_indiv_df(run_tables, vcf_path_indiv, args, consequences_path, formatted_datetime):
    """
    Returns a dataframe containing all important information, the VEP table and all relevant indices
                    Parameters :
                                    run_tables (str): directory to save the run tables
                                    vcf_path_indiv (str): path of the individual's VCF
                                    args (argparse.Namespace): object containing filtering parameters
                                    consequences_path (str): path of the worst consequences file (WIP)
                    Returns:
                                    df_indiv (pd.DataFrame): dataframe of the individual
                                    vep_table_indiv (pd.DataFrame): dataframe of the VEP information
                                    vep_indices (object): object containing all relevant indices
    """
    # check file extension to select the appropriate parsing function
    if vcf_path_indiv.split(".")[-1] == "vcf":
        df_indiv, vep_indices = parsing_functions.vcf_vep_parser(vcf_path_indiv)
    else:
        df_indiv, vep_indices = parsing_functions.gzvcf_vep_parser(vcf_path_indiv)
    # get read counts info from patient column
    df_indiv, subset = get_read_counts(
        df_indiv, vcf_path_indiv, args.min_ad, args.min_gq, args.base_length
    )
    # filter out rows based on DP min and max
    df_indiv, subset = filter_on_depth(df_indiv, subset, args.min_dp, args.max_dp)
    # convert heterozygote to homozygote if above threshold
    df_indiv = convert(df_indiv, subset, args.homozygosity_thr)
    # parse VEP information and add to dataframe
    df_indiv, vep_table_indiv = parsing_functions.vep_infos_parser(
        run_tables, df_indiv, vep_indices, vcf_path_indiv, args, formatted_datetime
    )
    if args.wc:
        vep_conseq_infos = parsing_functions.worst_consequences_parser(
            consequences_path
        )
        # same chr format for merge
        vep_conseq_infos["CHROM"] = vep_conseq_infos["CHROM"].str.replace("chr", "")
        vep_conseq_infos["POS"] = vep_conseq_infos["POS"].astype(int)
        df_indiv = pd.merge(
            df_indiv, vep_conseq_infos, how="inner", on=["CHROM", "POS"]
        )
        df_indiv["aa_vcf"] = df_indiv["aa_REF"] + "/" + df_indiv["aa_ALT"]
        df_indiv[["aa_REF", "aa_ALT"]] = df_indiv[["aa_REF_vep", "aa_ALT_vep"]]
        print("Running with XXX consequences")
    else:
        # remove above vep_infos_parser ?
        # parse all consequences VEP information and add to dataframe
        # df_indiv = parsing_functions.vep_infos_parser(df_indiv,aa_vep_index)
        print("Running with all consequences")
    # update dataframe with aa ref and aa alt from VEP info
    df_indiv = get_aa_indiv(df_indiv)
    # decide on multiple gnomad
    # print(df_indiv[df_indiv["gnomADe_AF"].str.len()>1])
    df_indiv = df_indiv[df_indiv["gnomADe_AF"].str.len()==1]
    df_indiv["gnomADe_AF"] = df_indiv["gnomADe_AF"].apply(lambda x:x[0])
    df_indiv = df_indiv[df_indiv["gnomADe_AF"]!=""]
    df_indiv["gnomADe_AF"] = df_indiv["gnomADe_AF"].astype(float)
    df_indiv = filter_on_gnomad_af(df_indiv,args.gnomad_af)
    # remove unnecessary columns
    df_indiv = clean_df(df_indiv, vcf_path_indiv)
    return df_indiv, vep_table_indiv, vep_indices


def merge_dfs(df_donor, df_recipient, orientation, imputation):
    """
    Returns a merged dataframe based on the given orientation and some helpers
                    Parameters :
                                    df_donor (pd.DataFrame): dataframe of the donor
                                    df_recipient (pd.DataFrame): dataframe of the recipient
                                    orientation (str): orientation of the merge
                    Returns :
                                    merged_dataframe of the individuals, based on the provided orientation
                                    side (str): reference when looking at duplicated indices within the merged dataframe
                                    opposite (str): opposite to the reference
    """
    if orientation == "dr":
        side = "x"
        opposite = "y"
    else:
        side = "y"
        opposite = "x"
    # merge dataframe from 2 indivs on chosen columns according to the imputation mode
    imputation_type = {
        "imputation": "outer",
        "no-imputation": "inner"
    }[imputation]
    return (
        pd.merge(
            df_donor,
            df_recipient,
            how=imputation_type,
            on=["CHROM", "POS", "REF", "ALT", "aa_REF", "aa_ALT"],
        ),
        side,
        opposite,
    )


def keep_alt(merged_df, side, opposite):
    """
    Returns a merged dataframe where the REF/REF vs NaN positions were removed
                    Parameters :
                                    merged_dataframe of the individuals
                                    side (str): reference when looking at duplicated indices within the merged dataframe
                                    opposite (str): opposite to the reference
                    Returns :
                                    merged_dataframe of the individuals after filtration
    """
    ###############
    # legacy v1.2 #
    ###############
    # remove positions where REF/REF and NaN (we assume it is also REF/REF)
    merged_df = merged_df[
        ~((merged_df[f"GT_{side}"] == "0/0") & (merged_df[f"GT_{opposite}"].isna()))
    ]
    return merged_df


def count_mismatches(merged_df, orientation):
    """
    Returns a merged dataframe and the mismatch count
                    Parameters :
                                    merged_df (pd.DataFrame): merged_dataframe of the individuals
                                    orientation (str): orientation of the merge
                    Returns :
                                    merged_df (pd.DataFrame): dataframe of the individuals
                                    mismatch (int): mismatch count
    """
    ###############
    # legacy v1.2 #
    ###############
    # fill nan with REF for individuals
    merged_df.loc[merged_df["GT_x"].isna(), "aa_indiv_x"] = merged_df["aa_REF"]
    merged_df.loc[merged_df["GT_y"].isna(), "aa_indiv_y"] = merged_df["aa_REF"]
    merged_df[["aa_indiv_x", "aa_indiv_y"]] = merged_df[
        ["aa_indiv_x", "aa_indiv_y"]
    ].fillna("")
    # create a diff column, oriented based on the cmd line argument (kidney : dr, hsc : rd)
    if orientation == "dr":
        merged_df["diff"] = [
            ",".join(set(a.split(",")) - set(b.split(",")))
            for a, b in zip(merged_df["aa_indiv_x"], merged_df["aa_indiv_y"])
        ]
    if orientation == "rd":
        merged_df["diff"] = [
            ",".join(set(b.split(",")) - set(a.split(",")))
            for a, b in zip(merged_df["aa_indiv_x"], merged_df["aa_indiv_y"])
        ]
    merged_df["diff"] = merged_df["diff"].replace("", np.nan)
    merged_df = merged_df.dropna(subset=["diff"]).copy()
    if merged_df.empty:
        mismatch = 0
    else: 
        merged_df["mismatch"] = merged_df["diff"].str.split(",").str.len()
        # create a mismatch column counting the aminoacids in the diff column
        mismatches_df = merged_df.dropna(subset=["diff"])
        mismatches_df["mismatch"] = mismatches_df["diff"].str.split(",").str.len()
        mismatch = mismatches_df["mismatch"].astype(int).sum()
    return merged_df, mismatch


def mismatch_type(merged_df):
    """
    Returns a merged dataframe where mismatches are annotated
                    Parameters :
                                    merged_df (pd.DataFrame): dataframe of the individuals
                    Returns :
                                    merged_df (pd.DataFrame): dataframe of the individuals with mismatch annotation
    """
    ###############
    # legacy v1.2 #
    ###############
    # tag each mismatch in merged_df dataframe with homozygous or heterozygous
    merged_df["mismatch_type"] = np.where(
        ((merged_df["TYPE_x"] == "homozygous") & (merged_df["TYPE_y"] == "homozygous")),
        "homozygous",
        "heterozygous",
    )
    return merged_df


# filter based on fisher test results and binomial
def create_noisy_ams_regions_bed(intervals_table, window, alpha):
    # filter based on number of pairs and p-values
    intervals_table = intervals_table[
        ~(
            (intervals_table["npairs"] <= 15)
            | (intervals_table["pval_GVHa"] < alpha)
            | (intervals_table["pval_GVHc"] < alpha)
            | (intervals_table["pval_GVH"] < alpha)
            | (intervals_table["pval_noGVH"] < alpha)
        )
    ]
    # split the interval column
    intervals_table[["CHROM", "pos_start"]] = intervals_table["interval"].str.split(
        "-", expand=True
    )
    # convert to int type
    intervals_table["pos_start"] = intervals_table["pos_start"].astype(int)
    # get the accurate position numbers
    intervals_table["pos_start"] = intervals_table["pos_start"] * window
    intervals_table["pos_end"] = intervals_table["pos_start"] + window
    # sort ascend
    intervals_table["CHR"] = np.select(
        [
            intervals_table["CHROM"].str.isdigit(),
            intervals_table["CHROM"] == "X",
            intervals_table["CHROM"] == "Y",
        ],
        [intervals_table["CHROM"], "23", "24"],
    )
    intervals_table["CHR"] = intervals_table["CHR"].astype(int)
    intervals_table = intervals_table[
        ["CHROM", "pos_start", "pos_end", "CHR"]
    ].sort_values(by=["CHR", "pos_start"])
    # add chr
    intervals_table["CHROM"] = "chr" + intervals_table["CHROM"]
    # drop CHR column
    intervals_table = intervals_table.drop("CHR", axis=1)
    # save in bed format
    intervals_table.to_csv("noisy_regions.bed", sep="\t", index=False, header=False)
