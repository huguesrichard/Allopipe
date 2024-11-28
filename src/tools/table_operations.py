# coding:utf-8
"""
This script contains functions to help handling tables
"""
import pandas as pd
import os
import glob
from pathlib import Path
import re
import sys
import numpy as np
import sklearn.linear_model as lm


################################################################################
############################## AMS table creation ##############################
################################################################################


def save_mismatch(run_ams, args, mismatch, formatted_datetime, df_donor_file, df_recipient_file, mismatches_file):
    """
    Returns the path of the saved file containing the AMS value and related pair
                    Parameters :
                                    run_ams (str): path of the run
                                    args (argparse.Namespace): object containing run parameters
                                    mismatch (int): AMS value
                                    formatted_datetime (str): timestamp if overwrite, empty otherwise
                    Returns :
                                    ams_exp_path (str): path of AMS table
    """
    save_ams = pd.DataFrame([[args.pair,
                              args.donor.split('/')[-1].split('.')[0],
                              args.recipient.split('/')[-1].split('.')[0],
                              args.orientation,
                              mismatch,
                              df_donor_file,
                              df_recipient_file,
                              mismatches_file]],
                            columns=["pair", "donor", "recipient", "orientation", "ams",
                                     "donor_table", "recipient_table", "mismatches_table"])
    
    # change
    ams_exp_path = os.path.join(
        run_ams,
        "{}_AMS_{}_{}_{}_{}_{}_{}_{}".format(
            args.run_name,
            args.min_dp,
            args.max_dp,
            args.min_ad,
            args.homozygosity_thr,
            args.min_gq,
            args.orientation,
            args.base_length,
        ),
    )
    Path(ams_exp_path).mkdir(parents=True, exist_ok=True)
    AMS_fname = os.path.join(
        ams_exp_path,
        "AMS_{}_{}_{}_{}_{}_{}_{}_{}_{}".format(
            args.run_name,
            args.pair,
            args.min_dp,
            args.max_dp,
            args.min_ad,
            args.homozygosity_thr,
            args.min_gq,
            args.orientation,
            args.base_length
        ),
    )
    save_ams.to_pickle(AMS_fname+formatted_datetime)
    save_ams.to_csv(AMS_fname+formatted_datetime+".csv",index=False)
    return ams_exp_path


# create dataframe of mismatch in given directory
def create_AMS_df(ams_exp_path):
    lines = [
        pd.read_pickle(file) for file in glob.glob(os.path.join(ams_exp_path, "*"))
        if ".csv" not in file
    ]
    print(lines)
    df = pd.concat(lines)
    print(df)
    df["value"] = df["pair"].str.split("R|P").str.join("").astype(int)
    df = df.sort_values(by=["value"])
    df = df.drop(["value"], axis=1)
    df.to_csv(os.path.join(ams_exp_path, "AMS_df.tsv"), sep="\t", index=False)


####################################################################
######################## Transcripts tables ########################
####################################################################

# build transcripts table for an individual
def build_transcripts_table_indiv(vep_table_path, merged_ams, vep_indices, indiv_type):
    """
    Returns the vep dataframe of the individual containing the transcripts information
                    Parameters :
                                    vep_table_path (str): path of the table containing vep information
                                    merged_ams (pd.DataFrame): merged dataframe of the individuals
                                    vep_indices (object): object containing all relevant vep indices
                                    indiv_type (str): "donor" or "recipient"
                    Returns :
                                    vep_table (pd.DataFrame): vep table containing transcript information
    """
    # read the VEP table
    vep_table = pd.read_pickle(vep_table_path)
    # merged_ams = pd.read_csv(merged_ams_path,sep="\t")
    if indiv_type == "donor":
        indiv = "x"
    elif indiv_type == "recipient":
        indiv = "y"
    else:
        raise Exception("Error: format donor/recipient incorrect")
        
    # get the positions and the aminoacids associated to the indiv, from the AMS table
    positions = merged_ams[
        [
            "CHROM",
            "POS",
            "aa_ref_indiv_{}".format(indiv),
            "aa_alt_indiv_{}".format(indiv),
            "diff",
        ]
    ].copy()
    # merge the VEP table and the subset from the AMS table
    vep_table = pd.merge(vep_table, positions, how="inner", on=["CHROM", "POS"])
    # explode the dataframe on the INFO column (VEP information only)
    vep_table = vep_table.explode("INFO", ignore_index=True)
    # keep the columns of interest
    vep_table[
        [
            "Consequence",
            "Gene_id",
            "Transcript_id",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "gnomADe_AF"
        ]
    ] = vep_table["INFO"].str.split("|", expand=True)[
        [
            vep_indices.consequence,
            vep_indices.gene,
            vep_indices.transcript,
            vep_indices.cdna,
            vep_indices.cds,
            vep_indices.prot,
            vep_indices.aa,
            vep_indices.codons,
            vep_indices.gnomad
        ]
    ]
    # return the dataframe without the INFO column
    return vep_table.drop("INFO", axis=1)


# build transcripts table for a pair
def build_transcripts_table(transcripts_donor, transcripts_recipients, merged_ams):
    """
    Returns the vep dataframe of the individual containing the transcripts information
                    Parameters :
                                    transcripts_donor (pd.DataFrame): vep table containing transcript information for the donor
                                    transcripts_recipient (pd.DataFrame): vep table containing transcript information for the recipient
                                    merged_ams (pd.DataFrame): merged dataframe of the individuals
                    Returns :
                                    transcripts_pair (pd.DataFrame): all transcripts of the pair
    """
    # merge transcripts donor dataframe with transcripts recipient dataframe
    transcripts_pair = pd.merge(
        transcripts_donor,
        transcripts_recipients,
        how="outer",
        on=[
            "CHROM",
            "POS",
            "Consequence",
            "Gene_id",
            "Transcript_id",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "gnomADe_AF",
            "diff",
        ],
    )
    # drop duplicates
    transcripts_pair = transcripts_pair.drop_duplicates()
    # group by transcript ID
    transcripts_pair.reset_index().groupby(
        [
            "CHROM",
            "POS",
            "Consequence",
            "Gene_id",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "gnomADe_AF",
            "diff",
        ],
        as_index=False,
    )["Transcript_id"].agg(lambda x: x.tolist())
    return transcripts_pair


####################################################################
######################### Table operations #########################
####################################################################

#
def get_relevant_clin_data(clinical_xls):
    clinical = pd.read_csv(clinical_xls, sep="\t")
    clinical = clinical[
        ["Id CryoStem", "GVH aigue", "GVH chronique", "Sexe Donneur", "Sexe Receveur"]
    ]
    clinical = clinical.dropna(subset=["Id CryoStem"])
    clinical = clinical.rename(
        {
            "Id CryoStem": "pair",
            "GVH aigue": "GVHa",
            "GVH chronique": "GVHc",
            "Sexe Donneur": "donor_sex",
            "Sexe Receveur": "recipient_sex",
        },
        axis=1,
    )
    clinical["donor_sex"] = np.where(
        (clinical["donor_sex"] == "Homme"), "Male", "Female"
    )
    clinical["recipient_sex"] = np.where(
        (clinical["recipient_sex"] == "Homme"), "Male", "Female"
    )
    clinical["GVHa"] = np.where((clinical["GVHa"] == 1), "yes", "no")
    clinical["GVHc"] = clinical["GVHc"].fillna("yes")
    clinical["GVHc"] = np.where((clinical["GVHc"] == 1), "yes", "no")
    clinical["pair_sex"] = clinical["donor_sex"] + "-" + clinical["recipient_sex"]
    clinical["sex_mismatch"] = np.where(
        (clinical["donor_sex"] == clinical["recipient_sex"]), "no", "yes"
    )
    clinical["GVH"] = np.where(
        (clinical["GVHa"] == "yes") | (clinical["GVHc"] == "yes"), "yes", "no"
    )
    clinical["GVH_status"] = np.select(
        [
            (clinical["GVHa"] == "yes") & (clinical["GVHc"] == "no"),
            (clinical["GVHa"] == "no") & (clinical["GVHc"] == "yes"),
            (clinical["GVHa"] == "yes") & (clinical["GVHc"] == "yes"),
            (clinical["GVHa"] == "no") & (clinical["GVHc"] == "no"),
        ],
        ["acute_only", "chronic_only", "both_gvh", "no_gvh"],
    )
    clinical = clinical[
        [
            "pair",
            "GVHa",
            "GVHc",
            "GVH",
            "GVH_status",
            "donor_sex",
            "recipient_sex",
            "pair_sex",
            "sex_mismatch",
        ]
    ]
    clinical.to_csv(
        "../../output/general_tables/csh_2021_12_9_relevant_clinical_data.tsv",
        sep="\t",
        index=False,
    )


# Get ref ratio
def get_ref_ratio_pair(donor_df, recipient_df):
    merged = pd.merge(donor_df, recipient_df, how="outer", on=["CHROM", "POS"])
    common_ref = len(merged[(merged["GT_x"] == "0/0") & (merged["GT_y"] == merged["GT_x"])])
    total_ref = len(merged[(merged["GT_x"] == "0/0") | (merged["GT_y"] == "0/0")])
    if total_ref == 0:
        print("Warning: Normalization cannot be computed - division by 0")
        return (common_ref, total_ref, None)
    else:
        return (common_ref, total_ref, common_ref / total_ref)


#
def get_ref_ratio(
    run_name,
    pair_name,
    run_path,
    ams_exp_path,
    donor_path,
    recipient_path,
    min_DP,
    max_DP,
    min_AD,
    homozygosity_threshold,
    min_GQ,
    orientation,
    base_length,
):

    donors = [
        file
        for file in glob.glob(os.path.join(run_path, "**/*.tsv"))
        if ("D0_" in file and "transcripts" not in file and "mismatches" not in file)
    ]
    recipients = [
        file
        for file in glob.glob(os.path.join(run_path, "**/*.tsv"))
        if ("R0_" in file and "transcripts" not in file and "mismatches" not in file)
    ]
    ams_folder = os.path.join(ams_exp_path, "AMS_{}_{}_{}_{}_{}_{}_{}_{}_{}").format(
        run_name, pair_name, min_DP, max_DP, min_AD, homozygosity_threshold, min_GQ, orientation, base_length
    )
    ams = [
    file
    for file in glob.glob(ams_exp_path + "/*")
    if (".csv" not in file and ".tsv" not in file)
    ]
    donors.sort()
    recipients.sort()
    ams.sort()

    infos = []
    for donor, recipient, ams_pair in zip(donors, recipients, ams):
        print(
            ams_pair.split("/")[-1].split("_")[0],
            donor.split("/")[-1].split("_")[2],
            recipient.split("/")[-1].split("_")[2],
        )
        donor_df, recipient_df = pd.read_csv(
            donor, sep="\t", dtype={"CHROM": str}
        ), pd.read_csv(recipient, sep="\t", dtype={"CHROM": str})
        common_ref, total_ref, ref_ratio = get_ref_ratio_pair(donor_df, recipient_df)
        # Add new columns to pickle
        ams_pair_df = pd.read_pickle(ams_pair)
        ams_index = ams_pair_df.columns.get_loc("ams")
        ams_pair_df.insert(ams_index + 1, "common_ref", common_ref)
        ams_pair_df.insert(ams_index + 2, "total_ref", total_ref)
        ams_pair_df.insert(ams_index + 3, "ref_ratio", ref_ratio)
        infos.append(ams_pair_df)
    if ref_ratio is not None:
        ams_df = pd.concat(infos).reset_index(drop=True)
        return (ams_df, ams_exp_path, ref_ratio)
    else:
        return (None, ams_exp_path, ref_ratio)




def add_norm(ams_df, ams_exp_path):
    x, y = ams_df["ref_ratio"].values.reshape(-1, 1), ams_df["ams"].values.reshape(-1, 1)
    reg = lm.LinearRegression().fit(x, y)
    ref_ratio_mean = ams_df["ref_ratio"].mean()
    ams_df["ams_norm"] = (
        float(reg.coef_) * (ref_ratio_mean - ams_df["ref_ratio"]) + ams_df["ams"]
    )
    # reorder "ams_norm" col after "ref_ratio" col
    ams_df.insert(ams_df.columns.get_loc("ref_ratio") + 1, "ams_norm", ams_df.pop("ams_norm"))
    ams_df = ams_df.rename({"ams": "ams_giab"}, axis=1)
    ams_df["ref_ratio"] = round(ams_df["ref_ratio"], 3)
    ams_df["ams_norm"] = ams_df["ams_norm"].astype(int)
    ams_df.to_csv(os.path.join(ams_exp_path, "AMS_df.tsv"), sep="\t", index=False)


#
def merge_clin_ams(
    clinical_df_path,
    run_path,
    ams_folder,
    min_DP,
    max_DP,
    min_AD,
    homozygosity_threshold,
    orientation,
):
    ams = pd.read_csv(
        os.path.join(
            ams_folder,
            "AMS_df.tsv".format(
                min_DP, max_DP, min_AD, homozygosity_threshold, orientation
            ),
        ),
        sep="\t",
    )
    clinical = pd.read_csv(clinical_df_path, sep="\t")
    ams_clin = pd.merge(ams, clinical, how="inner", on="pair")
    ams_clin.to_csv(os.path.join(run_path, "csh_clinical_ams.csv"), index=False)


################################################################################
#################### Shuffle AMS experiment table creation #####################
################################################################################


# save dataframe of shuffle experiment
def shuffle_mismatch(
    run_path,
    run_name,
    file_indiv1,
    file_indiv2,
    mismatch,
    min_DP,
    max_DP,
    min_AD,
    homozygosity_threshold,
):
    pair = file_indiv2.split("/")[-2]
    ind1, ind2 = (
        file_indiv1.split("/")[-1].split(".")[0],
        file_indiv2.split("/")[-1].split(".")[0],
    )
    new_dir = "AMS_shuffle"
    dir_path = os.path.join(run_path, new_dir)
    Path(dir_path).mkdir(parents=True, exist_ok=True)
    save_ams = pd.DataFrame(
        [[ind1, ind2, mismatch]], columns=["donor", "recipient", "AMS"]
    )
    shuffle_ams_path = os.path.join(
        dir_path,
        "{}_AMS_{}_{}_{}_{}_{}".format(
            run_name, min_DP, max_DP, min_AD, homozygosity_threshold, orientation
        ),
    )
    Path(shuffle_ams_path).mkdir(parents=True, exist_ok=True)
    AMS_fname = os.path.join(
        shuffle_ams_path,
        "AMS_{}_{}_{}_{}_{}_{}_{}_{}".format(
            run_name, pair, ind1, ind2, min_DP, max_DP, min_AD, homozygosity_threshold
        ),
    )
    save_ams.to_pickle(AMS_fname)
    return shuffle_ams_path


# create dataframe of all pairs for shuffle experiment
def create_AMS_df_shuffle(shuffle_ams_path):
    files = os.scandir(shuffle_ams_path)
    lines = [pd.read_pickle(file.path) for file in files]
    df = pd.concat(lines)
    df["value"] = df["recipient"].str.extract(r"-*(P\d+|R\d+)-")
    df["value"] = (
        df["recipient"].str.split("R|P").str.join("").str.split("-").str[0].astype(int)
    )
    df = df.sort_values(by="value")
    df = df.drop("value", axis=1)
    df.to_csv(
        os.path.join(shuffle_ams_path, "shuffle_AMS_2021.tsv"), sep="\t", index=False
    )


# !!!!! Must change saving system to the above one !!!!!!


def split_dataframes(shuffle_df):
    recipients = list(shuffle_df["recipient"].unique())
    for r in recipients:
        shuffle_df[shuffle_df["recipient"] == r].to_csv(
            "shuffle_{}_AMS_2021.tsv".format(r), sep="\t"
        )


# sort_shuffle_dataframe(,"../../output/indiv_vcf/joint_genotyping/hard-filtered/AMS_folder/AMS_20_400_5_0.8_rd/AMS_df_sorted.tsv")
def sort_shuffle_dataframe(shuffle_df_path, AMS_df_path):
    shuffle_df = pd.read_csv(shuffle_df_path, sep="\t")
    ams_df = pd.read_csv(AMS_df_path, sep="\t")
    shuffle_df["recipient_pair"] = shuffle_df["recipient"].str.extract(
        r"-*(P\d+|R\d+)-"
    )
    sorted_pairs = ams_df["pair"].tolist()
    sorting_dict = {
        pair: ind for (pair, ind) in zip(sorted_pairs, list(range(len(sorted_pairs))))
    }
    shuffle_df["sorter"] = shuffle_df["recipient_pair"].map(sorting_dict)
    shuffle_df = shuffle_df.sort_values(by="sorter")
    shuffle_df = shuffle_df.drop(["recipient_pair", "sorter"], axis=1)
    f_name, extension = os.path.splitext(os.path.basename(shuffle_df_path))
    sorted_shuffle_df_path = os.path.join(
        os.path.dirname(shuffle_df_path), f_name + "_sorted_on_recipients" + extension
    )
    print(shuffle_df)
    shuffle_df.to_csv(sorted_shuffle_df_path, sep="\t", index=False)
    return sorted_shuffle_df_path


def reshape_dataframe(sorted_recipients_shuffle_df_path, sorted_donors_shuffle_df_path):
    df = pd.read_csv(sorted_recipients_shuffle_df_path, sep="\t")
    df_d = pd.read_csv(sorted_donors_shuffle_df_path, sep="\t")
    sorted_donors = df_d["donor"].unique().tolist()
    sorted_recipients = df["recipient"].unique().tolist()
    df = (
        df.reset_index()
        .groupby(["recipient", "donor"])["AMS"]
        .aggregate("first")
        .unstack()
    )
    df = df.reindex(sorted_recipients)
    df = df.reindex(sorted_donors, axis=1)
    print(df)
    return df


# create a table of all indiv positions after filters
def regroup_dataframes(path, min_DP, max_DP):
    glob_path = os.path.join(path, "**/*.tsv")
    donors = [d for d in glob.glob(glob_path) if "donor" in d or "D0-" in d]
    recipients = [r for r in glob.glob(glob_path) if "recipient" in r or "R0-" in r]
    ls_df = []
    for file in glob.glob(glob_path):
        if (str(min_DP) in file and str(max_DP) in file) and (
            "merged" not in file and "excluded" not in file
        ):
            pair = file.split("/")[-2]
            indiv = "-".join(file.split("/")[-1].split("-")[:3])
            df = pd.read_csv(file, sep="\t")
            df["pair"] = pair
            df["indiv"] = indiv
            if file in donors:
                df["role"] = "donor"
            elif file in recipients:
                df["role"] = "recipient"
            else:
                df["role"] = "undefined"
                print("the role is not defined")
            ls_df.append(df)
    table = pd.concat(ls_df)
    table = table[table.columns.tolist()[-3:] + table.columns.tolist()[:-3]]
    if "Unnamed: 0" in table.columns:
        table = table.drop("Unnamed: 0", axis=1)
    table["value"] = table["pair"].str.split("R|P").str.join("").astype(int)
    table["CHROM"] = table["CHROM"].astype(str)
    table = table.sort_values(by=["value", "CHROM"])
    table = table.drop(["value"], axis=1)
    table.to_csv(
        "../../intermediate_files/table_all_indiv_{}_{}_{}.tsv".format(
            path.split("/")[-1], min_DP, max_DP
        ),
        index=False,
        sep="\t",
    )


# def add_sorted_columns(general_df_path,df_path):
#     gen_df = pd.read_csv(general_df_path,sep="\t")
#     df_to_add = pd.read_csv(,sep="\t")
#     gen_df_merged = pd.merge(gen_df,df_to_add,how="inner",on=["pair"])


def gvh_vs_sex_mismatch(clinical_df_path):
    df = pd.read_csv(clinical_df_path)
    gvha = pd.crosstab(df["pair_sex"], df["GVHa"]).rename_axis(None, axis=1)
    gvhc = pd.crosstab(df["pair_sex"], df["GVHc"]).rename_axis(None, axis=1)
    return (gvha, gvhc)


# !!!!! Must change saving system to the above one !!!!!!


def write_ams_bedfiles(
    directory, min_DP, max_DP, min_AD, homozygosity_threshold, window
):
    files = [
        file
        for file in glob.glob(os.path.join(directory, "**/*.tsv"))
        if min_DP in file
        and max_DP in file
        and min_AD in file
        and homozygosity_threshold in file
        and "bed_merged" in file
        and "TRANSCRIPTS" in file
    ]
    if len(files) != 48:
        print(files)
        sys.exit()
    for file in files:
        # bedfile = open(os.path.join(directory,"ams_bedfile.bed"))
        df = pd.read_csv(file, sep="\t")
        pos_df = df[["CHROM", "POS"]].copy()
        pos_df["CHROM"] = pos_df["CHROM"].astype(str)
        pos_df["CHROM"] = "chr" + pos_df["CHROM"]
        pos_df["POS"] = pos_df["POS"].astype(int)
        pos_df["chromStart"] = pos_df["POS"] - int(window / 2)
        pos_df["chromEnd"] = pos_df["POS"] + int(window / 2)
        pos_df = pos_df.drop("POS", axis=1)
        pos_df.to_csv(
            os.path.join(Path(file).parent, "ams_bedfile_window_10000.bed"),
            sep="\t",
            index=False,
            header=False,
        )
        pos_df.to_csv(
            os.path.join(Path(file).parent, "ams_bedfile_window_10000_headers.bed"),
            sep="\t",
            index=False,
        )
        print(file)


directory = "../../output/indiv_vcf/joint_genotyping/hard-filtered/"
min_DP = "20_"
max_DP = "400"
min_AD = "5"
homozygosity_threshold = "0.8"
window = 10000
# write_ams_bedfiles(directory,min_DP,max_DP,min_AD,homozygosity_threshold,window)


def write_general_ams(directory):
    files = [
        pd.read_csv(file, sep="\t")
        for file in glob.glob(os.path.join(directory, "**/*.bed"))
        if "headers" in file and "10000" in file
    ]
    df = pd.concat(files)
    df = df.drop_duplicates()
    df.to_csv(
        os.path.join(directory, "general_ams_bedfile_10000.bed"),
        sep="\t",
        index=False,
        header=False,
    )


# write_general_ams("../../output/indiv_vcf/joint_genotyping/hard-filtered")


def filter_vcf_regions_ams(directory):
    gz_vcf_files = [
        file
        for file in glob.glob(os.path.join(directory, "**/*.vcf.gz"))
        if "bed" in file
    ]
    # for vcf in gz_vcf_files:
    #     os.system("gzip -dk {}".format(vcf))
    d_vcf_files = [
        file
        for file in glob.glob(os.path.join(directory, "**/*.vcf"))
        if "bed" in file and "D0" in file and "regions" not in file
    ]
    r_vcf_files = [
        file
        for file in glob.glob(os.path.join(directory, "**/*.vcf"))
        if "bed" in file and "R0" in file and "regions" not in file
    ]
    bed_files = [
        file
        for file in glob.glob(os.path.join(directory, "**/*.bed"))
        if "headers" not in file and "_3." in file
    ]
    regex = re.compile(r"/((R|P)\d+)/")
    dirs = [regex.search(file).group(1) for file in bed_files]
    dirs = sorted(dirs, key=lambda x: int("".join(re.split("R|P", x))))
    d_vcf_files = sorted(
        d_vcf_files, key=lambda x: int("".join(re.split("R|P", x.split("/")[-2])))
    )
    r_vcf_files = sorted(
        r_vcf_files, key=lambda x: int("".join(re.split("R|P", x.split("/")[-2])))
    )
    bed_files = sorted(
        bed_files, key=lambda x: int("".join(re.split("R|P", x.split("/")[-2])))
    )
    print(dirs)
    for i in range(len(dirs)):
        d_file = Path(d_vcf_files[i])
        r_file = Path(r_vcf_files[i])
        os.system(
            "bedtools intersect -header -a {} -b {} > {}".format(
                d_vcf_files[i],
                bed_files[i],
                Path.joinpath(
                    d_file.parent, d_file.stem + "_ams_regions" + d_file.suffix
                ),
            )
        )
        os.system(
            "bedtools intersect -header -a {} -b {} > {}".format(
                r_vcf_files[i],
                bed_files[i],
                Path.joinpath(
                    r_file.parent, r_file.stem + "_ams_regions" + r_file.suffix
                ),
            )
        )
        print("vcf filtered for pair {}".format(dirs[i]))
        print(d_file, r_file, bed_files[i])


# filter_vcf_regions_ams("../../output/indiv_vcf/joint_genotyping/hard-filtered")

# concat merge files
def concat_AMS_positions(directory, min_DP, max_DP, min_AD, homozygosity_threshold):
    files = [
        file
        for file in glob.glob(os.path.join(directory, "**/*.tsv"))
        if min_DP in file
        and max_DP in file
        and min_AD in file
        and homozygosity_threshold in file
        and "merged" in file
    ]
    merged = []
    for file in files:
        df = pd.read_csv(file, sep="\t")
        df["pair"] = file.split("/")[-2]
        merged.append(df)
    AMS_pos = pd.concat(merged)
    AMS_pos.to_csv(
        os.path.join(
            directory,
            "csh_19_11_2021_geno_merged_AMS_positions_{}_{}_{}_{}.tsv".format(
                min_DP, max_DP, min_AD, homozygosity_threshold
            ),
        ),
        sep="\t",
        index=False,
    )


# concat_AMS_positions(directory,min_DP,max_DP,min_AD,homozygosity_threshold)


# function to apply to common_intervals dataframe created in R to get IGV coordinates
def write_igv_coordinates(x, window):
    coords = (
        "chr"
        + x.split("-")[0]
        + ":"
        + "{:,}-{:,}".format(
            int(x.split("-")[1]) * window, int(x.split("-")[1]) * window + window
        )
    )
    return coords


#####################################################################
######################## Peptide union check ########################
#####################################################################


# !!!!! Must change saving system to the above one !!!!!!

"""
unused
"""
def get_peptides(directory):
    files = [
        file for file in glob.glob(os.path.join(directory, "**/*.pkl")) if "pep" in file
    ]
    to_concat = []
    for file in files:
        df = pd.read_pickle(file)
        df = df.drop_duplicates(["POS", "peptide"])
        df = df[df["peptide"].str.len() >= 9]
        df["indiv"] = file.split("/")[-2]
        to_concat.append(df)
    concatenated = pd.concat(to_concat)
    concatenated.to_csv(
        os.path.join(directory, "merged_peptide_dataframe.csv"), index=False
    )


#############################################################
######################## AAMS tables ########################
#############################################################


def create_AAMS_df(aams_path):
    lines = [
        pd.read_pickle(file)
        for file in glob.glob(os.path.join(aams_path, "*"))
        if "csv" not in file
    ]
    print(lines)
    df = pd.concat(lines)
    print(df)
    df["value"] = df["pair"].str.split("R|P").str.join("").astype(int)
    df = df.sort_values(by=["value"])
    df = df.drop(["value"], axis=1)
    df.to_csv(os.path.join(aams_path, "AAMS_df.tsv"), sep="\t", index=False)


def process_geno_hla(hla_dataframe):
    hla_dict = hla_dataframe.to_dict(orient="list")
    geno_hla_dict = {}
    for key in hla_dict.keys():
        pair = key.split("-")[0]
        if pair not in geno_hla_dict.keys():
            geno_hla_dict[pair] = hla_dict[key][0]
    geno_hla = list(geno_hla_dict.items())
    geno_hla = sorted(geno_hla,key=lambda x:int(''.join(filter(str.isdigit,x[0]))))
    geno_hla = [tup[1] for i,tup in enumerate(geno_hla)]
    return geno_hla

