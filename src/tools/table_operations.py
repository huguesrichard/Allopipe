# coding:utf-8
"""
This script contains functions to help handling tables
"""
import pandas as pd
import os
import glob
from pathlib import Path
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
    save_ams = pd.DataFrame([[args.pair if args.pair else "-",
                              args.donor.split('/')[-1].split('.')[0],
                              args.recipient.split('/')[-1].split('.')[0],
                              args.orientation,
                              mismatch,
                              df_donor_file,
                              df_recipient_file,
                              mismatches_file]],
                            columns=["pair", "donor", "recipient", "orientation", "ams",
                                     "donor_table", "recipient_table", "mismatches_table"])
    
    ams_parts = [
        args.run_name if args.run_name else None,  # Include only if not empty
        "AMS",
        args.min_dp,
        args.max_dp,
        args.min_ad,
        args.homozygosity_thr,
        args.min_gq,
        args.orientation,
        args.base_length,
    ]
    # Join only the non-None parts with "_"
    ams_exp_path = os.path.join(
        run_ams,
        "_".join(str(part) for part in ams_parts if part is not None),
    )
    
    Path(ams_exp_path).mkdir(parents=True, exist_ok=True)
    AMS_parts = [
        args.run_name,
        args.pair if args.pair else None,  # include only if not empty
        args.min_dp,
        args.max_dp,
        args.min_ad,
        args.homozygosity_thr,
        args.min_gq,
        args.orientation,
        args.base_length,
    ]
    # join only the non-None parts with "_"
    AMS_fname = os.path.join(
        ams_exp_path,
        "AMS_{}".format("_".join(str(part) for part in AMS_parts if part is not None)),
    )
    save_ams.to_pickle(AMS_fname+formatted_datetime+".pkl")
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
def build_transcripts_table(transcripts_donor, transcripts_recipients):
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
######################### Normalization ############################
####################################################################

def get_ref_ratio_pair(donor_df, recipient_df):
    merged = pd.merge(donor_df, recipient_df, how="outer", on=["CHROM", "POS"])
    common_ref = len(merged[(merged["GT_x"] == "0/0") & (merged["GT_y"] == merged["GT_x"])])
    total_ref = len(merged[(merged["GT_x"] == "0/0") | (merged["GT_y"] == "0/0")])
    if total_ref == 0:
        print("Warning: Normalization cannot be computed - division by 0")
        return (common_ref, total_ref, None)
    else:
        return (common_ref, total_ref, common_ref / total_ref)


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
    print(
        f"Normalizing: "
        f"{donor_path.split('/')[-1].split('.')[0]} "
        f"{recipient_path.split('/')[-1].split('.')[0]}"
    )
    for donor, recipient, ams_pair in zip(donors, recipients, ams):
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
        ams_df = pd.concat(infos).reset_index(drop=True)
    return (ams_df, ams_exp_path, ref_ratio)


def add_norm(ams_df, ams_exp_path, ref_ratio):
    if ref_ratio is None:
        ams_df["ams_norm"] = "NA"
        ams_df["ref_ratio"] = "NA"
    x, y = ams_df["ref_ratio"].values.reshape(-1, 1), ams_df["ams"].values.reshape(-1, 1)
    # compute normalization only if no division by 0
    if ref_ratio is not None:
        reg = lm.LinearRegression().fit(x, y)
        ref_ratio_mean = ams_df["ref_ratio"].mean()
        ams_df["ams_norm"] = (
            float(reg.coef_) * (ref_ratio_mean - ams_df["ref_ratio"]) + ams_df["ams"]
        )
        ams_df["ref_ratio"] = round(ams_df["ref_ratio"], 3)
        ams_df["ams_norm"] = ams_df["ams_norm"].astype(int)
    # reorder "ams_norm" col after "ref_ratio" col
    ams_df.insert(ams_df.columns.get_loc("ref_ratio") + 1, "ams_norm", ams_df.pop("ams_norm"))
    ams_df = ams_df.rename({"ams": "ams_giab"}, axis=1)
    ams_df.to_csv(os.path.join(ams_exp_path, "AMS_df.tsv"), sep="\t", index=False)
