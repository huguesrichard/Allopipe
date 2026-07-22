# coding:utf-8
"""
This script contains functions to help handling tables
"""
import pandas as pd
import os
import glob
import re
from pathlib import Path


################################################################################
############################## AMS table creation ##############################
################################################################################


def save_mismatch(run_ams, args, mismatch, df_donor_file, df_recipient_file, mismatches_file):
    """
    Returns the path of the saved file containing the AMS value and related pair
                    Parameters :
                                    run_ams (str): path of the run
                                    args (argparse.Namespace): object containing run parameters
                                    mismatch (int): AMS value
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
    save_ams.to_pickle(AMS_fname+".pkl")
    save_ams.to_csv(AMS_fname+".csv",index=False)
    return ams_exp_path


def pair_sort_value(pair):
    match = re.search(r"\d+$", str(pair))
    if match:
        return (0, int(match.group()))
    return (1, str(pair))


def create_score_df(score_exp_path, pickle_pattern, output_name):
    lines = [
        pd.read_pickle(file)
        for file in glob.glob(os.path.join(score_exp_path, pickle_pattern))
    ]
    if not lines:
        return None

    df = pd.concat(lines, ignore_index=True)
    df = df.sort_values("pair", key=lambda col: col.map(pair_sort_value))
    output_path = os.path.join(score_exp_path, output_name)
    df.to_csv(output_path, sep="\t", index=False)
    return output_path


# create dataframe of mismatch in given directory
def create_AMS_df(ams_exp_path):
    return create_score_df(ams_exp_path, "*.pkl", "AMS_df.tsv")


def create_AAMS_df(aams_exp_path):
    return create_score_df(aams_exp_path, "*AAMS_df*.pkl", "AAMS_df.tsv")


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
    info_split = vep_table["INFO"].str.split("|", expand=True)
    base_indices = [
        vep_indices.consequence,
        vep_indices.gene,
        vep_indices.transcript,
        vep_indices.cdna,
        vep_indices.cds,
        vep_indices.prot,
        vep_indices.aa,
        vep_indices.codons,
        vep_indices.gnomad,
    ]
    base_cols = [
        "Consequence",
        "Gene_id",
        "Transcript_id",
        "cDNA_position",
        "CDS_position",
        "Protein_position",
        "Amino_acids",
        "Codons",
        "gnomADe_AF",
    ]
    vep_table[base_cols] = info_split[base_indices]
    # Frameshift_sequence is optional
    if vep_indices.frameshift is not None:
        vep_table["Frameshift_sequence"] = info_split[vep_indices.frameshift]
    else:
        vep_table["Frameshift_sequence"] = ""
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
            "Frameshift_sequence",
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
            "Frameshift_sequence",
        ],
        as_index=False,
    )["Transcript_id"].agg(lambda x: x.tolist())
    return transcripts_pair
