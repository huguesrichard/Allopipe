#coding:utf-8
"""
This file contains all the helpers required to run the aams pipeline
"""

import os
import glob
import re
from pathlib import Path
import numpy as np
import pandas as pd
import tools.parsing_functions as parsing


def create_aams_dependencies(ams_run_directory):
    """
    Returns directory paths after completing the creation of the required directories
        Parameters :
                ams_run_directory (str): name of the ams run directory
        Returns:
                aams_run_tables (str): path to the aams_run_tables dir
                netmhc_dir (str): path to the netMHCpan_out dir
                aams_path (str): path to the AAMS dir
    """
    aams_run_tables = f"../output/runs/{ams_run_directory}/aams_run_tables"
    Path(aams_run_tables).mkdir(parents=True, exist_ok=True)
    netmhc_dir = f"../output/runs/{ams_run_directory}/netMHCpan_out"
    Path(netmhc_dir).mkdir(parents=True, exist_ok=True)
    aams_path = f"../output/runs/{ams_run_directory}/AAMS"
    Path(aams_path).mkdir(parents=True, exist_ok=True)
    return aams_run_tables,netmhc_dir,aams_path

def dict_to_df(peptides):
    peptides_ensembl = pd.DataFrame(peptides.items(), columns=["Peptide_id", "INFO"])
    peptides_ensembl[
        ["Gene_id", "Coordinates", "Transcript_id", "Sequence_aa"]
    ] = pd.DataFrame(peptides_ensembl["INFO"].tolist(), index=peptides_ensembl.index)
    peptides_ensembl = peptides_ensembl.drop(["INFO"], axis=1)
    peptides_ensembl[["CHROM", "pos_start", "pos_end", "side"]] = peptides_ensembl[
        "Coordinates"
    ].str.split(":", expand=True)
    peptides_ensembl = peptides_ensembl[
        (peptides_ensembl["CHROM"].str.isdigit())
        | (peptides_ensembl["CHROM"] == "X")
        | (peptides_ensembl["CHROM"] == "Y")
    ]
    return peptides_ensembl

def contributing_ams_transcripts(merged_pair, ensembl_transcripts):
    merged_pair[["transcripts_x", "transcripts_y"]] = merged_pair[
        ["transcripts_x", "transcripts_y"]
    ].fillna("")
    donor_transcripts = {
        a for b in merged_pair["transcripts_x"].str.split(",").tolist() for a in b
    }
    recipient_transcripts = {
        a for b in merged_pair["transcripts_y"].str.split(",").tolist() for a in b
    }
    transcripts = list(set(donor_transcripts) | set(recipient_transcripts))
    print(f"Total potentially contributing transcripts : {len(transcripts)}")
    ams_transcripts = {
        key: value for key, value in ensembl_transcripts.items() if key in transcripts
    }
    return ams_transcripts

def filter_on_refseq(ams_transcripts, refseq_transcripts_path):
    transcripts_df = pd.DataFrame(
        ams_transcripts.items(), columns=["Transcript_id", "Sequence_nt"]
    )
    refseq_transcripts = pd.read_csv(refseq_transcripts_path, sep="\t")
    print(f"REFSEQ transcripts : {len(refseq_transcripts)}")
    refseq_transcripts = refseq_transcripts.rename(
        {"transcript_stable_id": "Transcript_id"}, axis=1
    )
    transcripts_df = pd.merge(
        transcripts_df, refseq_transcripts["Transcript_id"], how="inner"
    ).drop_duplicates()
    return transcripts_df

def intersect_positions(transcripts_pair, transcripts_df):
    positions = transcripts_pair[["CHROM", "POS"]].copy().drop_duplicates()
    print(len(transcripts_pair), len(transcripts_df))
    transcripts_pair = pd.merge(
        transcripts_df, transcripts_pair, how="inner", on="Transcript_id"
    )
    positions_refseq = transcripts_pair[["CHROM", "POS"]].copy().drop_duplicates()
    merged_positions = pd.merge(positions, positions_refseq, how="outer")
    print(len(positions), len(positions_refseq), len(merged_positions))
    return transcripts_pair

def aa_ref(transcripts_pair):
    transcripts_pair["aa_REF"] = np.select(
        [
            transcripts_pair["aa_ref_indiv_x"].notna()
            & transcripts_pair["aa_ref_indiv_y"].isna(),
            transcripts_pair["aa_ref_indiv_y"].notna()
            & transcripts_pair["aa_ref_indiv_x"].isna(),
            transcripts_pair["aa_ref_indiv_x"].notna()
            & transcripts_pair["aa_ref_indiv_y"].notna(),
            transcripts_pair["aa_ref_indiv_x"].isna()
            & transcripts_pair["aa_ref_indiv_y"].isna(),
        ],
        [
            transcripts_pair["aa_ref_indiv_x"],
            transcripts_pair["aa_ref_indiv_y"],
            transcripts_pair["aa_ref_indiv_x"],
            np.nan,
        ],
    )
    return transcripts_pair

def add_pep_seq(transcripts_pair, peptides_ensembl):
    transcripts_pair["CHROM"] = transcripts_pair["CHROM"].astype(str)
    peptides_ensembl["CHROM"] = peptides_ensembl["CHROM"].astype(str)
    transcripts_pair = pd.merge(
        transcripts_pair,
        peptides_ensembl,
        how="inner",
        on=["Gene_id", "Transcript_id", "CHROM"],
    )
    transcripts_pair = transcripts_pair[
        [
            "CHROM",
            "POS",
            "cDNA_position",
            "Protein_position",
            "Consequence",
            "Gene_id",
            "Transcript_id",
            "Sequence_nt",
            "Peptide_id",
            "Sequence_aa",
            "Amino_acids",
            "aa_ref_indiv_x",
            "aa_alt_indiv_x",
            "aa_ref_indiv_y",
            "aa_alt_indiv_y",
            "aa_REF",
            "diff",
        ]
    ]
    return transcripts_pair

def substitution_insertion(x, base, position, pep_size):
    if position - pep_size < 0:
        pep = (
            x[0 : position - 1] + base + x[position + len(base) - 1 : position + pep_size - 1]
        )
    elif position + pep_size - 1 > len(x):
        pep = x[position - pep_size : position - 1] + base + x[position + len(base) - 1 :]
    else:
        pep = (
            x[position - pep_size : position - 1]
            + base
            + x[position + len(base) - 1 : position + pep_size - 1]
        )
    return pep


def deletion(x, base, position, pep_size):
    if position - pep_size < 0:
        pep = x[0:position] + base + x[position : position + pep_size - 1]
    elif position + pep_size - 1 > len(x):
        pep = x[position - pep_size + 1 : position] + base + x[position:]
    else:
        pep = x[position - pep_size + 1 : position] + base + x[position : position + pep_size - 1]
    return pep


def stop(x, position, pep_size):
    if position - pep_size < 0:
        pep = x[0 : position - 1]
    else:
        pep = x[position - pep_size + 1 : position - 1]
    return pep


def mutation_process(x, base, position, pep_size):
    position = int(position)
    if "-" in base:
        pep = deletion(x, base, position, pep_size)
    elif "*" in base:
        pep = stop(x, position, pep_size)
    else:
        pep = substitution_insertion(x, base, position, pep_size)
    return pep


def peptide_seg(peptide, pep_size):
    hla_peptides = [peptide[i : i + pep_size] for i in range(len(peptide) - pep_size + 1)]
    return hla_peptides

def get_peptides_ref(transcripts_pair, pep_size):
    transcripts_pair["aa_REF"] = transcripts_pair["aa_REF"].str.split(",")
    transcripts_pair = transcripts_pair.explode("aa_REF")
    transcripts_pair["diff"] = transcripts_pair["diff"].str.split(",")
    transcripts_pair = transcripts_pair.explode("diff")
    # filter aa longer than 3 aminoacids
    transcripts_pair = transcripts_pair[
        (transcripts_pair["diff"].str.len() <= 3)
        & (transcripts_pair["aa_REF"].str.len() <= 3)
    ]
    transcripts_pair = transcripts_pair.dropna(subset=["Protein_position"])
    # exclude frameshift and stop variants
    transcripts_pair = transcripts_pair[
        ~(
            (transcripts_pair["Consequence"].str.contains("frameshift"))
            | (transcripts_pair["Consequence"].str.contains("stop"))
            | (transcripts_pair["Consequence"].str.contains("splice"))
        )
    ]
    transcripts_pair["Protein_position"] = transcripts_pair["Protein_position"].astype(
        str
    )
    transcripts_pair = transcripts_pair[
        ~(transcripts_pair["Protein_position"].str.contains("-"))
    ]
    transcripts_pair["Protein_position"] = transcripts_pair["Protein_position"].astype(
        float
    )
    transcripts_pair["peptide_REF"] = transcripts_pair.apply(
        lambda x: mutation_process(
            x["Sequence_aa"], x["aa_REF"], x["Protein_position"], pep_size
        ),
        axis=1,
    )
    transcripts_pair["peptide"] = transcripts_pair.apply(
        lambda x: mutation_process(
            x["Sequence_aa"], x["diff"], x["Protein_position"], pep_size
        ),
        axis=1,
    )
    transcripts_pair = transcripts_pair.drop_duplicates()
    transcripts_reduced = transcripts_pair.groupby(
        ["CHROM", "POS", "Gene_id", "peptide_REF", "peptide"], as_index=False
    )[["Transcript_id", "Peptide_id"]].agg("first")
    transcripts_reduced["hla_peptides_REF"] = transcripts_reduced["peptide_REF"].apply(
        lambda x: peptide_seg(x, pep_size)
    )
    transcripts_reduced["hla_peptides"] = transcripts_reduced["peptide"].apply(
        lambda x: peptide_seg(x, pep_size)
    )

    return transcripts_reduced

def create_header_fasta(x):
    if x["peptide"] == x["peptide_REF"]:
        header = (
            f">{x['Gene_id']}:{x['Transcript_id']}:"
            f"{x['Peptide_id']}:{x['CHROM']}:{x['POS']}:REF-"
            f"{x['peptide_REF']}:DIFF-{x['peptide']}:0"
        )
    else:
        header = (
            f">{x['Gene_id']}:{x['Transcript_id']}:"
            f"{x['Peptide_id']}:{x['CHROM']}:{x['POS']}:REF-"
            f"{x['peptide_REF']}:DIFF-{x['peptide']}:1"
        )
    return header

def write_pep_fasta(file_path, transcripts_pair):
    with open(file_path, "w", encoding="utf-8") as file:
        transcripts_pair["header"] = transcripts_pair.apply(
            lambda x: create_header_fasta(x), axis=1
        )
        peptides = dict(
            zip(
                transcripts_pair["header"],
                zip(
                    transcripts_pair["hla_peptides_REF"],
                    transcripts_pair["hla_peptides"],
                ),
            )
        )
        stock = False
        for head, hla_peptides in peptides.items():
            for i, peps in enumerate(hla_peptides):
                if stock:
                    stock = False
                    continue
                file.write(head[:-2] + f":{i+1}")
                file.write("\n")
                for hla_peptide in peps:
                    file.write(hla_peptide)
                    file.write("\n")
                if head[-1] == "0":
                    stock = True

def write_netmhc_fasta(pep_fasta_path, netmhc_fasta_file_name):
    with open(pep_fasta_path, "r", encoding="utf-8") as rf:
        dir_path = str(Path(pep_fasta_path).parent)
        fasta_path = os.path.join(dir_path, netmhc_fasta_file_name)
        with open(
            fasta_path, "w", encoding="utf-8"
        ) as wf:
            for line in rf:
                if ">" in line:
                    start_name = line.split(":")[0]
                    count = 1
                else:
                    wf.write(start_name + f":{count}")
                    wf.write("\n")
                    wf.write(line)
                    count += 1
    return fasta_path

def is_float(element):
    try:
        float(element)
        return True
    except ValueError:
        return False


def get_ams_params(mismatches_path):
    min_dp, max_dp, min_ad, gq, homozygosity_threshold, base_length = list(
        filter(lambda x: is_float(x), str(Path(mismatches_path).stem).split("mismatches")[1].split("_"))
    )
    str_params = (f"mindp_{min_dp}_maxdp_{max_dp}_minad"
    f"_{min_ad}_gq_{gq}_thresh_{homozygosity_threshold}_bl_{base_length}")
    return str_params

def build_peptides(aams_run_tables,str_params,args):
    # get ensembl transcripts from ensembl database
    ens_transcripts = parsing.read_fasta(args.ensembl_transcripts)
    print(f"Ensembl transcripts in dictionary : {len(ens_transcripts.keys())}")
    # get proteins from ensembl database
    proteins = parsing.read_pep_fa(args.peptides)
    print(f"Ensembl proteins in dictionary : {len(proteins.keys())}")
    peptides_ensembl = dict_to_df(proteins)
    print(f"{len(peptides_ensembl)} peptides selected in ensembl refseq")

    # get mismatches file containing ams positions
    mismatches_df = pd.read_csv(args.merged, sep="\t")
    # get list of transcripts in AMS and intersect it with ensembl transcripts
    ams_transcripts = contributing_ams_transcripts(mismatches_df, ens_transcripts)
    print(
        f"Potential contributing transcripts after ensembl filtering : {len(ams_transcripts)}"
    )
    # filtering to keep transcripts present in refseq table
    transcripts_df = filter_on_refseq(ams_transcripts, args.refseq)
    print(f"Transcripts after REFSEQ filter : {len(transcripts_df)}")
    # transcripts long format
    transcripts_pair = pd.read_csv(args.transcripts, sep="\t")
    # intersect filtered transcripts and long format to get a long format with all info
    transcripts_pair = intersect_positions(transcripts_pair, transcripts_df)
    transcripts_pair = aa_ref(transcripts_pair)
    # get pep seq on long format table
    transcripts_pair = add_pep_seq(transcripts_pair, peptides_ensembl)
    transcripts_pair.to_csv(os.path.join(aams_run_tables,f"{args.pair}_{args.run_name}_full.tsv"), sep="\t", index=False)
    transcripts_reduced = get_peptides_ref(transcripts_pair, args.length)
    # orientation = "rd"
    # print(transcripts_reduced)
    pep_indiv_path = os.path.join(aams_run_tables, f"{args.pair}_{args.run_name}_pep_df_{str_params}.pkl")
    transcripts_reduced.to_pickle(
        pep_indiv_path
    )
    write_pep_fasta(
        os.path.join(aams_run_tables, f"{args.pair}_{args.run_name}_kmers.fa"), transcripts_reduced
    )
    fasta_path = write_netmhc_fasta(
        os.path.join(aams_run_tables, f"{args.pair}_{args.run_name}_kmers.fa"),
        f"{args.pair}_{args.run_name}_netmhc_fasta.fa",
    )

    return fasta_path,pep_indiv_path


def run_netmhcpan_class_1(fasta_path,netmhc_dir,args):
    print("Entering netMHCpan handler : running netMHCpan will last a long time")
    netmhc_out = os.path.join(netmhc_dir,args.pair + ".out")
    netmhc_run_output = os.path.join(netmhc_dir,args.pair + "_full_run_information.txt")
    print("class 1")
    os.system(f"netMHCpan -BA -f {fasta_path} -inptype 0 -l {args.length} -xls -xlsfile {netmhc_out} -a {args.hla_typing} > {netmhc_run_output}")
    return netmhc_out
    
def run_netmhcpan_class_2(fasta_path,netmhc_dir,args):
    print("Entering netMHCIIpan handler : running netMHCIIpan will last a long time")
    netmhc_out = os.path.join(netmhc_dir,args.pair + ".out")
    netmhc_run_output = os.path.join(netmhc_dir,args.pair + "_full_run_information.txt")
    print("class 2")
#    os.system(f"netMHCIIpan -BA -f {fasta_path} -inptype 0 -length {args.length} -xls -xlsfile {netmhc_out} -a {args.hla_typing} > {netmhc_run_output}")
    return netmhc_out

# get the peptides table ready to merge
def clean_pep_df(netmhc_table, pep_path,args):
    """
    Returns
    Parameters :
            netmhc_table (pd.DataFrame): NetMHCpan handled table (check netmhc_t ables_handler.py)
            pep_path (str): path to the peptides pickle file 
            args (ArgumentParser object): parser object
    Return: 
    """

    # remove unwanted columns
    columns_to_drop = {
        1: ["Pos", "core", "icore", "Ave"],
        2: ["Pos", "Core", "Ave"]
    }
    netmhc_table = netmhc_table.drop(columns_to_drop.get(args.class_type, []), axis=1)
    # rename columns to match defined nomenclature
    netmhc_table = netmhc_table.rename({"ID": "Gene_id", "Peptide": "hla_peptides"}, axis=1)
    # read peptides indiv file
    pep_df = pd.read_pickle(pep_path)
    # drop unused columns
    pep_df = pep_df.drop(["peptide_REF", "hla_peptides_REF"], axis=1)
    # get the peptides per position by exploding the list of peptides into multiple lines
    pep_df = pep_df.explode("hla_peptides")
    return (netmhc_table, pep_df)


def merge_netmhc(netmhc_df, pep_df, mismatches_path, ELR_thr,pair,aams_path,aams_run_tables,str_params,class_type):
    # merge both pep and netmhc dataframes
    merged = pd.merge(netmhc_df, pep_df, how="inner", on=["Gene_id", "hla_peptides"])
    # group by peptide
    merged = merged.groupby(["peptide"]).agg("first")
    # remove duplicate positions
    merged = merged.drop_duplicates(["CHROM", "POS"])
    ams_df = pd.read_csv(mismatches_path, sep="\t")
    if (
        (
            "X" not in ams_df["CHROM"].unique().tolist()
            and "Y" not in ams_df["CHROM"].unique().tolist()
        )
        and (
            "X" not in pep_df["CHROM"].unique().tolist()
            and "Y" not in pep_df["CHROM"].unique().tolist()
        )
        and (
            "X" not in merged["CHROM"].unique().tolist()
            and "Y" not in merged["CHROM"].unique().tolist()
        )
    ):
        merged[["CHROM", "POS"]] = merged[["CHROM", "POS"]].astype(int)
        pep_df[["CHROM", "POS"]] = pep_df[["CHROM", "POS"]].astype(int)
        ams_df[["CHROM", "POS"]] = ams_df[["CHROM", "POS"]].astype(int)
    merged_aams = pd.merge(ams_df, merged, how="inner", on=["CHROM", "POS"])
    merged_aams = merged_aams.drop_duplicates("hla_peptides")
    column_to_filter = {
        1: "EL_Rank",
        2: "Rank"
    }
    print(merged_aams)
    print(class_type)
    print(ELR_thr)
    merged_aams = merged_aams[merged_aams[column_to_filter[class_type] <= ELR_thr]]
    print(merged_aams)
    mismatch_count = merged_aams["mismatch"].sum()
    aams_df = pd.DataFrame([[pair, mismatch_count]], columns=["pair", "AAMS"])
    aams_dir = os.path.join(aams_path,f"AAMS_{str_params}")
    Path(aams_dir).mkdir(parents=True, exist_ok=True)
    aams_df.to_pickle(
        os.path.join(
            aams_dir,
            pair
            + f"_AAMS_df_{str_params}_el_{ELR_thr}"
        )
    )
    aams_df.to_csv(
        os.path.join(
            aams_dir,
            pair
            + f"_AAMS_df_{str_params}_el_{ELR_thr}.csv"
        )
    )
    merged_aams.to_csv(
        os.path.join(aams_run_tables, f"{pair}_mismatches_aams_EL_{ELR_thr}.tsv"), sep="\t", index=False
    )
    import sys
    sys.exit()
    return mismatch_count
