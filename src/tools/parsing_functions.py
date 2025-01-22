# coding:utf-8
"""
This script contains all the functions to help parsing the data needed to run the pipeline
"""
import os
import gzip
import sys
from pathlib import Path
import re
import pandas as pd
import numpy as np

#TODO: check better technique for file paths
FILE_CONSEQUENCE = "../data/consequences.txt"

class VepIndices:
    """
    A class used to store VEP indices

    Attributes :
                    consequence (int): index of the consequence field in VEP
                    gene (int): index of the gene field in VEP
                    transcript (int): index of the transcript field in VEP
                    cdna (int): index of the cdna field in VEP
                    cds (int): index of the cds field in VEP
                    prot (int): index of the protein field in VEP
                    aa (int): index of the amino-acid field in VEP
                    codons (int): index of the codons field in VEP
                    gnomad (float): index of the frequency of existing variant in gnomAD exomes combined population in VEP
    """
    def __init__(self, consequence, gene, transcript, cdna, cds, prot, aa, codons, gnomad):
        self.consequence = consequence
        self.gene = gene
        self.transcript = transcript
        self.cdna = cdna
        self.cds = cds
        self.prot = prot
        self.aa = aa
        self.codons = codons
        self.gnomad = gnomad


def vcf_vep_parser(vcf_path):
    """
    Returns a dataframe of the parsed VCF file containing VEP information and the VepIndices object containing the indices
    Parameters :
                    vcf_path (str): path of the VCF file of the individual
    Returns :
                    (pd.DataFrame): dataframe of the individual
                    vep_indices (VepIndices object): object containing the indices for the VEP field parsing
    """
    with open(vcf_path, "r", encoding="utf-8") as file:
        count = 0
        # find line with headers in vcf file
        regex = re.compile('##INFO=<ID=CSQ.*Format:([^"]+)')
        matched = False
        for line in file:
            match = regex.search(line)
            # check if match
            if match:
                # get the indices for desired information
                consequence_index = match.group(1).split("|").index("Consequence")
                gene_index = match.group(1).split("|").index("Gene")
                transcript_index = match.group(1).split("|").index("Feature")
                cdna_index = match.group(1).split("|").index("cDNA_position")
                cds_index = match.group(1).split("|").index("CDS_position")
                prot_index = match.group(1).split("|").index("Protein_position")
                aa_index = match.group(1).split("|").index("Amino_acids")
                codons_index = match.group(1).split("|").index("Codons")
                gnomad_index = match.group(1).split("|").index("gnomADe_AF")
                matched = True
            # get row number and break the loop if column names are found
            if "#CHROM" in line:
                header_index = count
                break
            # increase count otherwise
            count += 1
    if not matched:
        raise ValueError("The provided VCF file does not contain the VEP information !")
    vep_indices = VepIndices(
        consequence_index,
        gene_index,
        transcript_index,
        cdna_index,
        cds_index,
        prot_index,
        aa_index,
        codons_index,
        gnomad_index
    )
    return (
        pd.read_csv(vcf_path, header=header_index, dtype="str", sep="\t"),
        vep_indices,
    )


def gzvcf_vep_parser(vcf_path):
    """
    Returns a dataframe of the parsed gzipped VCF file containing VEP information and the VepIndices object containing the indices
    Parameters :
                    vcf_path (str): path of the gzipped VCF file of the individual
    Returns :
                    (pd.DataFrame): dataframe of the individual
                    vep_indices (VepIndices object): object containing the indices for the VEP field parsing
    """
    file = gzip.open(vcf_path, "rb")
    count = 0
    # find line with headers in vcf file
    regex = re.compile(b'##INFO=<ID=CSQ.*Format:([^"]+)')
    for line in file:
        match = regex.search(line)
        if match:
            consequence_index = (
                match.group(1).decode("utf-8").split("|").index("Consequence")
            )
            gene_index = match.group(1).decode("utf-8").split("|").index("Gene")
            transcript_index = (
                match.group(1).decode("utf-8").split("|").index("Feature")
            )
            cdna_index = (
                match.group(1).decode("utf-8").split("|").index("cDNA_position")
            )
            cds_index = match.group(1).decode("utf-8").split("|").index("CDS_position")
            prot_index = (
                match.group(1).decode("utf-8").split("|").index("Protein_position")
            )
            aa_index = match.group(1).decode("utf-8").split("|").index("Amino_acids")
            codons_index = match.group(1).decode("utf-8").split("|").index("Codons")
            gnomad_index = match.group(1).decode("utf-8").split("|").index("gnomADe_AF")
        # get row number and break the loop if column names are found
        if b"#CHROM" in line:
            header_index = count
            break
        # increase count otherwise
        count += 1
    file.close()
    vep_indices = VepIndices(
        consequence_index,
        gene_index,
        transcript_index,
        cdna_index,
        cds_index,
        prot_index,
        aa_index,
        codons_index,
        gnomad_index,
    )
    return (
        pd.read_csv(vcf_path, header=header_index, dtype="str", sep="\t"),
        vep_indices,
    )


def extract_aa_from_vep(df_infos, vep_indices):
    """
    Returns a dataframe containing the amino-acid information of VEP
    Parameters :
                    df_infos (pd.DataFrame): dataframe containing VEP information
                    vep_indices (VepIndices object): object containing the indices for the VEP field parsing
    Returns :
                    selected_df (pd.DataFrame): dataframe containing the amino-acid information
    """
    # get the selected consequences from file
    with open(FILE_CONSEQUENCE, encoding="utf-8") as file:
        consequences = file.readline().split("\t")
    for conseq in consequences:
        # count number of times consequence appears in row and add count to the corresponding column
        df_infos[conseq] = df_infos["INFO"].str.count(conseq)
    # might want to create a sublist of consequences that we will not consider
    # select a subset based on the unwanted consequences
    selected_df = df_infos[
        df_infos[consequences].idxmax(axis=1) != "synonymous_variant"
    ].copy()
    vep_split_df = selected_df["INFO"].str.split(pat=",", expand=True)
    len_expand = len(list(vep_split_df.columns))
    # get gene ENSG ID using found index in the vep line parser
    gene_subset = vep_split_df[list(range(len_expand))].apply(
        lambda x: x.str.split("|").str[vep_indices.gene]
    )
    gene_subset["gene"] = (
        gene_subset.stack()
        .groupby(level=0)
        .apply(lambda x: x.dropna().unique().tolist())
    )
    # get transcript ENST ID using found index in the vep line parser
    transcript_subset = vep_split_df[list(range(len_expand))].apply(
        lambda x: x.str.split("|").str[vep_indices.transcript]
    )
    transcript_subset["transcript"] = (
        transcript_subset.stack()
        .groupby(level=0)
        .apply(lambda x: x.dropna().unique().tolist())
    )
    # get aa ref and aa alt from VEP info
    aa_ref_subset = vep_split_df[list(range(len_expand))].apply(
        lambda x: x.str.split("|").str[vep_indices.aa].str.split("/").str[0]
    )
    aa_alt_subset = vep_split_df[list(range(len_expand))].apply(
        lambda x: x.str.split("|").str[vep_indices.aa].str.split("/").str[1]
    )
    aa_ref_subset = aa_ref_subset.replace("", np.nan)
    aa_ref_subset["aa"] = (
        aa_ref_subset.stack()
        .groupby(level=0)
        .apply(lambda x: x.dropna().unique().tolist())
    )
    # print(aa_ref_subset)
    aa_alt_subset["aa"] = (
        aa_alt_subset.stack()
        .groupby(level=0)
        .apply(lambda x: x.dropna().unique().tolist())
    )
    gnomad_subset = vep_split_df[list(range(len_expand))].apply(
        lambda x: x.str.split("|").str[vep_indices.gnomad])
    gnomad_subset["gnomADe_AF"] = (
        gnomad_subset.stack()
        .groupby(level=0)
        .apply(lambda x: x.dropna().unique().tolist())
    )
    selected_df["transcripts"] = transcript_subset["transcript"].str.join(",")
    selected_df["genes"] = gene_subset["gene"].str.join(",")
    selected_df["aa_REF"] = aa_ref_subset["aa"].str.join(",")
    selected_df["aa_ALT"] = aa_alt_subset["aa"].str.join(",")
    selected_df["gnomADe_AF"] = gnomad_subset["gnomADe_AF"]
    return selected_df

def vep_infos_parser(run_tables, df_indiv, vep_indices, vcf_path_indiv, args, formatted_datetime):
    """
    Returns the dataframe of the individual after adding VEP information and the path of the individual vep table
    Parameters :
                    run_tables (str): path of the directory containing the run tables
                    df_indiv (pd.DataFrame): dataframe of the individual
                    vep_indices (VepIndices object): object containing the indices for the VEP field parsing
                    vcf_path_indiv (str): path of the VCF file of the individual
                    args (argparse.Namespace): object containing arguments from the command line
                    formatted_datetime (str): timestamp if overwrite, empty otherwise
    Returns :
                    (pd.DataFrame): dataframe containing the VEP information
                    vep_table_indiv (str): path of the individual vep table
    """
    # redefine INFO column based on CSQ information
    df_indiv["INFO"] = df_indiv["INFO"].str.extract("CSQ=([^;]+)", expand=False)
    # subset of dataframe
    df_infos = df_indiv[["CHROM", "POS", "INFO"]].copy()
    save = df_infos.copy()
    save["INFO"] = save["INFO"].str.split(",")
    indiv = Path(vcf_path_indiv).name.split(".")[0]
    vep_table_indiv = os.path.join(
        run_tables,
        f"{indiv}_vep_infos_table_{args.min_dp}_"
        f"{args.max_dp}_{args.min_ad}_{args.homozygosity_thr}{formatted_datetime}.pkl",
    )
    save.to_pickle(vep_table_indiv)
    # return an inner merge from the gained infos from group_from_vep function and the df_indiv dataframe
    return (
        pd.merge(
            df_indiv,
            extract_aa_from_vep(df_infos, vep_indices),
            how="inner",
            on=["CHROM", "POS", "INFO"],
        ),
        vep_table_indiv,
    )


def worst_consequences_parser(consequences_path):
    """
    Returns the VEP dataframe with the worst consequences per position only
    Parameters :
                    consequences_path (str): path to the worst consequences VEP file
    Returns :
                    vep_infos (pd.DataFrame): dataframe containing the VEP information with the worst consequences
    """
    with open(consequences_path, "r", encoding="utf-8") as file:
        count = 0
        for line in file:
            if "#Uploaded_variation" in line:
                header_index = count
                break
            count += 1
    # read the file using the found header_index to get the column names
    df_conseq = pd.read_csv(
        consequences_path, header=header_index, dtype="str", sep="\t"
    )
    df_conseq = df_conseq.replace("-", np.nan)
    df_conseq = df_conseq.dropna(subset=["Amino_acids"])
    df_conseq[["CHROM", "POS", "nt"]] = df_conseq["#Uploaded_variation"].str.split(
        "_", expand=True
    )
    # check if filter on LOW IMPACT
    df_conseq = df_conseq[df_conseq["Consequence"] != "synonymous_variant"].copy()
    vep_infos = (
        df_conseq.groupby(["CHROM", "POS", "nt"])["Amino_acids"]
        .apply(lambda x: ",".join(x.unique()))
        .reset_index()
    )
    vep_subset = vep_infos["Amino_acids"].str.split(",", expand=True).fillna(np.nan)
    vep_subset_ref = vep_subset[list(range(len(list(vep_subset.columns))))].apply(
        lambda x: x.str.split("/").str[0]
    )
    vep_subset_ref = (
        vep_subset_ref.stack()
        .groupby(level=0)
        .apply(lambda x: x.dropna().unique().tolist())
    )
    vep_subset_alt = vep_subset[list(range(len(list(vep_subset.columns))))].apply(
        lambda x: x.str.split("/").str[1]
    )
    vep_subset_alt = (
        vep_subset_alt.stack()
        .groupby(level=0)
        .apply(lambda x: x.dropna().unique().tolist())
    )
    vep_infos["aa_REF_vep"] = vep_subset_ref.str.join(",")
    vep_infos["aa_ALT_vep"] = vep_subset_alt.str.join(",")
    # print(df_conseq["Location","CHROM","POS","Amino_acids"])
    return vep_infos


#############################################
#################### AAMS ###################
#############################################


def read_fasta(file):
    """
    Returns the dictionary of transcripts from the fasta file
    Parameters :
                    file (str): path of the fasta file
    Returns :
                    transcripts (dict): dictionary containing all fasta entries and sequences
    """
    # open file
    with open(file, "r", encoding="utf-8") as fasta_file:
        # create transcripts dictionary
        transcripts = {}
        # initiate bool to tag the first iteration
        first = True
        # loop through the file
        for line in fasta_file:
            # check if line contains a > (transcript_id) or \n (end of file)
            if ">" in line or line == "\n":
                # first iteration
                if first:
                    # get the transcript_id
                    transcript_id = line[1:].split(" ")[0].split(".")[0]
                    first = False
                # other iterations
                else:
                    # check if transcript from the previous iteration has an entry in the dictionary
                    if transcript_id not in transcripts:
                        # add the sequence to the dictionary
                        transcripts[transcript_id] = transcript_seq
                    else:
                        print("transcript already has a sequence")
                        sys.exit()
                    # save the transcript_id of the new transcript
                    transcript_id = line[1:].split(" ")[0].split(".")[0]
                # initiate the transcript sequence
                transcript_seq = ""
            else:
                # add the transcript sequence to the line from previous iteration
                transcript_seq += line[:-1]
        # ensure to add last transcript if last line of file is not "\n"
        if transcript_id not in transcripts:
            transcripts[transcript_id] = transcript_seq
    return transcripts


def read_pep_fa(file):
    """
    Returns the dictionary of transcripts from the peptide fasta file
    Parameters :
                    file (str): path of the peptide fasta file
    Returns :
                    transcripts (dict): dictionary containing all fasta entries and peptidic sequences
    """
    # open file
    with open(file, "r", encoding="utf-8") as pep_file:
        # create peptides dictionary
        peptides = {}
        # initiate bool to tag the first iteration
        first = True
        # bool to track if protein coding
        prot_coding = False
        # loop through file
        for line in pep_file:
            # check if line is header line or end of file line
            if ">" in line or line == "\n":
                if "protein_coding" in line:
                    prot_coding = True
                    regex = re.compile(
                        r">(.+)\..+GRCh3[7-8]:(.*:\d+:\d+:\-*\d+).+gene:(.+)\..+transcript:(.+)\.\d+ gene"
                    )
                    if first:
                        pep_id, coords, gene, transcript = regex.search(line).groups()
                        first = False
                    else:
                        if pep_id not in peptides:
                            peptides[pep_id] = [gene, coords, transcript, peptide_seq]
                        else:
                            print("peptide already has a sequence")
                            sys.exit()
                        pep_id, coords, gene, transcript = regex.search(line).groups()
                    peptide_seq = ""
                else:
                    prot_coding = False
            else:
                if prot_coding:
                    peptide_seq += line[:-1]
        # ensure to add last peptide if last line of file is not "\n"
        if pep_id not in peptides and prot_coding:
            peptides[pep_id] = [gene, coords, transcript, peptide_seq]
    return peptides
