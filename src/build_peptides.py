# coding:utf-8

import os
import glob
import re
import time
from pathlib import Path
import numpy as np
import pandas as pd
import tools.parsing_functions as parsing


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
    # donor_transcripts = {key:value for key,value in ensembl_transcripts.items() if key in donor_transcripts}
    # recipient_transcripts = {key:value for key,value in ensembl_transcripts.items() if key in recipient_transcripts}
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


def substitution_insertion(x, base, position, k):
    if position - k < 0:
        pep = (
            x[0 : position - 1] + base + x[position + len(base) - 1 : position + k - 1]
        )
    elif position + k - 1 > len(x):
        pep = x[position - k : position - 1] + base + x[position + len(base) - 1 :]
    else:
        pep = (
            x[position - k : position - 1]
            + base
            + x[position + len(base) - 1 : position + k - 1]
        )
    return pep


def deletion(x, base, position, k):
    if position - k < 0:
        pep = x[0:position] + base + x[position : position + k - 1]
    elif position + k - 1 > len(x):
        pep = x[position - k + 1 : position] + base + x[position:]
    else:
        pep = x[position - k + 1 : position] + base + x[position : position + k - 1]
    return pep


def stop(x, position, k):
    if position - k < 0:
        pep = x[0 : position - 1]
    else:
        pep = x[position - k + 1 : position - 1]
    return pep


def mutation_process(x, base, position, k):
    position = int(position)
    if "-" in base:
        pep = deletion(x, base, position, k)
    elif "*" in base:
        pep = stop(x, position, k)
    else:
        pep = substitution_insertion(x, base, position, k)
    return pep


def peptide_seg(peptide, k):
    hla_peptides = [peptide[i : i + k] for i in range(len(peptide) - k + 1)]
    return hla_peptides


def get_peptides_ref(transcripts_pair, k):
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
    # frameshift = transcripts_pair[transcripts_pair["Protein_position"].str.contains("-")].copy()
    # frameshift[["b1","b2"]] = frameshift["Protein_position"].str.split("-",expand=True)
    # frameshift[["b1","b2"]] = frameshift[["b1","b2"]].astype(int)
    # frameshift["peptide"] = frameshift.apply(lambda x : mutation_process(x["Sequence_aa"],x["aa_REF"],x["b1"],k),axis=1)
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
            x["Sequence_aa"], x["aa_REF"], x["Protein_position"], k
        ),
        axis=1,
    )
    transcripts_pair["peptide"] = transcripts_pair.apply(
        lambda x: mutation_process(
            x["Sequence_aa"], x["diff"], x["Protein_position"], k
        ),
        axis=1,
    )
    print(transcripts_pair[["CHROM", "POS", "aa_REF", "diff"]])
    print(transcripts_pair.dtypes)
    transcripts_pair = transcripts_pair.drop_duplicates()
    # frameshift = frameshift.drop_duplicates()
    # frameshift = frameshift.drop(["b1","b2"],axis=1)
    # transcripts_pair = transcripts_pair.append(frameshift,ignore_index=True)
    transcripts_reduced = transcripts_pair.groupby(
        ["CHROM", "POS", "Gene_id", "peptide_REF", "peptide"], as_index=False
    )[["Transcript_id", "Peptide_id"]].agg("first")
    transcripts_reduced["hla_peptides_REF"] = transcripts_reduced["peptide_REF"].apply(
        lambda x: peptide_seg(x, k)
    )
    transcripts_reduced["hla_peptides"] = transcripts_reduced["peptide"].apply(
        lambda x: peptide_seg(x, k)
    )

    return transcripts_reduced


# def select_different_kmers(donor_reduced,recipient_reduced,orientation):
# 	test = pd.merge(transcripts_reduced_donor,transcripts_reduced_recipient,how="outer",on=["CHROM","POS","Gene_id","peptide","Transcript_id","Peptide_id"])
# 	if orientation == "dr":
# 		test["different_kmers"] = [','.join(set(a.split(','))-set(b.split(','))) for a,b in zip(test["kmers_x"], test["kmers_y"])]
# 	else:
# 		test["different_kmers"] = [','.join(set(b.split(','))-set(a.split(','))) for a,b in zip(test["kmers_x"], test["kmers_y"])]

# 	print(test)
# 	return(different_kmers)


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
        with open(
            os.path.join(dir_path, netmhc_fasta_file_name), "w", encoding="utf-8"
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


#############################################################
##################### OLD SAVING SYSTEM #####################
#############################################################
# def main():
# 	start = time.time()
# 	directories = [os.path.join("../output/indiv_vcf/joint_genotyping/hard-filtered",x) for x in os.listdir("../output/indiv_vcf/joint_genotyping/hard-filtered") if os.path.isdir(os.path.join("../output/indiv_vcf/joint_genotyping/hard-filtered",x)) and re.match("R\d+",x)]
# 	print(directories)
# 	# done = ['R908','R383','R1048','R6764']
# 	print(len(directories))
# 	for direct in directories:
# 		if direct.split("/")[-1] not in done:
# 			# get cdna transcripts from ensembl database
# 			ensembl_transcripts = parsing.read_fasta("../output/Ensembl/Homo_sapiens.GRCh38.cdna.all.103.fa")
# 			# get cds transcripts from ensembl database
# 			ensembl_transcripts_cds = parsing.read_fasta("../output/Ensembl/Homo_sapiens.GRCh38.cds.all.fa")
# 			print("Ensembl transcripts in dictionary : {}".format(len(ensembl_transcripts.keys())))
# 			# get proteins from ensembl database
# 			proteins = parsing.read_pep_fa("../output/Ensembl/Homo_sapiens.GRCh38.pep.all.103.fa")
# 			print("Ensembl proteins in dictionary : {}".format(len(proteins.keys())))
# 			peptides_ensembl = dict_to_df(proteins)
# 			print("{} peptides selected in ensembl refseq".format(len(peptides_ensembl)))
# 			# print(peptides_ensembl)

# 			# get merged file containing ams positions
# 			merged_pair = pd.read_csv(os.path.join(direct,"bed_merged_df_20_400_5_0.8_20_TRANSCRIPTS_GENES_STOPLOST.tsv"),sep="\t")
# 			# get list of transcripts in AMS and intersect it with ensembl transcripts
# 			ams_transcripts = contributing_ams_transcripts(merged_pair,ensembl_transcripts)
# 			print("Potential contributing transcripts after ensembl filtering : {}".format(len(ams_transcripts)))
# 			refseq_transcripts_path = "../output/Ensembl/Homo_sapiens.GRCh38.103.refseq.tsv"
# 			# filtering to keep transcripts present in refseq table
# 			transcripts_df = filter_on_refseq(ams_transcripts,refseq_transcripts_path)
# 			print("Transcripts after REFSEQ filter : {}".format(len(transcripts_df)))
# 			# transcripts long format
# 			transcripts_pair = pd.read_csv(os.path.join(direct,"transcripts_pair_codons_20_400_5_0.8_20_STOPLOST.tsv"),sep="\t")
# 			# intersect filtered transcripts and long format to get a long format with all info
# 			transcripts_pair = intersect_positions(transcripts_pair,transcripts_df)
# 			transcripts_pair = aa_REF(transcripts_pair)

# 			# get pep seq on long format table
# 			transcripts_pair = add_pep_seq(transcripts_pair,peptides_ensembl)
# 			transcripts_pair.to_csv("ugly.tsv",sep="\t",index = False)
# 			k=9
# 			transcripts_reduced = get_peptides_ref(transcripts_pair,k)
# 			orientation = "rd"
# 			# print(transcripts_reduced)
# 			transcripts_reduced.to_pickle(os.path.join(direct,"pep_df_20_400_5_0.8_20_TRANSCRIPTS_GENES_STOPLOST.pkl"))
# 			write_pep_fasta(os.path.join(direct,"test_kmers.fa"),transcripts_reduced)
# 			write_netmhc_fasta(os.path.join(direct,"test_kmers.fa"),"netmhc_fasta.fa")
# 	end = time.time()
# 	print((end-start)/60)


# considered inputs : run name
def main(run_name,ens_transcripts_path,ens_transcripts_cds_path,ens_pep_path,refseq_transcripts_path,pep_size):
    # start = time.time()
    # directories = [
    #     os.path.join("../output/indiv_vcf/joint_genotyping/hard-filtered", x)
    #     for x in os.listdir("../output/indiv_vcf/joint_genotyping/hard-filtered")
    #     if os.path.isdir(
    #         os.path.join("../output/indiv_vcf/joint_genotyping/hard-filtered", x)
    #     )
    #     and re.match(r"R\d+", x)
    # ]
    # directories = sorted(
    #     directories, key=lambda x: int("".join(re.split(Path("R").name, x)[1]))
    # )
    # print(directories)
    # # done = ['R908','R383','R1048','R6764']
    # print(len(directories))
    run_tables = f"../output/runs/{run_name}/run_tables"
    aams_run_tables = f"../output/runs/{run_name}/aams_run_tables"
    merged_files = glob.glob(f"{run_tables}/*_bed_merged*")
    transcripts_files = glob.glob(f"{run_tables}/*transcripts*")
    merged_files = sorted(
        merged_files,
        key=lambda x: int(
            "".join(re.split("R|P", str(Path(x).name).split("_", maxsplit=1)[0])[1])
        ),
    )
    transcripts_files = sorted(
        transcripts_files,
        key=lambda x: int(
            "".join(re.split("R|P", str(Path(x).name).split("_", maxsplit=1)[0])[1])
        ),
    )
    # get cdna transcripts from ensembl database
    ens_transcripts = parsing.read_fasta(ens_transcripts_path)
    # get cds transcripts from ensembl database
    ens_transcripts_cds = parsing.read_fasta(ens_transcripts_cds_path)
    print(f"Ensembl transcripts in dictionary : {len(ens_transcripts.keys())}")
    # get proteins from ensembl database
    proteins = parsing.read_pep_fa(ens_pep_path)
    print(f"Ensembl proteins in dictionary : {len(proteins.keys())}")
    peptides_ensembl = dict_to_df(proteins)
    print(f"{len(peptides_ensembl)} peptides selected in ensembl refseq")
    for i, merged in enumerate(merged_files):
        pair = str(Path(merged).name).split("_", maxsplit=1)[0]
        print(pair)
        # get merged file containing ams positions
        merged_pair = pd.read_csv(merged, sep="\t")
        # get list of transcripts in AMS and intersect it with ensembl transcripts
        ams_transcripts = contributing_ams_transcripts(merged_pair, ens_transcripts)
        print(
            f"Potential contributing transcripts after ensembl filtering : {len(ams_transcripts)}"
        )
        # filtering to keep transcripts present in refseq table
        transcripts_df = filter_on_refseq(ams_transcripts, refseq_transcripts_path)
        print(f"Transcripts after REFSEQ filter : {len(transcripts_df)}")
        # transcripts long format
        transcripts_pair = pd.read_csv(transcripts_files[i], sep="\t")
        # intersect filtered transcripts and long format to get a long format with all info
        transcripts_pair = intersect_positions(transcripts_pair, transcripts_df)
        transcripts_pair = aa_ref(transcripts_pair)
        # get pep seq on long format table
        transcripts_pair = add_pep_seq(transcripts_pair, peptides_ensembl)
        transcripts_pair.to_csv(os.path.join(aams_run_tables,f"{pair}_{run_name}_full.tsv"), sep="\t", index=False)
        transcripts_reduced = get_peptides_ref(transcripts_pair, pep_size)
        # orientation = "rd"
        # print(transcripts_reduced)
        transcripts_reduced.to_pickle(
            os.path.join(aams_run_tables, f"{pair}_{run_name}_pep_df_20_400_5_gq_20_0.8.pkl")
        )
        write_pep_fasta(
            os.path.join(aams_run_tables, f"{pair}_{run_name}_kmers.fa"), transcripts_reduced
        )
        write_netmhc_fasta(
            os.path.join(aams_run_tables, f"{pair}_{run_name}_kmers.fa"),
            f"{pair}_{run_name}_netmhc_fasta.fa",
        )

    return 0

    # 		transcripts_reduced.to_pickle(os.path.join(aams_run_tables,"pep_df_20_400_5_0.8_20_TRANSCRIPTS_GENES_STOPLOST.pkl"))
    # 		# write in AAMS_run folder ?
    # 		write_pep_fasta(os.path.join(direct,"test_kmers.fa"),transcripts_reduced)
    # 		write_netmhc_fasta(os.path.join(direct,"test_kmers.fa"),"netmhc_fasta.fa")

    # for direct in directories:
    # 		# print(peptides_ensembl)


if __name__ == "__main__":
    main()
