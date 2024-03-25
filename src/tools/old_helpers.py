#coding:utf-8

import sys
import numpy as np
import pandas as pd
import gzip
from pathlib import Path
import os
import re
import glob
import tools.parsing_functions as parsing_functions
import tools.ams_helpers as helpers

def convert(df_indiv,threshold):
    df_indiv["ratio"] = df_indiv["count_nt_1"].astype(int)/(df_indiv["count_nt_1"].astype(int)+df_indiv["count_nt_2"].astype(int))
    df_indiv.loc[df_indiv["ratio"]>=threshold,["count_nt_2","ratio"]] = [0,1]
    df_indiv.loc[df_indiv["ratio"]<=1-threshold,["count_nt_1","ratio"]] = [0,0]
    return(df_indiv)

# pandas.DataFrame -> pandas.DataFrame
def get_read_counts(df_indiv,indiv_file,min_DP,max_DP):
    df_indiv = df_indiv.rename(columns={"#CHROM":"CHROM"})
    df_indiv["POS"] = df_indiv["POS"].astype(int)
    df_indiv = df_indiv.drop_duplicates()
    df_indiv["CHROM"] = df_indiv["CHROM"].str.replace("chr","")
    # filter CHROM column
    df_indiv = df_indiv[(df_indiv["CHROM"].str.isdigit())|(df_indiv["CHROM"] == "X")|(df_indiv["CHROM"] == "Y")]
    indiv_path = "/".join(indiv_file.split("/")[:-1])
    indiv_pair = indiv_file.split("/")[-2]
    indiv_name = indiv_file.split("/")[-1].split(".")[0]
    if "filtered" in indiv_name:
        indiv_name = "-".join(indiv_name.split("-")[:-1])
    df_indiv["TYPE"] = df_indiv[indiv_name].str.split(pat=":").str[-1]
    
    sub_df_indiv = df_indiv[indiv_name].str.split(pat=',',expand=True)
    columns_range = list(range(2,len(list(sub_df_indiv.columns))))
    sub_df_indiv = sub_df_indiv.fillna(0)
    idx = list(sub_df_indiv[(sub_df_indiv[columns_range].T!=0).any()].index.values)
    excluded = df_indiv.loc[idx]
    # excluded.to_csv("{}/{}_{}_{}_{}_excluded_poly.tsv".format(indiv_path,indiv_pair,indiv_name,min_DP,max_DP),sep="\t")
    # print(sub_df_indiv[~(sub_df_indiv[columns_range].T==0).all()])
    sub_df_indiv = sub_df_indiv[(sub_df_indiv[columns_range].T==0).all()].drop(columns=columns_range)
    df_indiv = pd.concat([df_indiv,sub_df_indiv],join="inner",axis=1)
    df_indiv = df_indiv.fillna(value = np.nan)
    # voir ce qu'on fait pour les cas Mixture, chez AELTXQC alloval3 il y en a 75
    # insertion handle_mixture
    # handle_mixture(df_indiv)
    # voir si on conserve heterozygotie
    # insertion fonction convert
    # df_indiv = df_indiv[df_indiv[2].isnull()]
    df_indiv[0] = df_indiv[0].str.split(pat=':').str[1]
    df_indiv[1] = df_indiv[1].str.split(pat=':').str[0]
    # for i in [0,1]:
    #     df_indiv[["nt_{}".format(i+1),"count_nt_{}".format(i+1)]] = df_indiv[i].str.split(pat='=',expand=True)
    df_indiv["RO"] = np.select([df_indiv[0].str.split("=").str[0]==df_indiv["REF"],df_indiv[1].str.split("=").str[0]==df_indiv["REF"]],[df_indiv[0].str.split("=").str[1],df_indiv[1].str.split("=").str[1]],default=0)
    df_indiv["AO"] = np.select([df_indiv[0].str.split("=").str[0]==df_indiv["ALT"],df_indiv[1].str.split("=").str[0]==df_indiv["ALT"]],[df_indiv[0].str.split("=").str[1],df_indiv[1].str.split("=").str[1]],default=0)
    df_indiv = df_indiv.drop(["FORMAT",indiv_name,0,1],axis=1)
    return(df_indiv)

# pandas.DataFrame -> pandas.DataFrame
def filter_on_depth(df_indiv,min_DP,max_DP):
    df_indiv["count_nt_2"] = df_indiv["count_nt_2"].fillna(value = "0")
    df_indiv["DP"] = df_indiv["count_nt_1"].astype(int)+df_indiv["count_nt_2"].astype(int)
    df_indiv = df_indiv[df_indiv["DP"].between(min_DP,max_DP)]
    return(df_indiv)

# pandas.DataFrame,int -> pandas.DataFrame
def get_aa_indiv(df_indiv,threshold):
    # get aa ref and aa alt from VEP information filtered on read counts
    df_indiv["aa_ref_indiv"] = np.where(((df_indiv["count_nt_1"].astype(int)>=threshold)&(df_indiv["TAG_nt1"]=="REF"))|((df_indiv["count_nt_2"].astype(int)>=threshold)&(df_indiv["TAG_nt2"]=="REF")),df_indiv["aa_REF"],"")
    df_indiv["aa_alt_indiv"] = np.where(((df_indiv["count_nt_1"].astype(int)>=threshold)&(df_indiv["TAG_nt1"]=="ALT"))|((df_indiv["count_nt_2"].astype(int)>=threshold)&(df_indiv["TAG_nt2"]=="ALT")),df_indiv["aa_ALT"],"")
    # df_indiv["aa_indiv"] = df_indiv["aa_ref_indiv"] + df_indiv["aa_alt_indiv"]
    # print(df_indiv)
    df_indiv[["aa_ref_indiv","aa_alt_indiv"]] = df_indiv[["aa_ref_indiv","aa_alt_indiv"]].replace("",np.nan)
    df_indiv["aa_indiv"] = df_indiv[["aa_ref_indiv","aa_alt_indiv"]].apply(lambda x: ','.join(x.dropna().values.tolist()), axis=1)
    # print(df_indiv["aa_indiv"])
    # df_indiv["aa_indiv"] = df_indiv["aa_indiv"].apply(lambda x: list(set(x)))
    return(df_indiv)

# pandas.DataFrame, str -> pandas.DataFrame
def clean_df(df_indiv,vcf_path_indiv):
    indiv_name = vcf_path_indiv.split("/")[-1].split(".")[0]
    # remove unwanted known columns
    df_indiv = df_indiv.drop(["INFO"],axis=1)
    # remove columns containing only zeros
    df_indiv = df_indiv.loc[:, (df_indiv != 0).any(axis=0)]
    # remove rows where there is no aa in both indiv columns
    # print(df_indiv[["REF","ALT","nt_1","nt_2","aa_REF","aa_ALT","aa_ref_indiv","aa_alt_indiv"]])
    # check if useful
    df_indiv = df_indiv.dropna(subset =["aa_alt_indiv","aa_ref_indiv"],how="all")
    return(df_indiv)


def prepare_indiv_df(vcf_path_indiv,min_DP,max_DP,min_AD,homozygosity_threshold,wc,consequences_path):
    vcf_extension = vcf_path_indiv.split(".")[-1]
    if vcf_extension == "vcf":
        df_indiv,aa_vep_index,gene_vep_index = parsing_functions.vcf_vep_parser(vcf_path_indiv)
    else:
        df_indiv,aa_vep_index,gene_vep_index = parsing_functions.gzvcf_vep_parser(vcf_path_indiv)
    df_indiv = df_indiv.drop(["ID","FILTER","QUAL"],axis=1)
    df_indiv = get_read_counts(df_indiv,vcf_path_indiv,min_DP,max_DP)
    df_indiv = helpers.filter_on_depth(df_indiv,min_DP,max_DP)
    df_indiv = helpers.convert(df_indiv,homozygosity_threshold)
    df_indiv = parsing_functions.vep_infos_parser(df_indiv,aa_vep_index,gene_vep_index)
    if wc:
        vep_conseq_infos = parsing_functions.worst_consequences_parser(consequences_path)
        df_indiv = pd.merge(df_indiv,vep_conseq_infos,how="inner",on=["CHROM","POS"])
        df_indiv["aa_vcf"] = df_indiv["aa_REF"]+"/"+df_indiv["aa_ALT"]
        # print(df_indiv[(df_indiv["aa_REF"]!=df_indiv["aa_REF_vep"])|(df_indiv["aa_ALT"]!=df_indiv["aa_ALT_vep"])])
        df_indiv[["aa_REF","aa_ALT"]] = df_indiv[["aa_REF_vep","aa_ALT_vep"]]
        print("Worst consequences")
    else:
        print("All consequences")
        # df_indiv = convert(df_indiv,homozygosity_threshold)
    # update dataframe with aa ref and aa alt from VEP info
    df_indiv = helpers.get_aa_indiv(df_indiv,min_AD)
    df_indiv = helpers.clean_df(df_indiv,vcf_path_indiv)
    # print(df_indiv)
    return(df_indiv)

def merge_dfs(df_ind1,df_ind2,orientation):
    if orientation == "dr":
        merge_type = "left"
        side = "x"
        opposite = "y"
    else:
        merge_type = "right"
        side = "y"
        opposite = "x"
    # merge dataframe from 2 indivs on chosen columns, keeps rows in common
    return(pd.merge(df_ind1,df_ind2,how=merge_type,on=["CHROM","POS","REF","ALT","aa_REF","aa_ALT"]),side,opposite)

def keep_alt(merged_df,side,opposite):
    merged_df[["RO_{}".format(opposite),"AO_{}".format(opposite)]] = merged_df[["RO_{}".format(opposite),"AO_{}".format(opposite)]].fillna(0)
    merged_df = merged_df[~(((merged_df["RO_{}".format(side)]!=0)&(merged_df["AO_{}".format(side)]==0))&((merged_df["RO_{}".format(opposite)]==0)&(merged_df["AO_{}".format(opposite)]==0)))]
    # merged_df = merged_df[(merged_df["RO_{}".format(side)]==0)&((merged_df["RO_{}".format(opposite)]==0)&(merged_df))]
    return(merged_df)


def count_mismatches(merged_df,orientation):
    merged_df[["aa_indiv_x","aa_indiv_y"]] = merged_df[["aa_indiv_x","aa_indiv_y"]].fillna("")
    if orientation == "dr":
        merged_df["diff"] = [','.join(set(a.split(','))-set(b.split(','))) for a,b in zip(merged_df["aa_indiv_x"], merged_df["aa_indiv_y"])]
    if orientation == "rd":
        merged_df["diff"] = [','.join(set(b.split(','))-set(a.split(','))) for a,b in zip(merged_df["aa_indiv_x"], merged_df["aa_indiv_y"])]
    merged_df["diff"] = merged_df["diff"].replace("",np.nan)
    merged_df = merged_df.dropna(subset=["diff"]).copy()
    merged_df["mismatch"] = merged_df["diff"].str.split(",").str.len()
    mismatch = merged_df["mismatch"].astype(int).sum()
    # print(mismatches_df[mismatches_df["mismatch"]!=1])
    return(merged_df,mismatch)

def mismatch_type(merged):
    merged["mismatch_type"] = np.where(((merged["TYPE_x"]=="homozygous")&(merged["TYPE_y"]=="homozygous")),"homozygous","heterozygous")
    return(merged)

def rewrite_vcf(vcf_file):
    with open(vcf_file,'r') as read_file:
        write_file = open(".".join(vcf_file.split(".")[:-1])+"-filtered"+".vcf","w")
        for line in read_file:
            if "Mixture" in line:
                print(line)
            pat = "INDEL,Number=1"
            if pat in line:
                idx = line.index(pat)
                line = line[:pat+len(pat)-1]+"0"+line[pat+len(pat):]

            else:
                write_file.write(line)
        write_file.close()
def write_all_vcf(vcf_files):
    for vcf_file in vcf_files:
        rewrite_vcf(vcf_file)
# lst = glob.glob(os.path.dirname(os.path.realpath("../"))+"/output/indiv_vcf/cohort_GQPDOMB"+"/**/*.vcf")
# write_all_vcf(lst)

def filter_bed(cohort_path,giab_bed_path):
    files = [file for file in glob.glob(os.path.join(cohort_path,"**/*.vcf.gz"))]
    print(len(files))
    for file in files:
        path = os.path.splitext(os.path.splitext(file)[0])[0]+"_bed.vcf"
        os.system("bedtools intersect -a {} -b {} -header> {}".format(file,giab_bed_path,path))
        print(path)
        os.system("gzip {}".format(path))
