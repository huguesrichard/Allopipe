#coding:utf-8

import sys
import numpy as np
import pandas as pd
import gzip
import re
import os
from pathlib import Path
import tools.parsing_functions as parsing_functions

# python3 ams_pipeline.py ../output/Adaptive\ Sampling/5000_AS/pepper_r2.vcf.gz ../output/Adaptive\ Sampling/5000_AS/Hg38_5000_intersect_MS_VEP_illumina.vcf 0 400 0 0.95 0 dr -n AS_5000 -l 3 -f

def create_run_directory(run_name):
	Path("../output/runs").mkdir(parents=True,exist_ok=True)
	run_path = os.path.join("../output/runs",run_name)
	run_tables = os.path.join(run_path,"run_tables")
	run_plots = os.path.join(run_path,"plots")
	run_beds = os.path.join(run_path,"bedfiles")
	run_ams = os.path.join(run_path,"AMS_folder")
	Path(run_path).mkdir(parents=True, exist_ok=True)
	Path(run_tables).mkdir(parents=True,exist_ok=True)
	Path(run_plots).mkdir(parents=True,exist_ok=True)
	Path(run_ams).mkdir(parents=True,exist_ok=True)
	return(run_path,run_tables,run_plots,run_ams)

#############
# giab v1.3 #
#############

def fix_AD_field_for_np(dp,ad,gt):
	if "," not in ad:
		x = int(dp)-int(ad)
		if gt=="0/0" or gt=="0/1":
			new = "{},{}".format(ad,x)
		elif gt=="1/1":
			new = "{},{}".format(x,ad)
		else:
			print(gt,ad)
		return(new)
	else:
		return(ad)

# allow genotypes defined in multi sample vcf to be taken into account (other than 0/0 0/1 1/1)
def explode_GT_AD(df_indiv):
	# split AD field
	subset = df_indiv[["CHROM","POS","REF","ALT","GT","AD"]].copy()
	indiv_name = list(df_indiv.columns)[-1]
	if indiv_name=="Sample": 
		subset = df_indiv[["CHROM","POS","REF","ALT","GT","DP","AD"]].copy()
		subset["AD"] = subset.apply(lambda x: fix_AD_field_for_np(x["DP"],x["AD"],x["GT"]),axis=1)
	print(subset)
	subset["AD"] = subset["AD"].apply(lambda x:x.split(","))
	# explode dataframe -> duplicate lines except for the AD field which will be split on the lines
	subset = subset.explode("AD")
	# change type of AD field from object to int
	check_excluded = subset[subset["AD"]=="."]
	subset = subset[subset["AD"] != "."].copy()
	subset["AD"] = subset["AD"].astype(int)
	subset["ind"] = list(subset.index)
	# index counting the number of values in the AD field, ex : 36,2,5 -> 0 1 2 
	subset["ind"] = subset.groupby("ind").cumcount()
	subset = subset.reset_index()
	# proceed exactly the same for GT
	subset["GT"] = subset["GT"].apply(lambda x:x.split("/"))
	subset = subset.explode("GT",ignore_index=True)
	subset["GT"] = subset["GT"].astype(int)
	# keep only the right GT fields, ex : GT=0/2, AD=36,2,5 -> keep 0 = 0 and 2 = 2 (36 and 5 will be kept)
	subset = subset[subset["GT"]==subset["ind"]]
	# remove duplicate rows
	subset = subset.drop_duplicates()
	return(subset)

def parse_format(line):
	gt = line.index("GT")
	gq = line.index("GQ")
	ad = line.index("AD")
	dp = line.index("DP")
	# ft = line.index("FT")
	return([gt,gq,ad,dp])

def apply_format(line,indices):
	ch = ""
	for i in range(len(indices)):
		if i == len(indices)-1:
			ch+= line[indices[i]]
		else:
			ch+= line[indices[i]]+"_"
	return(ch)

# pandas.DataFrame -> pandas.DataFrame
def get_read_counts(df_indiv,indiv_file,min_AD,min_GQ,full=False):
	# rename the #CHROM column to CHROM
	df_indiv = df_indiv.rename(columns={"#CHROM":"CHROM"})
	# switch type of POS column to int
	df_indiv["POS"] = df_indiv["POS"].astype(int)
	df_indiv = df_indiv.drop_duplicates()
	# remove chr in front of the chromosome number
	df_indiv["CHROM"] = df_indiv["CHROM"].str.replace("chr","")
	###############
	# legacy v1.2 #
	###############
	# filter CHROM column
	df_indiv = df_indiv[(df_indiv["CHROM"].str.isdigit())|(df_indiv["CHROM"] == "X")|(df_indiv["CHROM"] == "Y")]
	# get the name of the indiv
	indiv_path = "/".join(indiv_file.split("/")[:-1])
	indiv_pair = indiv_file.split("/")[-2]
	# indiv_name = indiv_file.split("/")[-1].split(".")[0]
	indiv_name = list(df_indiv.columns)[-1]
	# indiv_name = "Sample"
	print(indiv_name)
	print(df_indiv["FORMAT"].unique().tolist())
	if "bed" in indiv_name or "VIP" in indiv_name:
		indiv_name = indiv_name.split("_")[0]
	format_infos = df_indiv["FORMAT"].tolist()
	# get the max columns of format field
	format_infos_max = format_infos[format_infos.index(max(format_infos,key=len))].split(":")
	format_length = len(format_infos_max)
	# split and expand the indiv field thanks to the format field column names
	df_indiv["FORMAT"] = df_indiv["FORMAT"].str.split(":")
	df_indiv["select"] = df_indiv[indiv_name].str.split(":")
	df_indiv["split_indices"] = df_indiv["FORMAT"].apply(lambda x:parse_format(x))
	df_indiv["select"] = df_indiv.apply(lambda x:apply_format(x["select"],x["split_indices"]),axis=1)
	df_indiv[["GT","GQ","AD","DP"]] = df_indiv["select"].str.split("_",expand=True)
	# df_indiv[format_infos_max] = df_indiv[indiv_name].str.split(pat=":",expand=True,n=format_length)
	df_indiv = df_indiv[~df_indiv["GT"].str.contains(r"\.")]

	# some positions are annotated 0 or 1 and seem to me 0/0 and 1/1
	# df_indiv["GT"] = df_indiv["GT"].replace("0","0/0")
	# df_indiv["GT"] = df_indiv["GT"].replace("1","1/1")
	df_indiv = df_indiv.dropna(subset =["GT"])
	# Genotype quality filtered above 20, defined value
	#############
	# giab v1.3 #
	# #############
	# df_indiv = df_indiv[df_indiv["GQ"].astype(float)>min_GQ]
	# Filter FT field
	###############
	# legacy v1.2 #
	###############
	# df_indiv = df_indiv[df_indiv["FT"]=="PASS"]
	# filter REF calls with higher length than 3, value should be selected in cmd line instead of written manually here
	###############
	# legacy v1.2 #
	###############
	df_indiv = df_indiv[df_indiv["REF"].str.len() <= 3]
	#### ADD CSV EXCLUDED
	###############
	# legacy v1.2 #
	###############
	df_indiv["phased"] = np.where((df_indiv["GT"].str.contains("|")),df_indiv["GT"],np.nan)
	# some GT are phased, remove this info in GT field
	df_indiv["GT"] = df_indiv["GT"].str.replace("|","/",regex=True)
	# when not phased, 1/0 is equivalent to 0/1
	df_indiv["GT"] = df_indiv["GT"].str.replace("1/0","0/1",regex=True)
	# remove DP field if present
	# if "DP" in df_indiv.columns:
	# 	df_indiv = df_indiv.drop("DP",axis=1)
	# explode df_indiv to get only the right AD values for the given GT
	# ###################################
	subset = explode_GT_AD(df_indiv)
	print("1st",df_indiv[df_indiv["POS"]==140247143])
	# filter AD
	###############
	# legacy v1.1 #
	###############
	subset["AD"] = subset["AD"].astype(int)
	subset = subset[subset["AD"]>=min_AD]
	# filter ALT longer than threshold (here 3), should also be selected in cmd line
	subset["ALT_length"] = np.where(subset["GT"]==0,0,subset.apply(lambda x: len(x["ALT"].split(",")[x["GT"]-1]),axis=1))
	#############################################
	sub = subset.groupby(["CHROM","POS"],as_index=False)["AD"].sum()
	sub = sub.rename({"AD":"sumAD"},axis=1)
	###############
	# legacy v1.2 #
	###############
	###################################
	subset = subset[subset["ALT_length"]<=3]
	# put the DP column for the exploded subset
	###################################
	subset = pd.merge(subset,sub,how="outer",on=["CHROM","POS"])
	subset = subset[(subset[["CHROM","POS"]].duplicated(keep=False))|(subset["AD"]==subset["sumAD"])]
	subset = subset.drop("sumAD",axis=1)
	sub = subset.groupby(["CHROM","POS"],as_index=False)["AD"].sum()
	sub = sub.rename({"AD":"DP"},axis=1)
	subset = pd.merge(subset,sub,how="outer",on=["CHROM","POS"])
	df_indiv = pd.merge(df_indiv,sub,how="outer",on=["CHROM","POS"])
	print(df_indiv)
	df_indiv = df_indiv.dropna(subset=["DP_x"])
	# remove
	df_indiv = df_indiv.drop("DP_y",axis=1)
	df_indiv = df_indiv.rename({"DP_x":"DP"},axis=1)
	if indiv_name == "Sample":
		subset = subset.drop("DP_y",axis=1)
		subset = subset.rename({"DP_x":"DP"},axis=1)
	print("reads counts complete")
	return(df_indiv,subset)

# pandas.DataFrame,int -> pandas.DataFrame
def filter_on_depth(df_indiv,subset,min_DP,max_DP):
	# filter df_indiv and subset on DP
	###############
	# legacy v1.1 #
	###############
	df_indiv["DP"] = df_indiv["DP"].astype(int)
	subset["DP"] = subset["DP"].astype(int)
	df_indiv = df_indiv[df_indiv["DP"].between(min_DP,max_DP)]
	subset = subset[subset["DP"].between(min_DP,max_DP)]
	print("filtering on depth complete")
	return(df_indiv,subset)

def convert(df_indiv,subset,homozygosity_threshold):
	###############
	# legacy v1.1 #
	###############
	# create TYPE field
	df_indiv["TYPE"] = np.where(df_indiv["GT"].str.split("/").str[0]==df_indiv["GT"].str.split("/").str[1],"homozygous","heterozygous")
	subset = pd.merge(subset,df_indiv[["CHROM","POS","DP","TYPE"]],how="outer",on=["CHROM","POS","DP"])
	# compute allelic ratio 
	subset["ratio"] = subset["AD"].astype(int)/subset["DP"]
	# keep only rows where the ratio is higher than the threshold
	subset = subset[subset["ratio"]>=1-homozygosity_threshold]
	# change the TYPE column of the rows with ratio < threshold to homozygous
	subset.loc[~(subset[["CHROM","POS"]].duplicated(keep=False))&(subset["ratio"].between(0.8,1,inclusive="left")),"TYPE"]="homozygous"
	subset.loc[(subset["ratio"]==1)&(subset["TYPE"]=="heterozygous"),"TYPE"] = "homozygous"
	type_subset = subset.groupby(["CHROM","POS"],as_index=False)["TYPE"].first()
	print(df_indiv[df_indiv.duplicated(subset=["CHROM","POS"],keep=False)])
	# WIP
	# if len(type_subset)==len(df_indiv):
	df_indiv["TYPE"]=type_subset["TYPE"]
	# else:
		# print("different length")
		# print(len(type_subset),len(df_indiv))
	# df_indiv = pd.merge(df_indiv,type_subset,how="outer",on=["CHROM","POS"])
	print("converting ambiguous heterozygous positions complete")
	return(df_indiv)

# pandas.DataFrame,int -> pandas.DataFrame
def get_aa_indiv(df_indiv,threshold):
	# get aa ref and aa alt from VEP information filtered on read counts
	df_indiv["aa_ref_indiv"] = np.where(df_indiv["GT"].str.contains("0"),df_indiv["aa_REF"],"")
	df_indiv["aa_alt_indiv"] = np.where((~(df_indiv["GT"].str.contains("0"))|(df_indiv["GT"].str.split("/").str[0]!=df_indiv["GT"].str.split("/").str[1])),df_indiv["aa_ALT"],"")
	df_indiv[["aa_ref_indiv","aa_alt_indiv"]] = df_indiv[["aa_ref_indiv","aa_alt_indiv"]].replace("",np.nan)
	df_indiv["aa_indiv"] = df_indiv[["aa_ref_indiv","aa_alt_indiv"]].apply(lambda x: ','.join(x.dropna().values.tolist()), axis=1)
	# print(df_indiv["aa_indiv"])
	# df_indiv["aa_indiv"] = df_indiv["aa_indiv"].apply(lambda x: list(set(x)))
	return(df_indiv)


# pandas.DataFrame, str -> pandas.DataFrame
def clean_df(df_indiv,vcf_path_indiv):
	indiv_name = vcf_path_indiv.split("/")[-1].split(".")[0]
	# remove unwanted known columns
	ls_drop = ["INFO"]
	if indiv_name in df_indiv.columns:
		ls_drop.append(indiv_name)
	df_indiv = df_indiv.drop(ls_drop,axis=1)
	# remove columns containing only zeros
	df_indiv = df_indiv.loc[:, (df_indiv != 0).any(axis=0)]
	# remove rows where there is no aa in both indiv columns
	# check if useful
	df_indiv = df_indiv.dropna(subset =["aa_alt_indiv","aa_ref_indiv"],how="all")
	return(df_indiv)

# build DataFrame from vcf file
# str,int,int -> pandas.DataFrame
def prepare_indiv_df(run_tables, vcf_path_indiv, args, consequences_path):
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
        run_tables, df_indiv, vep_indices, vcf_path_indiv, args
    )
    if args.wc:
        vep_conseq_infos = parsing_functions.worst_consequences_parser(
            consequences_path
        )
        df_indiv = pd.merge(
            df_indiv, vep_conseq_infos, how="inner", on=["CHROM", "POS"]
        )
        df_indiv["aa_vcf"] = df_indiv["aa_REF"] + "/" + df_indiv["aa_ALT"]
        df_indiv[["aa_REF", "aa_ALT"]] = df_indiv[["aa_REF_vep", "aa_ALT_vep"]]
        print("Worst consequences")
    else:
        # remove above vep_infos_parser ?
        # parse all consequences VEP information and add to dataframe
        # df_indiv = parsing_functions.vep_infos_parser(df_indiv,aa_vep_index)
        print("All consequences")
    # update dataframe with aa ref and aa alt from VEP info
    df_indiv = get_aa_indiv(df_indiv)
    # remove unnecessary columns
    df_indiv = clean_df(df_indiv, vcf_path_indiv)
    return df_indiv, vep_table_indiv, vep_indices




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
	###############
	# legacy v1.2 #
	###############
	# remove positions where REF/REF and NaN (we assume it is also REF/REF)
	merged_df = merged_df[~((merged_df["GT_{}".format(side)]=="0/0")&(merged_df["GT_{}".format(opposite)].isna()))]
	return(merged_df)


def count_mismatches(merged_df,orientation):
	###############
	# legacy v1.2 #
	###############
	# fill nan with REF for individuals
	merged_df.loc[merged_df["GT_x"].isna(),"aa_indiv_x"] = merged_df["aa_REF"]
	merged_df.loc[merged_df["GT_y"].isna(),"aa_indiv_y"] = merged_df["aa_REF"]
	merged_df[["aa_indiv_x","aa_indiv_y"]] = merged_df[["aa_indiv_x","aa_indiv_y"]].fillna("")
	# create a diff column, oriented based on the cmd line argument (kidney : dr, hsc : rd)
	if orientation == "dr":
		merged_df["diff"] = [','.join(set(a.split(','))-set(b.split(','))) for a,b in zip(merged_df["aa_indiv_x"], merged_df["aa_indiv_y"])]
	if orientation == "rd":
		merged_df["diff"] = [','.join(set(b.split(','))-set(a.split(','))) for a,b in zip(merged_df["aa_indiv_x"], merged_df["aa_indiv_y"])]
	merged_df["diff"] = merged_df["diff"].replace("",np.nan)
	merged_df = merged_df.dropna(subset=["diff"]).copy()
	merged_df["mismatch"] = merged_df["diff"].str.split(",").str.len()
	# create a mismatch column counting the aminoacids in the diff column
	mismatches_df = merged_df.dropna(subset=["diff"])
	mismatches_df["mismatch"] = mismatches_df["diff"].str.split(",").str.len()
	mismatch = mismatches_df["mismatch"].astype(int).sum()
	return(merged_df,mismatch)


def mismatch_type(merged):
	###############
	# legacy v1.2 #
	###############
	# tag each mismatch in merged dataframe with homozygous or heterozygous
	merged["mismatch_type"] = np.where(((merged["TYPE_x"]=="homozygous")&(merged["TYPE_y"]=="homozygous")),"homozygous","heterozygous")
	return(merged)

# WIP
def get_snps(run_path,donor_path,recipient_path,min_AD):
	df_d = pd.read_csv(donor_path,sep="\t")
	df_r = pd.read_csv(recipient_path,sep="\t")
	donor = Path(donor_path).stem.split("_")[0]
	recipient = Path(recipient_path).stem.split("_")[0]
	df_d["GT"] = df_d["GT"].replace("1/0","0/1")
	df_r["GT"] = df_r["GT"].replace("1/0","0/1")
	# remove 0/0 genotype (REF/REF) unless they are tagged as heterozygous
	# homozygous sometimes NaN ?
	##########################################################################
	df_d = df_d[~((df_d["GT"]=="0/0")&(df_d["TYPE"]=="homozygous"))]
	df_r = df_r[~((df_r["GT"]=="0/0")&(df_r["TYPE"]=="homozygous"))]
	# save 
	# df_d.to_csv(os.path.join(run_path,"SNPs_{}_minAD_{}.tsv").format(donor,minAD),sep="\t",index=False)
	# df_r.to_csv(os.path.join(run_path,"SNPs_{}_minAD_{}.tsv").format(donor,minAD),sep="\t",index=False)
	# merge both df 
	merged_pos = pd.merge(df_d,df_r,how="outer",on=["CHROM","POS","REF"])
	out = merged_pos[(merged_pos["ALT_x"].isna())|(merged_pos["ALT_y"].isna())]
	common = merged_pos[~((merged_pos["ALT_x"].isna())|(merged_pos["ALT_y"].isna()))][["CHROM","POS","REF","ALT_x","ALT_y","GT_x","GT_y","TYPE_x","TYPE_y","AD_x","AD_y"]]
	# snps = merged_pos[((merged_pos["GT_x"]!="0/0")|(merged_pos["AO_x"]>min_AD))&((merged_pos["GT_y"]!="0/0")|(merged_pos["AO_y"]>min_AD))]
	common_snps = snps[snps["GT_x"]==snps["GT_y"]]


# filter based on fisher test results and binomial
def create_noisy_ams_regions_bed(intervals_table,window,alpha):
	# filter based on number of pairs and p-values
	intervals_table = intervals_table[~((intervals_table["npairs"]<=15)|(intervals_table["pval_GVHa"]<alpha)|(intervals_table["pval_GVHc"]<alpha)|(intervals_table["pval_GVH"]<alpha)|(intervals_table["pval_noGVH"]<alpha))]
	# split the interval column
	intervals_table[["CHROM","pos_start"]] = intervals_table["interval"].str.split("-",expand=True)
	# convert to int type
	intervals_table["pos_start"] = intervals_table["pos_start"].astype(int)
	# get the accurate position numbers
	intervals_table["pos_start"] = intervals_table["pos_start"]*window
	intervals_table["pos_end"] = intervals_table["pos_start"]+window
	# sort ascend
	intervals_table["CHR"] = np.select([intervals_table["CHROM"].str.isdigit(),intervals_table["CHROM"]=="X",intervals_table["CHROM"]=="Y"],[intervals_table["CHROM"],"23","24"])
	intervals_table["CHR"] = intervals_table["CHR"].astype(int)
	intervals_table = intervals_table[["CHROM","pos_start","pos_end","CHR"]].sort_values(by=["CHR","pos_start"])
	# add chr
	intervals_table["CHROM"] = "chr"+intervals_table["CHROM"]
	# drop CHR column
	intervals_table = intervals_table.drop("CHR",axis=1)
	# save in bed format
	intervals_table.to_csv("noisy_regions.bed",sep="\t",index=False,header=False)


