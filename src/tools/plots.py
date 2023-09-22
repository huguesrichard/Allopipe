#coding:utf-8

import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import sys
import glob
import os
import table_operations as table_ops
import re
import seaborn as sns
from pathlib import Path
import matplotlib.ticker as ticker
import time
from datetime import date

plt.style.use("ggplot")

def plot_AMS(AMS_file,cohort,save=False):
	df = pd.read_excel(AMS_file).replace(0,np.nan)
	print(df)
	df = df.dropna(subset=["AMS_mesnard","AMS_oriented_CHR"])
	print(df)
	x,y = df["AMS_mesnard"].values.reshape(-1,1),df["AMS_oriented_CHR"].values.reshape(-1,1)
	reg = LinearRegression().fit(x,y)
	y_pred = reg.predict(x)
	plt.scatter(x,y,label="AMS pairs",c="c")
	plt.plot(x,y_pred,label="Linear Regression")
	plt.xlabel("AMS Mesnard et al.")
	plt.ylabel("AMS 2021")
	plt.title("Comparison of AMS")
	print(reg.score(x,y),reg.coef_)
	# plt.text(750,1250,"a = {}".format(reg.coeff_[0][0].round(4)))
	plt.legend()
	if save :
		plt.savefig("../../output/plots/AMS_comparison_{}.png".format(cohort))
	plt.show()

def plot_eGFR_corr(AMS_file,eGFR_file,cohort,save=False):
	"""
	# voir gestion np.nan
	# df1 = pd.read_excel(AMS_file)
	# df2 = pd.read_excel(eGFR_file)
	# merged_df = pd.merge(df1,df2,how="inner",on=["pair"])
	# print(merged_df)
	# x,y,z,t = merged_df["AMS_oriented_CHR"].values.reshape(-1,1),merged_df["MDRD eGFR 12 months"].values.reshape(-1,1),merged_df["MDRD eGFR 24 months"].values.reshape(-1,1),merged_df["MDRD eGFR 36 months"].values.reshape(-1,1)
	x,y,z,t = df1["AMS_oriented_CHR"].values.reshape(-1,1),df1["MDRD eGFR 12 months"].values.reshape(-1,1),df1["MDRD eGFR 24 months"].values.reshape(-1,1),df1["MDRD eGFR 36 months"].values.reshape(-1,1)
	regxy = LinearRegression().fit(x,y)
	regxz = LinearRegression().fit(x,z)
	regxt = LinearRegression().fit(x,t)
	print(regxy.score(x,y),regxy.coef_[0][0])
	print(regxz.score(x,z),regxz.coef_[0][0])
	print(regxt.score(x,t),regxt.coef_[0][0])
	y_pred = regxy.predict(x)
	z_pred = regxz.predict(x)
	t_pred = regxt.predict(x)
	plt.scatter(x,t,c="c")
	plt.plot(x,t_pred,label="Linear Regression")
	plt.xlabel("AMS 2021")
	plt.ylabel("eGFR 36 months")
	# plt.text(1250,20,"a = {}".format(regxy.coef_[0][0].round(4)),verticalalignment="bottom")
	plt.title("Correlation of AMS and eGFR 36 months")
	plt.legend()
	if save:
		plt.savefig("../../output/plots/corr_ams_36mo_egfr_{}.png".format(cohort))
	plt.show()
	"""
	df = pd.read_excel(AMS_file)
	subset1 = df[["pair","AMS_oriented_CHR","MDRD eGFR 12 months"]].dropna(subset=["MDRD eGFR 12 months"])
	subset2 = df[["pair","AMS_oriented_CHR","MDRD eGFR 24 months"]].dropna(subset=["MDRD eGFR 24 months"])
	subset3 = df[["pair","AMS_oriented_CHR","MDRD eGFR 36 months"]].dropna(subset=["MDRD eGFR 36 months"])
	x1,y = subset1["AMS_oriented_CHR"].values.reshape(-1,1),subset1["MDRD eGFR 12 months"].values.reshape(-1,1)
	x2,z = subset2["AMS_oriented_CHR"].values.reshape(-1,1),subset2["MDRD eGFR 24 months"].values.reshape(-1,1)
	x3,t = subset3["AMS_oriented_CHR"].values.reshape(-1,1),subset3["MDRD eGFR 36 months"].values.reshape(-1,1)
	regxy = LinearRegression().fit(x1,y)
	regxz = LinearRegression().fit(x2,z)
	regxt = LinearRegression().fit(x3,t)
	print(regxy.score(x1,y),regxy.coef_[0][0])
	print(regxz.score(x2,z),regxz.coef_[0][0])
	print(regxt.score(x3,t),regxt.coef_[0][0])
	y_pred = regxy.predict(x1)
	z_pred = regxz.predict(x2)
	t_pred = regxt.predict(x3)
	plt.scatter(x1,y,c="c")
	plt.plot(x1,y_pred,label="Linear Regression")
	plt.xlabel("AMS 2021")
	plt.ylabel("eGFR 12 months")
	# plt.text(1250,20,"a = {}".format(regxy.coef_[0][0].round(4)),verticalalignment="bottom")
	plt.title("Correlation of AMS and eGFR 12 months")
	plt.legend()
	if save:
		plt.savefig("../../output/plots/reg_all_ams_12mo_egfr_{}.png".format(cohort))
	plt.show()

def plot_acute(merged_df,save):
	merged_df["colors"] = np.where(merged_df["acute GVH"]==1,"r","b")
	merged_df["acute_gvh_status"] = np.where(merged_df["acute GVH"]==1,"yes","no")
	merged_df["AMS_all_conseq"].plot(kind="bar",stacked=True,color=merged_df["colors"])
	plt.xlabel("PAIRS")
	plt.ylabel("AMS_value")
	plt.title("AMS values for marrow data")
	labels = merged_df["acute_gvh_status"].unique()
	colors = {"yes":"r","no":"b"}
	handles = [plt.Rectangle((0,0),1,1,color=colors[lab])for lab in labels]
	plt.legend(handles,labels,title="acute GVH")
	# plt.legend()
	if save:
		plt.savefig("../../output/plots/acute_gvh.png")
	plt.show()

##### duplicated code, remove
def plot_chronic(merged_df,save):
	merged_df["colors"] = np.where(merged_df["GVH chronique"]==1,"r","b")
	merged_df["chronic_gvh_status"] = np.where(merged_df["GVH chronique"]==1,"yes","no")
	merged_df["AMS_all_conseq"].plot(kind="bar",stacked=True,color=merged_df["colors"])
	plt.xlabel("PAIRS")
	plt.ylabel("AMS_value")
	plt.title("AMS values for marrow data")
	labels = merged_df["chronic_gvh_status"].unique()
	colors = {"yes":"r","no":"b"}
	handles = [plt.Rectangle((0,0),1,1,color=colors[lab])for lab in labels]
	plt.legend(handles,labels,title="GVH chronique")
	# plt.legend()
	if save:
		plt.savefig("../../output/plots/chronic_gvh.png")
	plt.show()

# def plot_gvh_grades(merged_df,column,save):
# 	merged_df = merged_df[merged_df["acute GVH"]==1]
# 	grade_subset = merged_df[column].dropna()
# 	grade_subset.plot(kind="bar",stacked=True)
# 	plt.xlabel("PAIRS")
# 	plt.ylabel("AMS_value")
# 	plt.title("AMS values for marrow data")
# 	handles = [plt.Rectangle((0,0),1,1,color=colors[lab])for lab in labels]
# 	plt.legend(handles,labels,title="severity")
# 	if save:
# 		plt.save_fig("../../output/plots/AMS_{}.png".format(column))
# 	plt.show()



def plot_GVH(AMS_file,GVH_file,save):
	df_ams = pd.read_excel(AMS_file)
	df_gvh = pd.read_excel(GVH_file)
	df_gvh = df_gvh.drop([0])
	merged_df = pd.merge(df_ams,df_gvh,how="outer",on="Id CryoStem")
	# merged_df.to_csv("../../output/plots/AMS_table_moelle.tsv",sep="\t")
	merged_df.set_index("Id CryoStem",inplace=True)
	plot_chronic(merged_df,save)

	# plot_gvh_grades(merged_df,save)

def plot_allelic_ratio_vs_DP(df_indiv_file,save):
	df_indiv = pd.read_csv(df_indiv_file,sep="\t")
	x,y = df_indiv["DP"].astype(int).values.reshape(-1,1),df_indiv["ratio"].astype(float).values.reshape(-1,1)
	plt.scatter(x,y)
	plt.xlabel("Depth")
	plt.ylabel("Allelic ratio")
	plt.title("")
	if save:
		plt.savefig("../../output/plots/{}.png".format(df_indiv_file))
	plt.show()
	
def hetero_homo_values(df_indiv):
	heterozygosity = df_indiv1["TYPE"].value_counts()["heterozygous"]
	homozygosity = df_indiv1["TYPE"].value_counts()["homozygous"]
	heterozygosity_ratio = heterozygosity/(homozygosity+heterozygosity)
	return(homozygosity,heterozygosity)
def plot_smth(df_indiv_file_1,df_indiv_file_2):
	df_indiv1,df_indiv2 = pd.read_csv(df_indiv_file_1,sep="\t"),pd.read_csv(df_indiv_file_2,sep="\t")
	len_ind1,len_ind2 = len(df_indiv1),len(df_indiv2)
	hom1,het1 = hetero_homo_values(df_indiv1)
	hom2,het2 = hetero_homo_values(df_indiv2)
	print(len_ind1,len_ind2)
	print(hom1/(hom1+het1),hom2/(hom2+het2))

def plot_DP(df_donor_file,df_recipient_file,merged_df_file,pair,ams_score,n_bins,pdf):
	donor = df_donor_file.split("/")[-1].split(".")[0]
	recipient = df_recipient_file.split("/")[-1].split(".")[0]
	print(pair)
	if "P" not in pair and "R" not in pair:
		pair = "Pair"
	df_donor,df_recipient = pd.read_csv(df_donor_file,sep="\t"),pd.read_csv(df_recipient_file,sep="\t")
	merged_df = pd.read_csv(merged_df_file,sep="\t")
	fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(4,2,sharey="row",figsize=(6,8))
	w1 = np.array(df_donor["DP"].tolist())
	w2 = np.array(df_recipient["DP"].tolist())
	w3 = np.array(merged_df["DP_x"].tolist())
	w4 = np.array(merged_df["DP_y"].tolist())
	w5 = np.array(df_donor["ratio"].tolist())
	w6 = np.array(df_recipient["ratio"].tolist())
	w7 = np.array(merged_df["ratio_x"].tolist())
	w8 = np.array(merged_df["ratio_y"].tolist())
	merged_df[["DP_x","DP_y","ratio_x","ratio_y"]] = merged_df[["DP_x","DP_y","ratio_x","ratio_y"]].astype(float)
	df_donor[["DP","ratio"]] = df_donor[["DP","ratio"]].astype(float)
	df_recipient[["DP","ratio"]] = df_recipient[["DP","ratio"]].astype(float)
	ax1.hist(df_donor["DP"].values.reshape(-1,1),weights = np.zeros_like(w1)+1./w1.size,bins=n_bins)
	ax2.hist(df_recipient["DP"].values.reshape(-1,1),weights = np.zeros_like(w2)+1./w2.size,bins=n_bins)
	ax3.hist(merged_df["DP_x"].values.reshape(-1,1),weights = np.zeros_like(w3)+1./w3.size,bins=n_bins)
	ax4.hist(merged_df["DP_y"].values.reshape(-1,1),weights = np.zeros_like(w4)+1./w4.size,bins=n_bins)
	ax5.hist(df_donor["ratio"].values.reshape(-1,1),weights = np.zeros_like(w5)+1./w5.size,bins=n_bins)
	ax6.hist(df_recipient["ratio"].values.reshape(-1,1),weights = np.zeros_like(w6)+1./w6.size,bins=n_bins)
	ax7.hist(merged_df["ratio_x"].values.reshape(-1,1),weights = np.zeros_like(w7)+1./w7.size,bins=n_bins)
	ax8.hist(merged_df["ratio_y"].values.reshape(-1,1),weights = np.zeros_like(w8)+1./w8.size,bins=n_bins)
	# merged_df["DP_x"].plot(kind = "hist",weights = np.zeros_like(w1)+1./w1.size,bins=n_bins,ax=ax1)
	# df_donor["DP"].plot(kind = "hist",weights = np.zeros_like(w2)+1./w2.size,bins=n_bins,ax=ax2)
	# merged_df["DP_y"].plot(kind = "hist",weights = np.zeros_like(w3)+1./w3.size,bins=n_bins,ax=ax3)
	# df_recipient["DP"].plot(kind = "hist",weights = np.zeros_like(w4)+1./w4.size,bins=n_bins,ax=ax4)
	ax7.set_xlabel("{}".format("".join(donor.split("-filtered"))),fontsize=8)
	ax8.set_xlabel("{}".format("".join(recipient.split("-filtered"))),fontsize=8)
	ax1.set_ylabel("All positions DP")
	ax3.set_ylabel("Mismatches DP")
	ax5.set_ylabel("All positions AR")
	ax7.set_ylabel("Mismatches AR")
	fig.suptitle("{} depth and allelic ratio distributions, AMS = {}".format(pair,ams_score))
	# fig.text(0.3, 0.04, 'common xlabel', ha='center', va='center')
	# fig.text(0.75, 0.04, 'common xlabel', ha='center', va='center')
	# fig.text(0.06, 0.5, 'common ylabel', ha='center', va='center', rotation='vertical')
	# x,y = merged_df["DP_x"].astype(int).values.reshape(-1,1),merged_df["DP_y"].astype(int).values.reshape(-1,1)
	# # x,y = merged_df["ratio_x"].astype(int).values.reshape(-1,1),merged_df["ratio_y"].astype(int).values.reshape(-1,1)
	# # plt.scatter(x,y)
	# x = np.array(merged_df["DP"].tolist())
	# plt.hisintt(x,weights=np.zeros_like(x)+1./x.size)
	pdf.savefig()
	plt.close()



# def plot_ratio_():
# 	x = np.array(df["ratio"].tolist())
# 	plt.hist(x,weights=np.zeros_like(x)+1./x.size)
# 	plt.xlabel("Ratio")
# 	plt.ylabel("Frequency")
# 	plt.title("Allelic Ratio Distribution")
# 	plt.show()


def get_pdf_DP(cohort_path,ams_path,min_DP,max_DP,n_bins):
	df_ams = pd.read_csv(ams_path,sep ="\t")
	directories = [d.path for d in os.scandir(cohort_path) if d.is_dir() and "AMS" not in d.path]
	directories.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
	with PdfPages("../../output/plots/pairs_marrow_dp_ar_{}_{}_bins_{}_plots.pdf".format(min_DP,max_DP,n_bins)) as pdf:
		for directory in directories:
			merged_file = os.path.join(directory,"bed_merged_df_{}_{}.tsv".format(min_DP,max_DP))
			pair = merged_file.split("/")[-2]
			ams_score = int(df_ams["AMS"].loc[df_ams["PAIR"]==pair])
			# d_file,r_file = [os.path.join(directory,f) for f in os.listdir(directory) if "bed" in f and f.endswith("tsv") and "{}_{}".format(min_DP,max_DP) in f and "merged" not in f]
			# print([os.path.join(directory,f) for f in os.listdir(directory) if f.endswith("tsv") and "{}_{}".format(min_DP,max_DP) in f and "merged" not in f and "excluded" not in f])
			d_file,r_file = [os.path.join(directory,f) for f in os.listdir(directory) if f.endswith("tsv") and "{}_{}".format(min_DP,max_DP) in f and "merged" not in f and "excluded" not in f]

			if "donor" in r_file or "D0" in r_file:
				d_file,r_file = r_file,d_file
			plot_DP(d_file,r_file,merged_file,pair,ams_score,n_bins,pdf)

def plot_AMS_shuffle_per_recipient(shuffle_recipient_path,save):
	df = pd.read_csv(shuffle_recipient_path,sep="\t")
	# df["color"] = np.where(df["donor"].str.split("-").str[0]==df["recipient"].str.split("-").str[0],"r","b")
	df["color"] = np.where(df["donor"].str.extract("(P\d+)")==df["recipient"].str.extract("(P\d+)"),"r","b")
	df["donor_type"] = np.where(df["donor"].str.split("-").str[0]==df["recipient"].str.split("-").str[0],"pair donor","test donor")
	df["donor"] = df["donor"].str.split("-donor").str[0]
	# df["donor"] = df["donor"].str.split("-D0").str[0]
	df["recipient"] = df["recipient"].str.split("-recipient").str[0]
	# df["recipient"] = df["recipient"].str.split("-R0").str[0]
	recipient = df["recipient"][0]
	df.plot(kind="bar",x="donor",y="AMS",color=df["color"])
	pair_donor = mpatches.Patch(color='red', label='pair donor')
	test_donor = mpatches.Patch(color='blue', label='test donor')
	plt.xticks(fontsize=8)
	plt.title("AMS of Donors vs recipient {}".format(df["recipient"][0]))
	plt.legend(handles=[pair_donor,test_donor],prop={"size":6})
	if save:
		plt.savefig("../../output/plots/shuffle_kidney/shuffle_{}.png".format(recipient))
	plt.show()

def plot_AMS_shuffle_2_recipients(shuffle_recipient_path1,shuffle_recipient_path2,save):
	recipient1 = re.split("sorted_shuffle_|-filtered|\.",shuffle_recipient_path1.split("/")[-1])[1]
	recipient2 = re.split("sorted_shuffle_|-filtered|\.",shuffle_recipient_path2.split("/")[-1])[1]
	df_recipients = pd.merge(pd.read_csv(shuffle_recipient_path1,sep="\t"),pd.read_csv(shuffle_recipient_path2,sep="\t"),how="inner",on=["donor"])
	df_recipients["donor_pair"] = df_recipients["donor"].str.extract("(P\d+|R\d+)")
	df_recipients["donor_pair"] = df_recipients["donor_pair"].str.replace("R","P")
	print(df_recipients)
	x,y = df_recipients["AMS_x"].astype(int).values.reshape(-1,1),df_recipients["AMS_y"].astype(int).values.reshape(-1,1)
	plt.scatter(x,y,c="c")
	plt.xlabel("AMS Donors vs {}".format(recipient1))
	plt.ylabel("AMS Donors vs {}".format(recipient2))
	plt.title("Correlation of AMS obtained from all donors vs 2 recipients")
	for i in range(len(list(df_recipients["donor_pair"]))):
		if (df_recipients["donor_pair"][i]=="P3215") or (df_recipients["donor_pair"][i]=="P603") or (df_recipients["donor_pair"][i]=="P6241") or (df_recipients["donor_pair"][i]=="P874") or (df_recipients["donor_pair"][i]=="P3820") or (df_recipients["donor_pair"][i]=="P6938"):
			t = plt.text(x[i]*(1 + 0.01), y[i] , df_recipients["donor_pair"][i], fontsize=6)
		
	if save:
		plt.savefig("../../output/plots/marrow/shuffle_marrow/{}_{}_shuffle_scatter.png".format(recipient1,recipient2))
	plt.show()

def plot_SNPs(SNPs_df_path,AMS_path,save):
	SNPs_df = pd.read_csv(SNPs_df_path,sep="\t")
	AMS_df = pd.read_csv(AMS_path,sep="\t")
	# step merge df
	merged_df = pd.merge(SNPs_df,AMS_df,on="pair")
	x,y=merged_df["nb_SNPs_donor"].astype(int).values.reshape(-1,1),merged_df["AMS"].astype(int).values.reshape(-1,1)
	plt.scatter(x,y)
	plt.xlabel("SNPs in donors")
	plt.ylabel("AMS of the pair")
	plt.title("Impact of number of SNPs on AMS")
	for i in range(len(list(SNPs_df["pair"]))):
		if (merged_df["pair"][i]=="P22") or (merged_df["pair"][i]=="P24") or (merged_df["pair"][i]=="P6") or (merged_df["pair"][i]=="P9"):
			t = plt.text(x[i]*(1 + 0.006), y[i] , merged_df["pair"][i], fontsize=8)
	plt.show()
	if save:
		plt.savefig("../../output/plots/marrow/SNPs_vs_AMS.png")

def heatmap_AMS_shuffle(shuffle_df_path,shuffle_donor_df_path,save,AMS_df_path=""):
	if AMS_df_path != "":
		shuffle_df_path = table_ops.sort_shuffle_dataframe(shuffle_df_path,AMS_df_path)
	df = table_ops.reshape_dataframe(shuffle_df_path,shuffle_donor_df_path)
	plt.figure(figsize=(20,20))
	if "GQPDOMB" in shuffle_df_path:
		cols = [re.compile("-(P\d+-donor)").search(col).group(1) for col in list(df.columns)]
		indices = [re.compile("-(P\d+-recipient)").search(col).group(1) for col in list(df.index)]
		df.columns = cols
		df.index = indices
	sns.heatmap(df,xticklabels=True,yticklabels=True,cmap="YlOrBr",annot=True,annot_kws={"size":35/np.sqrt(len(df))},fmt="g")
	plt.title("AMS of all possible donor/recipient pairs")
	if save:
		plt.savefig("../../output/plots/kidney/shuffle_kidney/heatmap_AMS_shuffle_all_pairs.png")
	plt.show()

def prepare_DP_SNPs_count(cohort_directory,filtered=False):
	if filtered:
		donors = [file for file in glob.glob(os.path.join(cohort_directory,"**/*.tsv")) if "side" in file and "filtered" in file and "bed" in file and "donor" in file]
		recipients = [file for file in glob.glob(os.path.join(cohort_directory,"**/*.tsv")) if "side" in file and "filtered" in file and "bed" in file and "recipient" in file]
	else:
		donors = [file for file in glob.glob(os.path.join(cohort_directory,"**/*.tsv")) if "side" in file and "D0" in file and "bed" not in file]
		recipients = [file for file in glob.glob(os.path.join(cohort_directory,"**/*.tsv")) if "side" in file and "R0" in file and "bed" not in file]
	# SORT
	donors.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
	recipients.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
	ls_d = []
	ls_r = []
	for d,r in zip(donors,recipients):
		df_d = pd.read_csv(d,sep="\t")
		df_r = pd.read_csv(r,sep="\t")
		sub_d = pd.DataFrame(df_d["DP"].value_counts().sort_index())
		sub_d = sub_d.rename({"DP":"count"},axis=1)
		sub_d["DP"] = sub_d.index
		sub_d = sub_d.reset_index(drop=True)
		sub_d["r_ecdf"] = sub_d["count"].sum() - sub_d["count"].cumsum()
		sub_d["pair"] = "_".join(Path(d).stem.split("_")[:1])
		sub_d = sub_d[["pair","DP","count","r_ecdf"]]
		sub_r = pd.DataFrame(df_r["DP"].value_counts().sort_index())
		sub_r = sub_r.rename({"DP":"count"},axis=1)
		sub_r["DP"] = sub_r.index
		sub_r = sub_r.reset_index(drop=True)
		sub_r["r_ecdf"] = sub_r["count"].sum() - sub_r["count"].cumsum()
		sub_r["pair"] = "_".join(Path(r).stem.split("_")[:1])
		sub_r = sub_r[["pair","DP","count","r_ecdf"]]
		ls_d.append(sub_d)
		ls_r.append(sub_r)
	pd.concat(ls_d).to_csv(os.path.join("{}".format(cohort_directory),"SNPs_donors_reverse_ecdf.csv"),index=False)
	pd.concat(ls_r).to_csv(os.path.join("{}".format(cohort_directory),"SNPs_recipients_reverse_ecdf.csv"),index=False)

		# fig = plt.figure(figsize=(20,10))
		
		# ax1 = fig.add_subplot(221)
		# xd,yd1,yd2 = sub_d["DP"].values.reshape(-1,1),sub_d["count"].values.reshape(-1,1),sub_d["r_ecdf"].values.reshape(-1,1)
		# ax1.plot(xd,yd1)
		# ax1.title.set_text("{} DP count".format(d.split("/")[-1].split("_")[0]))
		# plt.xlabel("DP")
		# plt.ylabel("count")
		


		# ax2 = fig.add_subplot(222)
		# ax2.plot(xd,yd2)
		# ax2.title.set_text("{} DP cumulative function".format(d.split("/")[-1].split("_")[0]))
		# plt.xlabel("DP")
		# plt.ylabel("cumsum DP")

		# ax3 = fig.add_subplot(223)
		# xr,yr1,yr2 = sub_r["DP"].values.reshape(-1,1),sub_r["count"].values.reshape(-1,1),sub_r["r_ecdf"].values.reshape(-1,1)
		# ax3.plot(xr,yr1)
		# ax3.title.set_text("{} DP count".format(r.split("/")[-1].split("_")[0]))
		# plt.xlabel("DP")
		# plt.ylabel("count")
		

		# ax4 = fig.add_subplot(224)
		# ax4.plot(xr,yr2)
		# ax4.title.set_text("{} DP cumulative function".format(r.split("/")[-1].split("_")[0]))
		# plt.xlabel("DP")
		# plt.ylabel("cumsum DP")

		# plt.show()
		# break

def DP_2_count(cohort_directory):
	donors = [file for file in glob.glob(os.path.join(cohort_directory,"**/*.tsv")) if "side" in file and "D0" in file and ("R927" in file or "R2078" in file) and "bed" in file]
	recipients = [file for file in glob.glob(os.path.join(cohort_directory,"**/*.tsv")) if "side" in file and "R0" in file and ("R3298" in file or "R3820" in file) and "bed" in file]
	recipients = [recipients[1],recipients[0]]
	df_d_0 = pd.read_csv(donors[0],sep="\t")
	print(df_d_0.columns)
	df_r_0 = pd.read_csv(recipients[0],sep="\t")
	sub_d_0 = pd.DataFrame(df_d_0["DP"].value_counts().sort_index())
	sub_d_0 = sub_d_0.rename({"DP":"count"},axis=1)
	sub_d_0["DP"] = sub_d_0.index
	sub_d_0 = sub_d_0.reset_index(drop=True)
	sub_d_0["r_ecdf"] = sub_d_0["count"].sum() - sub_d_0["count"].cumsum()
	sub_r_0 = pd.DataFrame(df_r_0["DP"].value_counts().sort_index())
	sub_r_0 = sub_r_0.rename({"DP":"count"},axis=1)
	sub_r_0["DP"] = sub_r_0.index
	sub_r_0 = sub_r_0.reset_index(drop=True)
	sub_r_0["r_ecdf"] = sub_r_0["count"].sum() - sub_r_0["count"].cumsum()
	df_d_1 = pd.read_csv(donors[1],sep="\t")
	df_r_1 = pd.read_csv(recipients[1],sep="\t")
	sub_d_1 = pd.DataFrame(df_d_1["DP"].value_counts().sort_index())
	sub_d_1 = sub_d_1.rename({"DP":"count"},axis=1)
	sub_d_1["DP"] = sub_d_1.index
	sub_d_1 = sub_d_1.reset_index(drop=True)
	sub_d_1["r_ecdf"] = sub_d_1["count"].sum() - sub_d_1["count"].cumsum()
	sub_r_1 = pd.DataFrame(df_r_1["DP"].value_counts().sort_index())
	sub_r_1 = sub_r_1.rename({"DP":"count"},axis=1)
	sub_r_1["DP"] = sub_r_1.index
	sub_r_1 = sub_r_1.reset_index(drop=True)
	sub_r_1["r_ecdf"] = sub_r_1["count"].sum() - sub_r_1["count"].cumsum()
	fig = plt.figure(figsize=(20,10))
	ax1 = fig.add_subplot(221)
	xd_0,yd_0 = sub_d_0["DP"].values.reshape(-1,1),sub_d_0["count"].values.reshape(-1,1)
	xd_1,yd_1 = sub_d_1["DP"].values.reshape(-1,1),sub_d_1["count"].values.reshape(-1,1)
	ax1.plot(xd_0,yd_0,label="{}, approx SNPs: {}, DP_mean : {:.2f}".format(donors[0].split("/")[-1].split("_")[0],len(df_d_0.index),df_d_0["DP"].mean()))
	ax1.plot(xd_1,yd_1,label="{}, approx SNPs: {}, DP_mean : {:.2f}".format(donors[1].split("/")[-1].split("_")[0],len(df_d_1.index),df_d_1["DP"].mean()))
	ax1.title.set_text("{} and {} DP count".format(donors[0].split("/")[-1].split("_")[0],donors[1].split("/")[-1].split("_")[0]))
	plt.xlabel("DP")
	plt.ylabel("count")
	plt.legend()
	ax2 = fig.add_subplot(222)
	xr_0,yr_0 = sub_r_0["DP"].values.reshape(-1,1),sub_r_0["count"].values.reshape(-1,1)
	xr_1,yr_1 = sub_r_1["DP"].values.reshape(-1,1),sub_r_1["count"].values.reshape(-1,1)
	ax2.plot(xr_0,yr_0,label="{}, approx SNPs: {}, DP_mean : {:.2f}".format(recipients[0].split("/")[-1].split("_")[0],len(df_r_0.index),df_r_0["DP"].mean()))
	ax2.plot(xr_1,yr_1,label="{}, approx SNPs: {}, DP_mean : {:.2f}".format(recipients[1].split("/")[-1].split("_")[0],len(df_r_1.index),df_r_1["DP"].mean()))
	ax2.title.set_text("{} and {} DP count".format(recipients[0].split("/")[-1].split("_")[0],recipients[1].split("/")[-1].split("_")[0]))
	plt.xlabel("DP")
	plt.ylabel("count")
	plt.legend()
	ax3 = fig.add_subplot(223)
	yd_01,yd_02 = sub_d_0["r_ecdf"].values.reshape(-1,1),sub_d_1["r_ecdf"].values.reshape(-1,1)
	ax3.plot(xd_0,yd_01,label="{}, approx SNPs: {}, DP_mean : {:.2f}".format(donors[0].split("/")[-1].split("_")[0],len(df_d_0.index),df_d_0["DP"].mean()))
	ax3.plot(xd_1,yd_02,label="{}, approx SNPs: {}, DP_mean : {:.2f}".format(donors[1].split("/")[-1].split("_")[0],len(df_d_1.index),df_d_1["DP"].mean()))
	ax3.title.set_text("{} and {} DP inverse cumulative function".format(donors[0].split("/")[-1].split("_")[0],donors[1].split("/")[-1].split("_")[0]))
	plt.xlabel("DP")
	plt.ylabel("inverse cumsum DP")
	plt.legend()
	ax4 = fig.add_subplot(224)
	yr_01,yr_02 = sub_r_0["r_ecdf"].values.reshape(-1,1),sub_r_1["r_ecdf"].values.reshape(-1,1)
	ax4.plot(xr_0,yr_01,label="{}, approx SNPs: {}, DP_mean : {:.2f}".format(recipients[0].split("/")[-1].split("_")[0],len(df_r_0.index),df_r_0["DP"].mean()))
	ax4.plot(xr_1,yr_02,label="{}, approx SNPs: {}, DP_mean : {:.2f}".format(recipients[1].split("/")[-1].split("_")[0],len(df_r_1.index),df_r_1["DP"].mean()))
	ax4.title.set_text("{} and {} DP inverse cumulative function".format(recipients[0].split("/")[-1].split("_")[0],recipients[1].split("/")[-1].split("_")[0]))
	plt.xlabel("DP")
	plt.ylabel("inverse cumsum DP")
	plt.legend()
	plt.show()

def windows_AMS(merged,window_size,ls_chr,pair,pdf):
	for chrom in ls_chr:
		subset = merged[merged["CHROM"]==chrom]
		pos = subset["POS"].tolist()
		intervals = np.arange(pos[0],pos[-1]+window_size,window_size)
		nintervals = [str(i) for i in intervals]
		counts = subset.groupby(pd.cut(subset["POS"],list(intervals),include_lowest=True))["CHROM"].count().reset_index()
		counts = counts.rename({"POS":"intervals","CHROM":"counts"},axis=1)
		print(counts)
		fig,ax=plt.subplots(figsize=(12,8))
		ax.bar(x=nintervals[:-1],height=counts["counts"].values,tick_label=nintervals[:-1],align="edge")
		# ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
		# ax.xaxis.set_major_locator(plt.MaxNLocator(5))
		ax.set_xticklabels(nintervals[:-1])
		plt.xticks(rotation=90,fontsize=6)
		plt.title("Barplot of positions contributing to the AMS accross the chromosome {} for pair {}".format(chrom,pair))
		plt.xlabel("positions")
		plt.ylabel("counts")
		pdf.savefig()
		plt.close()

def get_pdf_windows(path,window_size,min_DP,max_DP,min_AD,heterozygosity_thr):
	directories = [d.path for d in os.scandir(path) if d.is_dir() and "AMS" not in d.path]
	directories.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
	with PdfPages("../../output/plots/marrow/AMS_contributions_window_{}.pdf".format(window_size)) as pdf:
		for directory in directories:
			merged_file = os.path.join(directory,"bed_merged_df_bed_{}_{}_{}_{}_GROUPBY.tsv".format(min_DP,max_DP,min_AD,heterozygosity_thr))
			pair = merged_file.split("/")[-2]
			merged = pd.read_csv(merged_file,sep="\t")
			ls_chr = list(merged["CHROM"].drop_duplicates())
			windows_AMS(merged,window_size,ls_chr,pair,pdf)

def prepare_counts(counts,pairs_dict):
	counts["pair_num"] = counts["pair"].str.split("P|R").str.join("")
	counts["sorter"] = counts["pair_num"].map(pairs_dict)
	counts["CHR"] = np.select([counts["CHROM"].str.isdigit(),counts["CHROM"]=="X",counts["CHROM"]=="Y"],[counts["CHROM"],"23","24"])
	counts["CHR"] = counts["CHR"].astype(int)
	return(counts)

def AMS_contributions_dataframe(directory,ams_path,window_size,min_DP,max_DP,min_AD,heterozygosity_thr,thresh):
	directories = [d.path for d in os.scandir(directory) if d.is_dir() and "AMS" not in d.path]
	# sort directories
	directories.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
	dfs = []
	for folder in directories:
		merged_file = os.path.join(folder,"bed_merged_df_{}_{}_{}_{}_TRANSCRIPTS_GENES_GQ_20.tsv".format(min_DP,max_DP,min_AD,heterozygosity_thr))
		pair = merged_file.split("/")[-2]
		merged = pd.read_csv(merged_file,sep="\t")
		merged["pair"] = pair
		merged["CHROM"] = merged["CHROM"].astype(str)
		indivs = [col for col in merged.columns if pair in col]
		merged = merged.drop(indivs,axis=1)
		dfs.append(merged)
	full_merged = pd.concat(dfs).reset_index(drop=True)
	df = pd.read_csv(ams_path,sep="\t")
	pairs_list = df["pair"].tolist()
	pairs_dict = {"".join(pair.split("R")):ind for (pair,ind) in zip(pairs_list,list(range(len(pairs_list))))}
	if window_size != 0:
		full_merged["interval"] = full_merged["POS"]//window_size
		counts = full_merged[["pair","CHROM","interval"]].value_counts(sort=False).reset_index()
		counts = prepare_counts(counts,pairs_dict)
		sorted_int_df = counts.sort_values(by=["CHR","interval"]).reset_index(drop=True)
		sorted_int_df["interval"] = sorted_int_df["CHROM"]+"-"+sorted_int_df["interval"].astype(str)
		sorted_intervals = sorted_int_df["interval"].unique().tolist()
		counts = counts.sort_values(by=["sorter","CHR","interval"]).reset_index(drop=True)
		counts["interval"] = counts["CHROM"]+"-"+counts["interval"].astype(str)
		counts = counts.drop(["pair_num","sorter","CHR","CHROM"],axis=1)
		counts = counts.rename({0:"counts"},axis=1)
		counts = counts.fillna(0)
		counts = counts[counts["counts"]>=thresh]
		counts.to_csv(os.path.join(directory,"csh_{}_counts_chr_long_{}_{}_{}_{}_{}_{}.csv".format(date.today().strftime("%Y_%m_%d"),window_size,min_DP,max_DP,min_AD,heterozygosity_thr,thresh)),index=False)
		sorted_pairs = counts["pair"].unique().tolist()
		counts = counts.groupby(["pair","interval"])["counts"].aggregate("first").unstack()
		counts = counts.reindex(sorted_pairs)
		counts = counts.reindex(sorted_intervals,axis=1)
		counts.to_csv(os.path.join(directory,"csh_{}_counts_heatmap_{}_{}_{}_{}_{}_{}.csv").format(date.today().strftime("%Y_%m_%d"),window_size,min_DP,max_DP,min_AD,heterozygosity_thr,thresh))
		fig,ax = plt.subplots(figsize=(20,20))
		sns.heatmap(counts,yticklabels=True,cmap="YlGnBu")
		plt.title("Heatmap of AMS per window of size {} bases with {} local AMS min".format(window_size,thresh))
		plt.xlabel("position in {} bases, format : chr-POS/{}".format(window_size,window_size))
	else:
		counts = full_merged[["pair","CHROM"]].value_counts(sort=False).reset_index()
		counts = prepare_counts(counts,pairs_dict)
		counts = counts.sort_values(by=["sorter","CHR"]).reset_index(drop=True)
		sorted_chr = counts["CHROM"].unique().tolist()
		print(counts)
		counts = counts.drop(["pair_num","sorter","CHR"],axis=1)
		counts = counts.rename({0:"counts"},axis=1)
		counts = counts.fillna(0)
		counts = counts[counts["counts"]>=thresh]
		counts.to_csv(os.path.join(directory,"csh_{}_counts_chr_long_{}_{}_{}_{}_{}_{}.csv".format(date.today().strftime("%Y_%m_%d"),window_size,min_DP,max_DP,min_AD,heterozygosity_thr,thresh)),index=False)
		sorted_pairs = counts["pair"].unique().tolist()
		counts.to_csv(os.path.join(directory,"counts_heatmap_CHROM_long_{}_{}_{}_{}_{}.csv").format(min_DP,max_DP,min_AD,heterozygosity_thr,thresh))
		counts = counts.groupby(["pair","CHROM"])["counts"].aggregate("first").unstack()
		counts = counts.reindex(sorted_pairs)
		counts = counts.reindex(sorted_chr,axis=1)
		# counts.to_csv(os.path.join(directory,"counts_heatmap_CHROM_{}_{}_{}_{}_{}.csv").format(min_DP,max_DP,min_AD,heterozygosity_thr,thresh))
		# print(counts)
		fig,ax = plt.subplots(figsize=(20,20))
		sns.heatmap(counts,yticklabels=True,cmap="YlGnBu")
		plt.title("Heatmap of AMS per chromosome with {} local AMS min".format(thresh))
		plt.xlabel("Chromosome")

	plt.show()
	return(full_merged)

def main(args):
	# AMS_file,GVH_file = args[0],args[1]
	# df_indiv_file_1,df_indiv_file_2 = args[0],args[1]
	# df_donor_file,df_recipient_file,merged_df_file,cohort = args[0],args[1],args[2],args[3]
	# cohort_path,ams_path,min_DP,max_DP,n_bins = args[0],args[1],int(args[2]),int(args[3]),int(args[4])
	# plot_AMS(AMS_file,"nyc",False)
	# plot_eGFR_corr(AMS_file,eGFR_file,"nyc",True)
	# plot_GVH(AMS_file,GVH_file,False)
	# plot_allelic_ratio_vs_DP(df_indiv_file_1,False)
	# plot_smth(df_indiv_file_1,df_indiv_file_2)
	# files = [f for f in glob.glob("../../output/indiv_vcf/joint_genotyping/hard-filtered/**/bed_merged_df_10_500.tsv")]
	# for f in files:
	# plot_DP_position("../../output/indiv_vcf/joint_genotyping/hard-filtered/R3215/R3215-R0-N1_10_500.tsv")
	# plot_DP_position("../../output/indiv_vcf/joint_genotyping/hard-filtered/R3215/bed_merged_df_10_500.tsv")
	# plot_DP_position("../../output/indiv_vcf/cohort_GQPDOMB/P16/bed_merged_df_10_1000.tsv")
	# plot_DP(df_donor_file,df_recipient_file,merged_df_file)
	# get_pdf_DP(cohort_path,ams_path,min_DP,max_DP,n_bins)
	save = False
	# shuffle_files = [file for file in glob.glob("../../output/indiv_vcf/joint_genotyping/hard-filtered/AMS_shuffle/**.tsv") if "sorted" in file]
	# for file in shuffle_files:
	# 	plot_AMS_shuffle_per_recipient("{}".format(file),save)
	# plot_AMS_shuffle_2_recipients(shuffle_files[0],shuffle_files[10],save)
	# heatmap_AMS_shuffle("../../output/indiv_vcf/cohort_GQPDOMB/AMS_shuffle/shuffle_AMS_2021_sorted_on_recipients.tsv","../../output/indiv_vcf/cohort_GQPDOMB/AMS_shuffle/shuffle_AMS_2021_sorted_on_donors.tsv",False)
	# heatmap_AMS_shuffle("../../output/indiv_vcf/joint_genotyping/hard-filtered/AMS_shuffle/shuffle_AMS_2021.tsv","../../output/indiv_vcf/joint_genotyping/hard-filtered/AMS_shuffle/shuffle_AMS_2021_sorted_on_donors.tsv",False)
	# DP_count("../../output/indiv_vcf/joint_genotyping/hard-filtered")
	# merged = pd.read_csv("../../output/indiv_vcf/joint_genotyping/hard-filtered/R2107/bed_merged_df_bed_20_400_5_0.8_GROUPBY.tsv",sep="\t")
	window_size = 10000
	thresh = 0
	# ls_chr = list(merged["CHROM"].drop_duplicates())
	# windows_AMS(merged,window_size,ls_chr,pair)
	d = AMS_contributions_dataframe("../../output/indiv_vcf/joint_genotyping/hard-filtered","../../output/indiv_vcf/joint_genotyping/hard-filtered/AMS_folder/AMS_20_400_5_0.8_rd/AMS_df_sorted.tsv",window_size,20,400,5,0.8,thresh)
if __name__ == '__main__':
	main(sys.argv[1:])


files = glob.glob("**/*.HC.annot.vcf.gz")
print(len(files))
ch = "vcf-merge"
for file in files:
	ch += " "+file
os.system("{} | bgzip -c > out.vcf.gz".format(ch))

