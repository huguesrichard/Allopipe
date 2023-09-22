# coding : utf-8
# command line help : python3 comparison.py [-h]
# python3 comparison.py ../output/indiv_vcf/cohort_GQPDOMB/P1/BNAGPIP-paris-SCE240246-P1-donor-filtered.bed.vcf.gz ../output/indiv_vcf/cohort_GQPDOMB/P1/ONDAMJK-paris-SCE260544-P1-recipient-filtered.bed.vcf.gz 10 500 5 0.95 dr
# python3 comparison.py ../output/indiv_vcf/joint_genotyping/hard-filtered/R383/R383-D0-N1.vcf.gz ../output/indiv_vcf/joint_genotyping/hard-filtered/R383/R383-R0-N1.vcf.gz 10 500 5 0.95 rd
# python3 comparison.py ../output/indiv_vcf/cohort_GQPDOMB/P1/BNAGPIP-paris-SCE240246-P1-donor-filtered.bed.vcf.gz ../output/indiv_vcf/cohort_GQPDOMB/P1/ONDAMJK-paris-SCE260544-P1-recipient-filtered.bed.vcf.gz 10 500 5 0.95 dr -wc -wcd ../output/GQPDOMB-VEP.txt -wcr ../output/GQPDOMB-VEP.txt
# python3 comparison.py ../output/indiv_vcf/joint_genotyping/hard-filtered/R383/R383-D0-N1.vcf.gz ../output/indiv_vcf/joint_genotyping/hard-filtered/R383/R383-R0-N1.vcf.gz 10 500 5 0.95 dr -wc -wcd ../output/indiv_vcf/joint_genotyping/hard-filtered/allsamples.jointgenotyping.hard-filtered.most_severe.txt -wcr ../output/indiv_vcf/joint_genotyping/hard-filtered/allsamples.jointgenotyping.hard-filtered.most_severe.txt
import sys
import numpy as np
import pandas as pd
import tools.old_helpers as helpers
import tools.arguments_handling as argh
import tools.table_operations as table_ops


def main():
	args = argh.arguments()
	donor_path = "/".join(args.donor.split("/")[0:-1])
	pair = donor_path.split("/")[-1]
	donor = args.donor.split("/")[-1].split(".")[0]
	if "bed" in args.donor:
		donor+="-bed"
	print(donor)
	# df_donor = helpers.prepare_indiv_df(args.donor,args.min_DP,args.max_DP,args.min_AD,args.homozygosity_thr,args.wc,args.wc_donor)
	# save
	# df_donor.to_csv("{}/{}_{}_{}_side.tsv".format(donor_path,donor,args.min_DP,args.max_DP),sep="\t")
	# df_donor.to_csv("{}/df1.tsv".format(donor_path),sep="\t")


	# df_recipient = helpers.prepare_indiv_df(args.recipient,args.min_DP,args.max_DP,args.min_AD,args.homozygosity_thr,args.wc,args.wc_recipient)
	recipient_path = "/".join(args.recipient.split("/")[0:-1])
	recipient = args.recipient.split("/")[-1].split(".")[0]
	if "bed" in args.recipient:
		recipient+="-bed"
	# save
	# df_recipient.to_csv("{}/{}_{}_{}_side.tsv".format(recipient_path,recipient,args.min_DP,args.max_DP),sep="\t")
	# df_recipient.to_csv("{}/df2.tsv".format(recipient_path),sep="\t")
	
	# Can run these lines if the above ones were run previously
	df_donor = pd.read_csv("{}/{}_{}_{}_side.tsv".format(donor_path,donor,args.min_DP,args.max_DP),sep="\t")
	df_recipient = pd.read_csv("{}/{}_{}_{}_side.tsv".format(recipient_path,recipient,args.min_DP,args.max_DP),sep="\t")

	merged_df,side,opposite = helpers.merge_dfs(df_donor,df_recipient,args.orientation)
	merged_df = helpers.keep_alt(merged_df,side,opposite)
	# print(merged_df)
	merged_df,mismatch = helpers.count_mismatches(merged_df,args.orientation)
	merged_df = helpers.mismatch_type(merged_df)
	subset = merged_df[["CHROM","POS","REF","ALT","aa_indiv_x","aa_indiv_y","aa_REF","aa_ALT","diff","mismatch","mismatch_type"]].sample(50).sort_index()
	# print(merged_df[["aa_indiv_x","aa_indiv_y","diff","mismatch"]][merged_df["mismatch"]==1])
	# print(subset)
	# print(mismatch)
	# save
	merged_df.to_csv("{}/bed_merged_df_{}_{}_side.tsv".format(donor_path,args.min_DP,args.max_DP),sep="\t")
	# save
	table_ops.save_mismatch(args.donor,mismatch,args.min_DP,args.max_DP,args.min_AD,args.homozygosity_thr)
	# table_ops.shuffle_mismatch(args.donor,args.recipient,mismatch,args.min_DP,args.max_DP,args.min_AD,args.homozygosity_thr)
	# table_ops.create_AMS_df(donor_path.split("/")[-2])

if __name__ == '__main__':
	main()
