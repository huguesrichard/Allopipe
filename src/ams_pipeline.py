# coding : utf-8
"""
The AMS pipeline estimates the mismatch between a donor and a recipient based on VCF information
command line help : python3 ams_pipeline.py [-h]
"""
import os
import sys
from tools import ams_helpers, arguments_handling, table_operations


def main():
    """
    The AMS pipeline estimates the mismatch between a donor and a recipient based on VCF information
    command line help : python3 ams_pipeline.py [-h]
    """
    # get args from cmd line
    args = arguments_handling.arguments()
    process = ams_helpers.handle_overwrite(args)
    if not process:
        sys.exit()
    # create run dependencies
    run_path, run_tables, run_plots, run_ams = ams_helpers.create_run_directory(
        args.run_name
    )
    # filter the donor file
    df_donor, vep_donor, vep_indices_donor = ams_helpers.prepare_indiv_df(
        run_tables, args.donor, args, args.wc_donor
    )
    # output : filtered donor file
    df_donor.to_csv(
        f"{run_tables}/{args.run_name}_"
        f"{args.donor.split('/')[-1].split('.')[0]}_"
        f"{args.min_dp}_{args.max_dp}_{args.min_ad}_"
        f"gq_{args.min_gq}_{args.homozygosity_thr}_"
        f"bl_{args.base_length}.tsv",
        sep="\t",
        index=False,
    )
    # same for recipient
    df_recipient, vep_recipient, vep_indices_recipient = ams_helpers.prepare_indiv_df(
        run_tables, args.recipient, args, args.wc_recipient
    )
    df_recipient.to_csv(
        f"{run_tables}/{args.run_name}_"
        f"{args.recipient.split('/')[-1].split('.')[0]}_"
        f"{args.min_dp}_{args.max_dp}_{args.min_ad}_"
        f"gq_{args.min_gq}_{args.homozygosity_thr}_"
        f"bl_{args.base_length}.tsv",
        sep="\t",
        index=False,
    )

    # merge dataframes based on the orientation of the comparison

    merged_df, side, opposite = ams_helpers.merge_dfs(
        df_donor, df_recipient, args.orientation
    )
    # remove REF/REF vs NaN positions
    merged_df = ams_helpers.keep_alt(merged_df, side, opposite)
    # count mismatches
    merged_df, mismatch = ams_helpers.count_mismatches(merged_df, args.orientation)
    if mismatch != 0:
        # estimate if a mismatch is heterozygous or homozygous
        merged_df = ams_helpers.mismatch_type(merged_df)
        # output : merged df with mismatch counts
        merged_df.to_csv(
            f"{run_tables}/{args.pair}_"
            f"{args.run_name}_mismatches_{args.min_dp}_"
            f"{args.max_dp}_{args.min_ad}_gq_{args.min_gq}_"
            f"{args.homozygosity_thr}_bl_{args.base_length}.tsv",
            sep="\t",
            index=False,
        )

        transcripts_donor = table_operations.build_transcripts_table_indiv(
            vep_donor, merged_df, vep_indices_donor
        )
        transcripts_recipient = table_operations.build_transcripts_table_indiv(
            vep_recipient, merged_df, vep_indices_recipient
        )
        merged_transcripts = table_operations.build_transcripts_table(
            transcripts_donor, transcripts_recipient, merged_df
        )
        merged_transcripts.to_csv(
            os.path.join(
                run_tables,
                f"{args.pair}_{args.run_name}_transcripts_pair_codons_"
                f"{args.min_dp}_{args.max_dp}_{args.min_ad}_gq_"
                f"{args.min_gq}_{args.homozygosity_thr}_bl_"
                f"{args.base_length}.tsv",
            ),
            sep="\t",
            index=False,
        )
    # create small df with pair and AMS
    table_operations.save_mismatch(run_ams, args, mismatch)
    print(mismatch)
    return 0


if __name__ == "__main__":
    main()
