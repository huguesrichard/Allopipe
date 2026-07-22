#coding:utf-8
"""
The AMS pipeline estimates the mismatch between a donor and a recipient based on VCF information
command line help : python3 ams_pipeline.py [-h]
"""
import os
from tools import ams_helpers, arguments_handling, table_operations, plot_hist, plot_pie


def main():
    """
    The AMS pipeline estimates the mismatch between a donor and a recipient based on VCF information
    command line help : python3 ams_pipeline.py [-h]
    """
    # get args from cmd line
    args = arguments_handling.arguments()
    # create run dependencies
    run_path, run_tables, run_plots, run_ams, run_logs = ams_helpers.create_run_directory(
        args.run_name, args.output_dir
    )
    # basic logging (to pass paramaters from AMS to AAMS)
    ams_helpers.write_log(run_logs, args)
    
    # avoids long time with no print
    pair_print = f"[{args.pair}] " if args.pair else ""
    print(f"{pair_print}Preparing donor and recipient tables...")
    
    # filter the donor file
    df_donor, vep_donor, vep_indices_donor = ams_helpers.prepare_indiv_df(
        run_tables, args.donor, args, args.wc_donor
    )
    # output : filtered donor file
    df_donor_file = (
        f"{run_tables}/"
        f"{args.pair + '_' if args.pair else ''}" # avoid file starts with "_"
        f"{args.run_name}_"
        f"D0_{args.donor.split('/')[-1].split('.')[0]}_"
        f"{args.min_dp}_{args.max_dp}_{args.min_ad}_"
        f"gq_{args.min_gq}_{args.homozygosity_thr}_"
        f"bl_{args.base_length}.tsv"
    )
    df_donor.to_csv(df_donor_file, sep="\t", index=False)
    # same for recipient
    df_recipient, vep_recipient, vep_indices_recipient = ams_helpers.prepare_indiv_df(
        run_tables, args.recipient, args, args.wc_recipient
    )
    df_recipient_file = (
        f"{run_tables}/"
        f"{args.pair + '_' if args.pair else ''}" # avoid file starts with "_"
        f"{args.run_name}_"
        f"R0_{args.recipient.split('/')[-1].split('.')[0]}_"
        f"{args.min_dp}_{args.max_dp}_{args.min_ad}_"
        f"gq_{args.min_gq}_{args.homozygosity_thr}_"
        f"bl_{args.base_length}.tsv"
    )
    df_recipient.to_csv(df_recipient_file, sep="\t", index=False)
    # merge dataframes based on the orientation of the comparison
    merged_df, side, opposite = ams_helpers.merge_dfs(
        df_donor, df_recipient, args.orientation, args.imputation
    )
    # remove REF/REF vs NaN positions
    merged_df = ams_helpers.keep_alt(merged_df, side, opposite)
    # count mismatches
    merged_df, mismatch = ams_helpers.count_mismatches(merged_df, args.orientation)
    if mismatch != 0:
        # estimate if a mismatch is heterozygous or homozygous
        merged_df = ams_helpers.mismatch_type(merged_df)
        # output : merged df with mismatch counts
        mismatches_file = (
            f"{run_tables}/"
            f"{args.pair + '_' if args.pair else ''}" # avoid file starts with "_"
            f"{args.run_name}_mismatches_{args.min_dp}_"
            f"{args.max_dp}_{args.min_ad}_gq_{args.min_gq}_"
            f"{args.homozygosity_thr}_bl_{args.base_length}.tsv"
        )
        merged_df.to_csv(mismatches_file, sep="\t", index=False)

        transcripts_donor = table_operations.build_transcripts_table_indiv(
            vep_donor, merged_df, vep_indices_donor, "donor"
        )
        transcripts_recipient = table_operations.build_transcripts_table_indiv(
            vep_recipient, merged_df, vep_indices_recipient, "recipient"
        )
        merged_transcripts = table_operations.build_transcripts_table(
            transcripts_donor, transcripts_recipient
        )
        transcripts_file = (
            os.path.join(
                run_tables,
                f"{args.pair + '_' if args.pair else ''}" # avoid file starts with "_"
                f"{args.run_name}_transcripts_pair_codons_"
                f"{args.min_dp}_{args.max_dp}_{args.min_ad}_gq_"
                f"{args.min_gq}_{args.homozygosity_thr}_bl_"
                f"{args.base_length}.tsv",
            )
        )
        merged_transcripts.to_csv(transcripts_file, sep="\t", index=False)
    else:
        raise ValueError("No mismatch found. No output table will be generated.")
    # create small df with pair and AMS
    ams_exp_path = table_operations.save_mismatch(
        run_ams, args, mismatch, df_donor_file, df_recipient_file, mismatches_file
    )
    print(f"Mismatches: {mismatch} — pair "
          f"{args.donor.split('/')[-1].split('.')[0]} "
          f"{args.recipient.split('/')[-1].split('.')[0]}"
          )

    # plots
    plot_hist.hist(ams_exp_path, run_plots)
    plot_pie.pie(run_tables, run_plots, args.pair, args.run_name)
    
    return 0


if __name__ == "__main__":
    main()
