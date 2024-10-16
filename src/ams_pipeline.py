# coding : utf-8
"""
The AMS pipeline estimates the mismatch between a donor and a recipient based on VCF information
command line help : python3 ams_pipeline.py [-h]
"""
import os
import sys
import datetime
from tools import ams_helpers, arguments_handling, table_operations


def main():
    """
    The AMS pipeline estimates the mismatch between a donor and a recipient based on VCF information
    command line help : python3 ams_pipeline.py [-h]
    """
    # get args from cmd line
    args = arguments_handling.arguments()
    process = ams_helpers.handle_overwrite(args)
    formatted_datetime=""
    # append timestamp instead of overwriting
    if process:
        current_datetime = datetime.datetime.now()
        timestamp = current_datetime.timestamp()
        datetime_object = datetime.datetime.fromtimestamp(timestamp)
        formatted_datetime = datetime_object.strftime("_%Y-%m-%d_%H-%M-%S")
        print("Warning: directory already exists and contains an AMS file with the same run parameters.\nTimestamp will be added to output files names instead of overwriting.")
    # create run dependencies
    run_path, run_tables, run_plots, run_ams = ams_helpers.create_run_directory(
        args.run_name
    )
    # filter the donor file
    df_donor, vep_donor, vep_indices_donor = ams_helpers.prepare_indiv_df(
        run_tables, args.donor, args, args.wc_donor, formatted_datetime
    )
    # output : filtered donor file
    df_donor_file = (
        f"{run_tables}/{args.pair}_{args.run_name}_"
        f"D0_{args.donor.split('/')[-1].split('.')[0]}_"
        f"{args.min_dp}_{args.max_dp}_{args.min_ad}_"
        f"gq_{args.min_gq}_{args.homozygosity_thr}_"
        f"bl_{args.base_length}{formatted_datetime}.tsv"
    )
    df_donor.to_csv(df_donor_file, sep="\t", index=False)
    # same for recipient
    df_recipient, vep_recipient, vep_indices_recipient = ams_helpers.prepare_indiv_df(
        run_tables, args.recipient, args, args.wc_recipient, formatted_datetime
    )
    df_recipient_file = (
        f"{run_tables}/{args.pair}_{args.run_name}_"
        f"R0_{args.recipient.split('/')[-1].split('.')[0]}_"
        f"{args.min_dp}_{args.max_dp}_{args.min_ad}_"
        f"gq_{args.min_gq}_{args.homozygosity_thr}_"
        f"bl_{args.base_length}{formatted_datetime}.tsv"
    )
    df_recipient.to_csv(df_recipient_file, sep="\t", index=False)
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
        mismatches_file = (
            f"{run_tables}/{args.pair}_"
            f"{args.run_name}_mismatches_{args.min_dp}_"
            f"{args.max_dp}_{args.min_ad}_gq_{args.min_gq}_"
            f"{args.homozygosity_thr}_bl_{args.base_length}{formatted_datetime}.tsv"
        )
        merged_df.to_csv(mismatches_file, sep="\t", index=False)

        transcripts_donor = table_operations.build_transcripts_table_indiv(
            vep_donor, merged_df, vep_indices_donor, "donor"
        )
        transcripts_recipient = table_operations.build_transcripts_table_indiv(
            vep_recipient, merged_df, vep_indices_recipient, "recipient"
        )
        merged_transcripts = table_operations.build_transcripts_table(
            transcripts_donor, transcripts_recipient, merged_df
        )
        transcripts_file = (
            os.path.join(
                run_tables,
                f"{args.pair}_{args.run_name}_transcripts_pair_codons_"
                f"{args.min_dp}_{args.max_dp}_{args.min_ad}_gq_"
                f"{args.min_gq}_{args.homozygosity_thr}_bl_"
                f"{args.base_length}{formatted_datetime}.tsv",
            )
        )
        merged_transcripts.to_csv(transcripts_file, sep="\t", index=False)
    else:
        raise ValueError("No mismatch found. No output table will be generated.")
    # create small df with pair and AMS
    ams_exp_path = table_operations.save_mismatch(
        run_ams, args, mismatch, formatted_datetime, df_donor_file, df_recipient_file, mismatches_file
    )
    print(mismatch)
    
    # score normalization (multiprocess_ams only)
    if args.norm_score:
        ams_df, ams_exp_path = table_operations.get_ref_ratio(
            args.run_name, args.pair, run_path, ams_exp_path, args.donor, args.recipient,
            args.min_dp, args.max_dp, args.min_ad, args.homozygosity_thr, args.min_gq,
            args.orientation, args.base_length
        )
        table_operations.add_norm(ams_df, ams_exp_path)
    
    return 0


if __name__ == "__main__":
    main()
