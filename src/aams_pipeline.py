#coding:utf-8
"""
The AAMS pipeline estimates the mismatch between a donor and a recipient after running NetMHCpan
command line help : python3 aams_pipeline.py [-h]
"""
import sys
from tools import netmhc_arguments, aams_helpers, netmhc_tables_handler, cleavage


def main():
    """
    The AAMS pipeline estimates the mismatch between a donor and a recipient after running NetMHCpan
    command line help : python3 aams_pipeline.py [-h]
    """
    args = netmhc_arguments.netmhc_arguments()
    str_params, mismatches_path = aams_helpers.get_ams_params(args.run_name, args.output_dir)
    aams_run_tables, netmhc_dir, aams_path, netchop_dir = aams_helpers.create_aams_dependencies(
        args.run_name, args.output_dir
    )
    log_file = aams_helpers.append_log(args)
    if args.cleavage == True:
        try:
            cleavage.validate_cleavage_compatibility(log_file)
        except ValueError as err:
            print(err)
            return 1
    fasta_path, pep_indiv_path, ens_transcripts, peptides_ensembl, refseq_file, pair_print = aams_helpers.build_peptides(
        aams_run_tables, str_params, args, log_file, mismatches_path, cleavage_mode=False
    )
    if args.cleavage == True:
        print(f"{pair_print}Entering NetChop handler: running NetChop may last a few minutes...")
        pep_paths = {}
        for sample in ("donor", "recipient"):
            sample_suffix = f"_{sample}"
            pickle_df = cleavage.pickle_parsing(str_params, args, log_file, sample)
            mismatches_df, transcripts_pair, peptides_ensembl, pair_print = aams_helpers.build_peptides(
                aams_run_tables, str_params, args, log_file, mismatches_path, mismatches_df=pickle_df, cleavage_mode=args.cleavage,
                ens_transcripts=ens_transcripts, peptides_ensembl=peptides_ensembl, refseq_file=refseq_file
            )
            chop_table, chop_table_path = cleavage.netchop_table_prep(
                mismatches_df, transcripts_pair, peptides_ensembl, args, netchop_dir, sample_suffix
            )
            chop_output = cleavage.run_netchop(chop_table, args, netchop_dir, sample_suffix)
            pep_paths[sample] = cleavage.postprocess_netchop(chop_output, chop_table_path, args, netchop_dir, sample_suffix)
        deduced_pep_path = cleavage.deduce_cleaved_peptides(pep_paths["donor"], pep_paths["recipient"], netchop_dir, args, log_file, pair_print)
        fasta_path, pep_indiv_path = cleavage.prepare_cleavage_netmhcpan_inputs(
            aams_run_tables, args, pep_indiv_path, deduced_pep_path
        )
        if cleavage.fasta_is_empty(fasta_path):
            print(f"{pair_print}: No cleaved peptides available for NetMHCpan; skipping affinity prediction.")
            return 0
    if args.dry_run == False:
        netmhc_out = aams_helpers.run_netmhcpan(fasta_path, netmhc_dir, args, pair_print)
        netmhc_table = netmhc_tables_handler.handle_netMHCpan(netmhc_out, args)
        netmhc_df, pep_df = aams_helpers.clean_pep_df(netmhc_table, pep_indiv_path, args)
        if args.cleavage == False:
            mismatch = aams_helpers.merge_netmhc(
                netmhc_df, pep_df, mismatches_path,
                aams_path,aams_run_tables,str_params, args
            )
            print(mismatch)
            return 0

if __name__ == '__main__':
    sys.exit(main())
