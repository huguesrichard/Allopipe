#coding:utf-8
"""
The AAMS pipeline estimates the mismatch between a donor and a recipient after running NetMHCpan
command line help : python3 aams_pipeline.py [-h]
"""
import os
from tools import netmhc_arguments, aams_helpers, netmhc_tables_handler


def main():
    """
    The AAMS pipeline estimates the mismatch between a donor and a recipient after running NetMHCpan
    command line help : python3 aams_pipeline.py [-h]
    """
    # get args from cmd line
    args = netmhc_arguments.netmhc_arguments()
    str_params = aams_helpers.get_ams_params(args.merged)
    aams_run_tables,netmhc_dir,aams_path = aams_helpers.create_aams_dependencies(args.run_name)
    fasta_path,pep_indiv_path = aams_helpers.build_peptides(aams_run_tables,str_params,args)
    netmhc_out = aams_helpers.run_netmhcpan_class_1(fasta_path,netmhc_dir,args)
    netmhc_table = netmhc_tables_handler.handle_netMHCpan(netmhc_out,args)
    netmhc_df, pep_df = aams_helpers.clean_pep_df(netmhc_table, pep_indiv_path)
    mismatch = aams_helpers.merge_netmhc(netmhc_df, pep_df, args.merged, args.el_rank,args.pair,aams_path,aams_run_tables,str_params)
    print(mismatch)
    return 0

if __name__ == '__main__':
    main()
