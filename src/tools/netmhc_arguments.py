#coding:utf-8
"""
This script was developed with the argparse module to make 
the aams pipeline command line easy to use and user friendly
"""
import os
import argparse
import re
from pathlib import Path
from tools import arguments_handling

def check_if_valid_k(arg):
    """
    Returns the peptide size if the value is valid

    Parameters:
            arg (int): positive integer between 8 and 14
    Returns:
            arg (int): positive integer between 8 and 14
    """
    try:
        ls_args = arg.split(",")
        if not all(key.isdigit() for key in ls_args):
            raise ValueError(f"The list of k sizes {arg} contains non integer values")
        if len(ls_args)==1:
            arg = int(arg)  
            if arg < 0:
                raise ValueError(f"{arg} is not a positive integer")
            if not 8 <= arg <= 14 :
                raise ValueError(f"A length of peptide of {arg} is not supported by netMHCpan")
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"{err}")
    return arg

def check_if_existing_path(parser,arg):
    """
    Returns the directory path if it exists

    Parameters:
            parser (CustomParser Object): parser object
            arg (str): path to directory or file
    Returns:
            arg (pathlib.PosixPath): Path object of the path to the directory or file
    """
    arg = Path(f"{arg}")
    if not arg.exists():
        parser.error(f"No such file or directory : {arg}")
    return arg

def check_if_existing_run_name(parser,arg):
    """
    Returns the directory name if it exists

    Parameters:
            parser (CustomParser Object): parser object
            arg (str): name of a run directory
    Returns:
            arg (str): name of a run directory
    """
    if not Path(os.path.join("../output/runs",f"{arg}")).exists():
        parser.error(f"No such file or directory : {arg}")
    return arg

def check_hla_i_format(parser,arg):
    """
    Returns the list of HLA genes if they are valid

    Parameters:
            parser (CustomParser Object): parser object
            arg (str): list of HLA genes
    Returns:
            arg (str): list of HLA genes
    """
    if not arg.isalnum() and not ("," in arg or "*" in arg):
        parser.error(f"{arg} contains non accepted characters")
    try:
        ls_args = arg.split(",")
        if len(ls_args) != 6:
            parser.error(f"The number of provided HLA class I genes is incorrect")
        if not all(re.match(r"HLA\-[A-C]\*?\d{2}\:\d{2}",key) for key in ls_args):
            raise ValueError(f"The list {arg} does not comply to the expected format")
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"{err}")
    if "*" in arg:
        arg = "".join(arg.split("*"))
    return arg

def check_if_valid_float(parser,arg):
    try:
        thresh = float(arg)
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"{err}")
    if not 0 < thresh <= 100:
        parser.error(f"{threshold} not between 0 and 100%")
    return thresh

def netmhc_arguments():
    """
    Returns the parsed arguments of the aams pipeline
    Returns:
                    args (argparse.Namespace): object containing all parameters
    """
    parser = arguments_handling.CustomParser(
        prog ="aams_pipeline.py",
        usage="python %(prog)s [options]",
        description="Compute the AAMS for a pair of individuals"
        )
    parser.add_argument("-M","--merged",
        help="path of the pair merged file",
        action=arguments_handling.UniqueStore,
        required=True,
        type=lambda x: check_if_existing_path(parser,x))
    parser.add_argument("-T","--transcripts",
        help="path of the pair transcripts file",
        action=arguments_handling.UniqueStore,
        required=True,
        type=lambda x: check_if_existing_path(parser,x))
    parser.add_argument("-E","--ensembl_transcripts",
        help="path of the ensembl transcripts file",
        action=arguments_handling.UniqueStore,
        required=True,
        type=lambda x: check_if_existing_path(parser,x))
    parser.add_argument("-P","--peptides",
        help="path of the ensembl peptides file",
        action=arguments_handling.UniqueStore,
        required=True,
        type=lambda x: check_if_existing_path(parser,x))
    parser.add_argument("-R","--refseq",
        help="path of the ensembl refseq transcripts file",
        action=arguments_handling.UniqueStore,
        required=True,
        type=lambda x: check_if_existing_path(parser,x))
    parser.add_argument(
        "-n",
        "--run_name",
        help="name of the ams pipeline ran previously",
        action=arguments_handling.UniqueStore,
        required=True,
        type=lambda x: check_if_existing_run_name(parser,x))
    parser.add_argument(
        "-p",
        "--pair",
        help="name of the pair",
        action=arguments_handling.UniqueStore,
        nargs="?",
        default="input_pair",
        const="input_pair",
        type=lambda x: arguments_handling.check_if_accepted_str(parser,x))
    # netMHCpan arg
    parser.add_argument("-l",
        "--length",
        help="peptide length, default is 9, a list of values separated by commas is accepted",
        nargs="?",
        default=9,
        const=9,
        type=lambda x: check_if_valid_k(x))
    # netMHCpan arg
    parser.add_argument("-a",
        "--hla_typing",
        help="comma separated list of HLA genes",
        required=True,
        type=lambda x: check_hla_i_format(parser,x)
        )
    parser.add_argument("-e","--el_rank",
        help=r"%%EL-rank filtration, all values above the given value are filtered out",
        default=100,
        type=lambda x: check_if_valid_float(parser,x)
        )
    args = parser.parse_args()
    return args
