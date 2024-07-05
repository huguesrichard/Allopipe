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

def check_if_valid_k(initial_args, arg):
    """
    Returns the peptide size if the value is valid

    Parameters:
            initial_args (ArgumentParser object): parser object
            arg (int): positive integer between 8 and 14 (class 1) or not less than 9 (class 2)
    Returns:
            arg (int): positive integer between 8 and 14 (class 1) or not less than 9 (class 2)
    """
    try:
        ls_args = arg.split(",")
        if not all(key.isdigit() for key in ls_args):
            raise ValueError(f"The list of k sizes {arg} contains non integer values")
        if len(ls_args)==1:
            arg = int(arg)  
            if arg < 0:
                raise ValueError(f"{arg} is not a positive integer")
            if initial_args.class_type == 1:
                if not 8 <= arg <= 14 :
                    raise ValueError(f"A length of peptide of {arg} is not supported by netMHCpan (should be between 8 and 14)")
            else:
                if not arg >= 9 :
                    raise ValueError(f"A length of peptide of {arg} is not supported by netMHCIIpan (should not be less than 9)")

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

def check_hla_format(initial_args, parser,arg):
    """
    Returns the list of HLA genes if they are valid

    Parameters:
            initial_args (ArgumentParser object): parser object
            parser (CustomParser Object): parser object
            arg (str): list of HLA genes
    Returns:
            arg (str): list of HLA genes
    """
    if initial_args.class_type == 1:
        if not arg.isalnum() and not ("," in arg or "*" in arg):
            parser.error(f"{arg} contains non accepted characters")
    else:
        if not arg.isalnum() and not ("," in arg or "_" in arg):
            parser.error(f"{arg} contains non accepted characters")
    try:
        ls_args = arg.split(",")
        if len(ls_args) > 6:
            parser.error(f"The number of provided HLA class I genes is incorrect")
        if initial_args.class_type == 1:
            if not all(re.match(r"HLA\-[A-C]\*?\d{2}\:\d{2}",key) for key in ls_args):
                raise ValueError(f"The list {arg} does not comply to the expected format")
        else:
            if not all(re.match(r"D[P-R][A-B][1-9]\_?\d{2}:\d{2}",key) for key in ls_args):
                raise ValueError(f"The list {arg} does not comply to the expected format")
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"{err}")
    
    if initial_args.class_type == 1:        
        if "*" in arg:
            arg = "".join(arg.split("*"))
    else:
        if "_" in arg:
            arg = "".join(arg.split(":"))

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
    parser.add_argument("-n", "--run_name",
        help="name of the ams pipeline ran previously",
        action=arguments_handling.UniqueStore,
        required=True,
        type=lambda x: check_if_existing_run_name(parser,x))
    parser.add_argument("-p", "--pair",
        help="name of the pair",
        action=arguments_handling.UniqueStore,
        nargs="?",
        default="input_pair",
        const="input_pair",
        type=lambda x: arguments_handling.check_if_accepted_str(parser,x))
        
    # netMHCpan class
    parser.add_argument("-c", "--class_type",
        help="class type: class 1 (netMHCpan); class 2 (netMHCIIpan)",
        choices=[1, 2],
        default=1,
        type=int,
        )
    
    # intermediate arg parsing to consider classes 1 & 2
    initial_args, _ = parser.parse_known_args()
        
    # netMHCpan arg
    parser.add_argument("-l", "--length",
        help="peptide length, default is 9 (class 1) or 15 (class 2), a list of values separated by commas is accepted",
        nargs="?",
        default=9 if initial_args.class_type == 1 else 15,
        const=  9 if initial_args.class_type == 1 else 15,
        type=lambda x: check_if_valid_k(initial_args, x))
    # netMHCpan arg
    parser.add_argument("-a", "--hla_typing",
        help="comma separated list of HLA genes",
        required=True,
        type=lambda x: check_hla_format(initial_args, parser,x)
        )
    # netMHCpan arg
    parser.add_argument("-e","--el_rank",
        help=r"%%EL-rank filtration, all values above the given value are filtered out",
        default=100,
        type=lambda x: check_if_valid_float(parser,x)
        )
    args = parser.parse_args()
    return args
