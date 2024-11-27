# coding:utf-8
"""
This script was developed with the argparse module to make 
the ams pipeline command line easy to use and user friendly
"""
import argparse
import sys
import os

# Custom Parser


class CustomParser(argparse.ArgumentParser):
    """
    A class to handle the errors of the parser
    Methods :
            error :
                    Returns the error message, prints the help and exits with exit status 2
    """

    def error(self, message):
        sys.stderr.write(f"error: {message}\n")
        self.print_help()
        sys.exit(2)


# Custom Actions


class UniqueStore(argparse.Action):
    """
    A class used to make sure there is a unique use of the optional argument
    Methods :
        __call__ :
            Returns the value set in the command line or raises an error if the argument was already called
    """

    def __call__(self, parser, namespace, values, option_string):
        if getattr(namespace, self.dest, self.default) is not self.default:
            parser.error(option_string + " appears several times.")
        setattr(namespace, self.dest, values)


# Error Handling
def check_file(parser, arg, wc=False, csv=False):
    """
    Returns the file name after checking if it exists and if the format is accepted

    Parameters:
                    parser (CustomParser Object): parser object
                    arg (str): file name
                    wc (bool): boolean that is True if the worst consequences are toggled
                    csv (bool): boolean that is True if the expected arg is a CSV file
    Returns:
                    arg (str): file name
    """
    try:
        file = open(arg, "r")
        file.close()
    except OSError as err:
        parser.error(f"OS error: {err}")
    else:
        # check file extension
        if csv:
            if not arg.endswith("csv") and not wc:
             parser.error(
                    f"{arg} has an invalid extension :\nsupported file extension : csv"
                )
            else:
                return arg
        else:
            if not arg.endswith("vcf") and not arg.endswith("vcf.gz") and not wc:
                parser.error(
                    f"{arg} has an invalid extension :\nsupported file extensions : vcf, vcf.gz"
                )
            else:
                return arg
    return None


def check_if_accepted_str(parser, arg):
    """
    Returns the run name after checking if it is valid

    Parameters:
                    parser (CustomParser Object): parser object
                    arg (str): run name
    Returns:
                    arg (str): file name
    """
    if not arg.isalnum() and not ("-" in arg or "_" in arg):
        parser.error(f"{arg} contains non accepted characters")
    return arg


# check if type of arg is int
def handle_type_error(arg):
    """
    Returns a positive integer after checking if it is one

    Parameters:
                    arg (str): int
    Returns:
                    arg (str): int
    """
    try:
        arg = int(arg)
        if arg < 0:
            raise ValueError(f"{arg} is not a positive integer")
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"{err}")
    return arg


def check_threshold_value(parser, arg, arg_name):
    """
    Returns the homozygosity threshold value after checking if it is a float between 0 and 1

    Parameters:
                    arg (str): float
    Returns:
                    arg (str): float
    """
    try:
        threshold = float(arg)
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"{err}")
    if arg_name == "gnomad_af":
        condition = 0 <= threshold < 1
        error_message = f"{threshold} not between 0 and 1 (0 included)"
    elif arg_name == "homozygosity_thr":
        condition = 0 < threshold < 1
        error_message = f"{threshold} not between 0 and 1 (0 excluded)"
    else:
        parser.error(f"Unknown argument: {arg}")
    if not condition:
            parser.error(error_message)
    return threshold
    
    
def check_workers_count(parser,arg):
    cpu_count = os.cpu_count()
    if cpu_count is None:
        raise ValueError("Failed to get CPU count!")
    try:
       workers = int(arg)
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"{err}")
    if not 1 <= workers <= cpu_count:
        parser.error(f"{workers} is not in available number of cores (1-{cpu_count})")
    return workers


def arguments(from_filePair : bool = False):
    """
    Returns the parsed arguments of the ams pipeline
    Returns:
                    args (argparse.Namespace): object containing all parameters
    """
    if from_filePair:
        parser = CustomParser(
            prog="multiprocess_ams.py",
            usage="python %(prog)s [options] file_pairs.csv orientation",
            description="Compute the AMS for a list of pairs of individuals",
            )
        parser.add_argument(
            "multi_vcf",
            help="multi-sample VCF file, accepted formats are vcf and vcf.gz",
            type=lambda x: check_file(parser, x)
            )
        parser.add_argument(
            "file_pairs",
            help="file with pairs of donors and recipients, accepted format is csv",
            type=lambda x: check_file(parser, x, False, True)
        ) 
    else:
        parser = CustomParser(
            prog="ams_pipeline.py",
            usage="python %(prog)s [options] donor recipient orientation",
            description="Compute the AMS for a pair of individuals",
        )
        parser.add_argument(
            "donor",
            help="donor file, accepted formats are vcf and vcf.gz",
            type=lambda x: check_file(parser, x),
        )
        parser.add_argument(
            "recipient",
            help="recipient file, accepted formats are vcf and vcf.gz",
            type=lambda x: check_file(parser, x),
        )
    parser.add_argument(
        "orientation",
        help="choose the orientation of the comparison of the pair",
        choices=["dr", "rd"],
    )
    parser.add_argument(
        "imputation",
        help="choose the imputation mode",
        choices=["imputation", "no-imputation"],
    )
    parser.add_argument(
        "--min_dp",
        help="minimal accepted depth per position",
        nargs="?",
        default=20,
        const=20,
        type=lambda x: handle_type_error(x),
    )
    parser.add_argument(
        "--max_dp",
        help="maximal accepted depth per position",
        nargs="?",
        default=400,
        const=400,
        type=lambda x: handle_type_error(x),
    )
    parser.add_argument(
        "--min_ad",
        help="minimal value for read counts per position",
        nargs="?",
        default=5,
        const=5,
        type=lambda x: handle_type_error(x),
    )
    parser.add_argument(
        "-t",
        "--homozygosity_thr",
        help="allelic ratio threshold for which a heterozygous position can be converted to homozygous",
        nargs="?",
        default=0.8,
        const=0.8,
        type=lambda x: check_threshold_value(parser, x, "homozygosity_thr"),
    )
    parser.add_argument("--gnomad_af",
        help="SNPs with AF in combined population below this value will be filtered",
        nargs="?",
        default=0.01,
        const=0.01,
        type=lambda x: check_threshold_value(parser,x, "gnomad_af")
    )
    parser.add_argument(
        "--min_gq",
        help="genotype quality, the higher the more reliable the predicted genotype",
        nargs="?",
        default=20,
        const=20,
        type=lambda x: handle_type_error(x),
    )
    parser.add_argument(
        "-l",
        "--base_length",
        help="maximal base length tolerated in REF and ALT fields",
        action=UniqueStore,
        nargs="?",
        default=3,
        const=3,
        type=lambda x: handle_type_error(x),
    )
    parser.add_argument(
        "-n",
        "--run_name",
        help="name of the run",
        action=UniqueStore,
        nargs="?",
        default="run",
        const="run",
        type=lambda x: check_if_accepted_str(parser, x),
    )
    parser.add_argument(
        "-p",
        "--pair",
        help="name of the pair",
        action=UniqueStore,
        nargs="?",
        default="input_pair",
        const="input_pair",
        type=lambda x: check_if_accepted_str(parser, x),
    )
    parser.add_argument(
        "-ns",
        "--norm_score",
        help="toggle score normalization (recommended for multiprocess_ams only)",
        action="store_true",
    )
    parser.add_argument(
        "-wc",
        help="toggle worst consequence annotations from Variant Effect Predictor",
        action="store_true",
    )
    parser.add_argument(
        "-wcd",
        "--wc_donor",
        help="donor file of the worst consequences per position predicted by Variant Effect Predictor",
        metavar="<wc_donor>",
        action=UniqueStore,
        required="-wc" in sys.argv,
        type=lambda x: check_file(parser, x, True),
    )
    parser.add_argument(
        "-wcr",
        "--wc_recipient",
        help="recipient file of the worst consequences per position predicted by Variant Effect Predictor",
        metavar="<wc_recipient>",
        action=UniqueStore,
        required="-wc" in sys.argv,
        type=lambda x: check_file(parser, x, True),
    )
    parser.add_argument(
        "-w",
        "--workers",
        help="number of workers (cores) for multiprocessing",
        default=os.cpu_count() // 2,
        type=lambda x: check_workers_count(parser, x)
    )
    args = parser.parse_args()
    if args.min_dp > args.max_dp:
        raise ValueError(
            "The minimal Depth should not be higher than the maximal Depth"
        )
    if args.min_ad > args.min_dp:
        raise ValueError(
            "The minimal Allelic Depth should not be higher than the minimal Depth"
        )
    return args
