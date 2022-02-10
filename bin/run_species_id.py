#!/usr/bin/env python

import argparse, sys, os
import logging
import shutil

#############################
## Argument Error Messages ##
#############################
out_prefix = "out"
log = os.path.join(out_prefix, "run.log")
logging_message = " "
## preset some arguments
# class Inputs:
# out_prefix = "out"
# log = os.path.join(out_prefix, "run.log")
# logging_message = " "

class ParserWithErrors(argparse.ArgumentParser):
    def error(self, message):
        print('\n{0}\n\n'.format(message))
        self.print_help()
        sys.exit(2)

    def is_valid_mash(self, parser, arg):
        base, ext = os.path.splitext(arg)
        if ext not in ('.msh') or ext in ('') :
            parser.error('This is not a file ending with .msh.\
 Did you generate the mash sketch and specified that file to be uploaded?')
        else:
            return arg

    def is_valid_fastq(self, parser, arg):
        base, ext = os.path.splitext(arg)
        if ext not in ('.gz', '.fastq', '.fastq.gz'):
            parser.error('This is not a file ending with either .fastq or \
.fastq.gz. This flag requires the input of a fastq file.')
        else:
            return arg

    def is_valid_distance(self, parser, arg):
        isFloat = True
        try:
            float(arg)
        except ValueError:
            isFloat = False
        if (isFloat == False) or (float(arg) < 0) :
            parser.error('%s is not a positive number (e.g., a float, aka a \
 number with a decimal point)' % arg)
        else:
            return arg

    def is_valid_int(self, parser, arg):
        isInt = True
        try:
            int(arg)
        except ValueError:
            isInt = False
        if isInt == False:
            parser.error("You input %s. This is NOT an integer." % arg)
        elif arg.isnumeric() == False:
            parser.error("You input %s. This is NOT a positive number." % arg)
        else:
            return arg

#########################
## ArgParser Arguments ##
#########################

def argparser():
    """
    Returns an argument parser for this script
    """
    description = """
    A script to parse the output from Mash and sorts the output to give you the\
    most similar species/isolate determined using various thresholds.
    """

    parser = ParserWithErrors(description = description)

    ## Define required and optional groups
    parser._action_groups.pop()
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')

    required.add_argument("--database", "-d", required=True,
                        help="Pre-built Mash Sketch",
                        type=lambda x: parser.is_valid_mash(parser, x))
    required.add_argument("--read1", "-r1", required=True,
                        help="Input Read 1 (forward) file",
                        type=lambda x: parser.is_valid_fastq(parser, x))
    required.add_argument("--read2", "-r2", required=True,
                        help="Input Read 2 (reverse) file",
                        type=lambda x: parser.is_valid_fastq(parser, x))
    optional.add_argument("--max_dist", "-m", default=0.05,
                        help="User specified mash distance (default: 0.05)",
                        type=lambda x: parser.is_valid_distance(parser, x))
    optional.add_argument("--min_kmer", "-k", default=2,
                        help="Minimum copies of kmer count to use (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    optional.add_argument("--num_threads", "-t", default=2,
                        help="Number of computing threads to use (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    optional.add_argument("--out", "-o", default=out_prefix,
                        help="Output folder name (default: %(default)s)",
                         required=False)
    return parser

def make_output_directory(arg):
    """
    Makes the output directory
    """
    global logging_message
    #print(os.path.isdir(out_prefix))                                           ## should say false if never created
    if os.path.isdir(arg):
        print("Output directory - %s - exists. Please remove or rename the\
 directory. Exiting." % arg)
        sys.exit(1)
    else:
        os.mkdir(out_prefix)
        logging_message += "New output directory created\n"

def configure_logger():
    try:
        print(logging.basicConfig(filename=log, filemode="w", level=logging.DEBUG,
                            format=f"[%(asctime)s | {out_prefix} ]  %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p"))
    except FileNotFoundError:
        print(
            f"The supplied location for the log file '{log}' doesn't exist. Please check if the location exists.")
        sys.exit(1)
    except IOError:
        print(
            f"I don't seem to have access to make the log file. Are the permissions correct or is there a directory with the same name?")
        sys.exit(1)


if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

log = os.path.join(out_prefix, "run.log")
logging_message = " "

inMash = args.database
inMaxDist = args.max_dist
inKmer = args.min_kmer
inThreads = args.num_threads
inRead1 = args.read1
inRead2 = args.read2
inOut = args.out

#set_inputs()
make_output_directory(inOut)
configure_logger()
logging.info("Starting preprocessing")
for line in logging_message.rstrip().split("\n"):
    logging.info(line)

print(inMash)
print("This is read 1: ", inRead1)
print("This is read 2: ", inRead2)
print("This is the max distance: ", inMaxDist)
print("This is min_kmer :", inKmer)
print("This is default number of threads: ", inThreads)
