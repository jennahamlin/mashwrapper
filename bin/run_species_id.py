#!/usr/bin/env python

import argparse, sys, os
import logging
import shutil

#############################
## Argument Error Messages ##
#############################

## create class of formatted error messages
## input variable is argparse.ArgumentParser, which initializes the parser
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

    ## use class to parse the arguments with formatted error message

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
    optional.add_argument("--out_folder", "-o", default="out",
                        help="Output folder name (default: %(default)s)",
                         required=False)
    return parser


def configure_logger():
    """
    Configures the logger
    """
    try:
        logging.basicConfig(filename=log, filemode="w", level=logging.DEBUG,
                            format=f"%(asctime)s | Result folder {out_prefix}: %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p")
        logging.info(" Starting the tool...")
    except IOError:
        print(
            "I don't seem to have access to make the log file. Are the \
permissions correct or is there a directory with the same name?")
        sys.exit(1)

def make_output_directory():
    """
    Makes the output directory
    """
    global logging_message                                                      ## make global so visible outside function

    if os.path.isdir(out_prefix):                                               ## false if no output folder of same name
        print("Output directory - %s - exists. Please remove or rename the\
 directory. Exiting." % out_prefix)
        sys.exit(1)
    else:
        os.mkdir(out_prefix)
        logging_message += "New output directory created,\
 called: %s\n" % out_prefix


def check_program(program_name=None):
    """Checks if the supplied program_name exists
    Parameters
    ----------
    program_name : str
        Name of the program to check if it exists
    Returns
    -------
    None
        Exits the program if a dependency doesn't exist
    """
    logging.info(" Checking for program %s" % program_name)
    path = shutil.which(program_name)
    if path is None:
            logging.critical(" Program %s not found! Cannot continue; dependency not fulfilled." % program_name)
            print(" Looks like the program Mash is not available")
            #if not Inputs.verbose:
            #    print(f"Program {program_name} not found! Cannot continue; dependency not fulfilled.")
            sys.exit(1)
    else:
        logging.info(" Great, the program %s is loaded." % program_name)

def check_files() -> None:
    """Checks if all the input files exists; exits if file not found or if file is a directory
    Returns
    -------
    None
        Exits the program if file doesn't exist
    """
    if inRead1 and not os.path.isfile(inRead1):
        logging.critical("Read file 1: %s doesn't exist. Exiting" % inRead1)
        #if not Inputs.verbose:
        #    print(f"Read file 1: '{Inputs.read1}' doesn't exist. Exiting")
        sys.exit(1)
    if inRead2 and not os.path.isfile(inRead2):
        logging.critical("Read file 2: %s doesn't exist. Exiting" % inRead2)
        #if not Inputs.verbose:
        #    print(f"Read file 2: '{Inputs.read2}' doesn't exist. Exiting")
        sys.exit(1)
    if inMash and not os.path.isfile(inMash):
        logging.critical("Mash database file: %s doesn't exist. Exiting" %inMash)
        #if not Inputs.verbose:
        #    print(f"Assembly file: '{Inputs.assembly}' doesn't exist. Exiting")
        sys.exit(1)

if __name__ == '__main__':
    ## parser is created from the function argparser
    ## parse the arguments
    parser = argparser()
    args = parser.parse_args()

inMash = args.database
inMaxDist = args.max_dist
inKmer = args.min_kmer
inThreads = args.num_threads
inRead1 = args.read1
inRead2 = args.read2
out_prefix = args.out_folder
log = os.path.join(out_prefix, "run.log")
logging_message = " "
req_programs=['mash']

make_output_directory()
configure_logger()
for line in logging_message.rstrip().split("\n"):
    logging.info(line)

logging.info(" Checking if all the prerequisite programs are installed")
for program in req_programs:
    check_program(program)
logging.info(" All prerequisite programs are accessible")

logging.info(" Checking if all the required input files exist")
check_files()
logging.info(" Input files are present")


print(inMash)
print("This is read 1: ", inRead1)
print("This is read 2: ", inRead2)
print("This is the max distance: ", inMaxDist)
print("This is min_kmer :", inKmer)
print("This is default number of threads: ", inThreads)
