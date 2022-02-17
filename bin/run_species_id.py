#!/usr/bin/env python

import argparse, sys, os
import logging
import shutil
import subprocess

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
    optional.add_argument("--min_kmer", "-k",
                        help="Minimum copies of kmer count to use (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    optional.add_argument("--num_threads", "-t",
                        help="Number of computing threads to use (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    optional.add_argument("--out_folder", "-o", default="out",
                        help="Output folder name (default: %(default)s)",
                         required=False)
    return parser

def make_output_w_log(log, inResults):
    """
    Makes the output directory and the log file

    Parameters
    ----------
    log : str
        Name of the log file
    inResults : str
        Name of the output folder where the log file and results are stored

    Returns
    -------
    None
        Exits the program if unable to make output directory
    """

    if os.path.isdir(inResults):                                                ## false if no output folder of same name
        print("Output directory - %s - exists. Please remove or rename the\
 directory. Exiting." % inResults)
        sys.exit(1)
    else:
        os.mkdir(inResults)
        logging.basicConfig(filename=log, filemode="w", level=logging.DEBUG,
                             format="%(asctime)s - %(message)s", \
                             datefmt="%m/%d/%Y %I:%M:%S %p")
        logging.info("Starting the tool...")
        logging.info("New output directory created - %s... " % inResults)
        logging.info("New log file created in output directory - %s... " % log)

def get_input(inRead1, inRead2, inMash, inMaxDis, inKmer, inThreads, inResults):
    """
    XXXX

    Parameters
    ----------
    XX : XX
        XXX

    Returns
    -------
    None
        XXX
    """

    logging.info("The input parameters... \n \n \
 * Read1: %s \n \
 * Read2: %s \n \
 * Mash Databse: %s \n \
 * Maximum Distance: %s \n \
 * Minimum Kmer Count: %s \n \
 * Number of Threads: %s \n \
 * Result Folder: %s \n "  % \
 (inRead1, inRead2, inMash, inMaxDis, inKmer, inThreads,\
  os.path.abspath(inResults)))

def check_files(inRead1, inRead2, inMash):
    """
    Checks if all the input files exists; exits if file not found or if file is
    a directory

    Parameters
    ----------
    XX : XX
        XXX

    Returns
    -------
    None
        Exits the program if file doesn't exist
    """

    if inMash and not os.path.isfile(inMash):
        logging.critical("The database - %s - doesn't exist. Exiting." % inMash)
        print("Your mash database - %s - does not exist. Exiting." % inMash)
        #if not Inputs.verbose:
        #    print(f"Assembly file: '{Inputs.assembly}' doesn't exist. Exiting")
        sys.exit(1)
    if inRead1 and not os.path.isfile(inRead1):
        logging.critical("Read file 1: %s doesn't exist. Exiting." % inRead1)
        print("Your read 1 file - %s - does not exist. Exiting." % inRead2)
        #if not Inputs.verbose:
        #    print(f"Read file 1: '{Inputs.read1}' doesn't exist. Exiting")
        sys.exit(1)
    if inRead2 and not os.path.isfile(inRead2):
        logging.critical("Read file 2: %s doesn't exist. Exiting." % inRead2)
        print("Your read 2 file - %s - does not exist. Exiting." % inRead2)
        #if not Inputs.verbose:
        #    print(f"Read file 2: '{Inputs.read2}' doesn't exist. Exiting")
        sys.exit(1)
    if inRead1 == inRead2:
        logging.critical("Read1 - %s" % inRead1)
        logging.critical("Read2 - %s" % inRead2)
        logging.critical("Looks like you entered the same read file twice. \
 Exiting.")
        print("It appears the read names are not unique, please enter Read1 \
(R1) and Read2 (R2). Exiting")
        sys.exit(1)

def check_program(program_name):
    """
    Checks if the supplied program_name exists

    Parameters
    ----------
    program_name : str
        Name of the program to check if it exists

    Returns
    -------
    None
        Exits the program if a dependency doesn't exist
    """

    logging.info("Checking for program %s..." % program_name)
    path = shutil.which(program_name)                                           ## assumes that programs are lower case
    if path is None:
            logging.critical("Program %s not found! Cannot continue; dependency\
 not fulfilled. Exiting." % program_name)
            print("Looks like the program Mash is not available. Exiting.")
            #if not Inputs.verbose:
            #    print(f"Program {program_name} not found! Cannot continue; dependency not fulfilled.")
            sys.exit(1)
    else:
        logging.info("Great, the program %s is loaded..." % program_name)

def cat_files(inResults, inRead1, inRead2):
    """
    XXXX

    Parameters
    ----------
    XX : XX
        XXX

    Returns
    -------
    None
        XXX
    """

    ## change into the results folder
    os.chdir(inResults)

    if inRead1 and inRead2 != None:
        with open(inRead1) as readFile:
            read1 = readFile.read()
        with open(inRead2) as readFile:
            read2 = readFile.read()
            read1 += read2
        with open('myCatFile', 'w') as readFile:
            readFile.write(read1)
            readFile.close()
        logging.info("Great, I was able to concatenate the files...")
    else:                                                                       ## this needs to be better message
        logging.critical("Hmm, I was unable to concatnate the files. Are the\
 permissions correct? Exiting.")
        print("Looks like I could not concatenate the files. Exiting.")
        sys.exit(1)

def minKmer(calculatedKmer, inKmer):
    """
    Determine the value of the kmers (-m flag); if less than 2, set as 2

    Parameters
    ----------
    calculatedKmer : int
        Value that is calculated based on genomeCoverage/3
    min_kmer : int, optional
        Input kmer value specified by user; to be used instead of calucated

    Returns
    ----------
    int
        integer value used for min_kmer (-m flag) with paired-end reads
    """

    if inKmer != None:
        logging.info("Okay, user specified a value for minimum kmer: %s " % inKmer)
        return inKmer
    elif (calculatedKmer < 2):
        logging.info("The calucated kmer is less than 2, so will use 2")
        calculatedKmer = 2
        return calculatedKmer
    else:
        logging.info("This is the calcuated kmer: %s " % calculatedKmer)
        return calculatedKmer

def cal_kmer():
    """
    XXXX

    Parameters
    ----------
    XX : XX
        XXX

    Returns
    -------
    None
        XXX
    """

    f = open('myCatFile', 'r')
    fastqCmd = ['mash', 'dist', '-r', inMash, 'myCatFile']

    outputFastq = subprocess.run(fastqCmd, capture_output=True,
    check=True, text=True)

    ## get genome size and coverage; will provide as ouput for user
    gSize = outputFastq.stderr.splitlines()[0]
    gSize = gSize[23:]
    logging.info("Genome Size: %s " % gSize)
    gCoverage = outputFastq.stderr.splitlines()[1]
    gCoverage = gCoverage[23:]
    logging.info("Genome coverage: %s "% gCoverage)

    minKmers = int(float(gCoverage))/3
    minKmers = int(float(minKmers))
    logging.info("Min. kmer is genome coverage divided by 3 = %s " % minKmers)

    ## this is used the calucate the minimum kmer copies to use (-m flag)
    mFlag = minKmer(minKmers, inKmer ) # returned as an integer



if __name__ == '__main__':
    ## parser is created from the function argparser
    ## parse the arguments
    parser = argparser()
    args = parser.parse_args()

inMash = args.database
inMaxDis = args.max_dist
inKmer = args.min_kmer
inThreads = args.num_threads
inRead1 = args.read1
inRead2 = args.read2
inResults= args.out_folder
log = os.path.join(inResults, "run.log")
req_programs="mash"

make_output_w_log(log, inResults)

get_input(inRead1, inRead2, inMash, inMaxDis, inKmer, inThreads, inResults)

logging.info("Checking if all the required input files exist...")
check_files(inRead1, inRead2, inMash)
logging.info("Input files are present...")

logging.info("Checking if all the prerequisite programs are installed...")
check_program(req_programs)
logging.info("All prerequisite programs are accessible...")

logging.info("First concatenating the fastq files...")
cat_files(inResults, inRead1, inRead2)

logging.info("Determining minimum kmer to use unless specified as input...")
cal_kmer()
logging.info("Completed the calculation. ")
