#!/usr/bin/env python

import argparse, sys, os
import logging
import shutil
import subprocess
import pandas as pd
from io import StringIO
from datetime import datetime
from tabulate import tabulate

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
    description = """ A script to parse the output from Mash into a table \
    listing the top five matches from the user specified pre-built Mash Database.
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
                        help="Minimum copies of kmer count (default: None)",
                        type=lambda x: parser.is_valid_int(parser, x))
    optional.add_argument("--num_threads", "-t",
                        help="Number of computing threads to use (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    optional.add_argument("--out_folder", "-o", default="resultsLog",
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
        logging.basicConfig(filename=log, filemode="a", level=logging.DEBUG,
                             format="%(asctime)s - %(message)s", \
                             datefmt="%m/%d/%Y %I:%M:%S %p")
        logging.info("New output directory created - %s... " % inResults)
        logging.info("New log file created in output directory - %s... " % log)
        logging.info("Starting the tool...")


def get_input(inRead1, inRead2, inMash, inMaxDis, inKmer, inThreads, inResults):
    """
    Prints the command line input to the log file

    Parameters
    ----------
    inRead1 : str
        XXX
    inRead2 : str
        XXX
    inMash : str
        XXX
    inMaxDis : str
        XXX
    inKmer : str
        XXX
    inThreads : str
        XXX
    inResults : str
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
    inRead1 : XX
        XXX
    inRead2 : XX
        XXX
    inMash : XX
        XXX
    Returns
    -------
    None
        Exits the program if file doesn't exist
    """

    if inMash and not os.path.isfile(inMash):
        logging.critical("The database - %s - doesn't exist. Exiting." % inMash)
        sys.exit(1)
    if inRead1 and not os.path.isfile(inRead1):
        logging.critical("Read file 1: %s doesn't exist. Exiting." % inRead1)
        sys.exit(1)
    if inRead2 and not os.path.isfile(inRead2):
        logging.critical("Read file 2: %s doesn't exist. Exiting." % inRead2)
        sys.exit(1)
    if inRead1 == inRead2:
        logging.critical("Read1 - %s" % inRead1)
        logging.critical("Read2 - %s" % inRead2)
        logging.critical("Looks like you entered the same read file twice. \
 Exiting.")
        sys.exit(1)
##TO DO:  - do i want to check if the beginning of the file name is a match?

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
            sys.exit(1)
    else:
        logging.info("Great, the program %s is loaded..." % program_name)

def cat_files(inResults, inRead1, inRead2):
    """
    XXXX

    Parameters
    ----------
    inResults : XX
        XXX
    inRead1 : XX
        XXX
    inRead2 : XX
        XXX

    Returns
    -------
    None
        XXX
    """

    ## change into the results folder
    print("Line 257 - path: ", os.getcwd())
    #os.chdir(inResults)
    print("Line 259 - New path: ", os.getcwd())
    #print(os.getcwd(inRead1))
    #gunzip inRead1
    #gunzip inRead2

    if inRead1 and inRead2 != None:
        with open(inRead1) as readFile:
            read1 = readFile.read()
        with open(inRead2) as readFile:
            read2 = readFile.read()
            read1 += read2
        with open('myCatFile', 'w') as readFile:
            readFile.write(read1)
            readFile.close()
    else:                                                                       ## this needs to be better message
        logging.critical("Hmm, I was unable to concatenate the files. Are the\
 permissions correct? Exiting.")
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
        logging.info("User specified a value for minimum kmer: %s " % inKmer)
        return inKmer
    elif (calculatedKmer < 2):
        logging.info("The calucated kmer is less than 2, so will use 2")
        calculatedKmer = 2
        return calculatedKmer
    else:
        logging.info("Min. kmer = genome coverage divided by 3..." )
        logging.info("This is the calcuated kmer: %s " % calculatedKmer)
        return calculatedKmer

def run_cmd(command):
    """
    XXXX

    Parameters
    ----------
    XX : XX
        XXX

    Returns
    -------
    XXX : XXX
        XXX
    """

    try:
        result = subprocess.run(command, capture_output=True,\
        check=True, text=True)
        logging.info("This is the command... \n %s " % command)
    except subprocess.CalledProcessError:
        logging.critical("CRITICAL ERROR. The following command had an improper\
 error: \n %s ." % command)
        sys.exit(1)
    return result

def cal_kmer():
    """
    XXXX

    Parameters
    ----------
    XX : XX
        XXX

    Returns
    -------
    mFlag : tuple, position 0
        XXX
    """

    f = open('myCatFile', 'r')
    fastqCmd1 = ['mash', 'dist', '-r', inMash, 'myCatFile']

    outputFastq1 = run_cmd(fastqCmd1)

    ## get genome size and coverage; will provide as ouput for user
    gSize = outputFastq1.stderr.splitlines()[0]
    gSize = gSize[23:]
    logging.info("Estimated Genome Size: %s " % gSize)
    gCoverage = outputFastq1.stderr.splitlines()[1]
    gCoverage = gCoverage[23:]
    logging.info("Estimated Genome coverage: %s "% gCoverage)

    minKmers = int(float(gCoverage))/3
    minKmers = int(float(minKmers))

    ## this is used the calucate the minimum kmer copies to use (-m flag)
    mFlag = minKmer(minKmers, inKmer) # returned as an integer
    return mFlag, gSize, gCoverage

def get_results(mFlag):
    fastqCmd2 = ['mash', 'dist', '-r', '-m', str(mFlag), inMash, 'myCatFile']
    outputFastq2 = run_cmd(fastqCmd2)
    os.remove('myCatFile')
    return outputFastq2

def isTie(df):
    """
    Determine if the kmers count value is a tie with the second top isolate; if
    so then indicate a tie was found in the results

    Parameters
    ----------
    df : pandas data frame

    Returns
    ----------
    str
        string of best genus either genus or a sentence
    str
        string of best species either species or a blank string

    """
    dfSort = df.sort_values('KmersCount', ascending=False)
    logging.info("Checking if matching kmers count is tied for top 2 results...")

## assumes based on position the first and second value are always
## column 11; row 1 and 2
    if dfSort.iloc[0:1, 10:11].equals(dfSort.iloc[1:2, 10:11]):
        bestGenus = "This was a tie, see the top 5 results below"
        bestSpecies = " "
        logging.info("The top two isolates have the same number of matching\
kmers, indicating a tie... ")
        return bestGenus, bestSpecies
    else:
        best = dfSort.head(1)
        bestGenus = best['Genus']
        bestGenus = bestGenus.str.cat(sep='\n')
        bestSpecies = best['Species']
        bestSpecies = bestSpecies.str.cat(sep='\n')
        logging.info("There was not a tie of kmers for the top two species...")
        return bestGenus, bestSpecies

def noResult(inFile, inMaxDis, bestG, bestS):
    """
    Determine if top hit mash distances are >= user specified mash distance

    Parameters
    ----------
    inFile : pandas core frame data frame
        Parsed results from running mash
    inDist : int
        User specified maximum mash distance as a cut-off

    Returns
    ----------
        Message to log file and replacement of best species with message
    """
    inMaxDis = float(inMaxDis)
    logging.info("Confirming that best match is less than user specfied distance...")

    if (inFile['Mash Dist'].values[0] < inMaxDis):
        logging.info("Okay, a best species match was found with mash distance \
less than %s..." % inMaxDis)
    else:
        bestG = "No matches found with mash distances < %s..." % inMaxDis
        bestS = " "
        logging.info("No matches found with mash distances < %s..." % inMaxDis)
    return bestG, bestS


def parseResults(cmd, inMaxDis):
    """
    run initial command and parse the results from mash

    Parameters
    ----------
    cmd : list
        Initial command to run for either fasta or fastq
    inMaxDis : XXX
        XXX

    Returns
    ----------
    bestGenus
        The most likely genus of the isolate tested
    bestSpecies
        The most likely species of the isolate tested
    dfTop
        The top five results from sorting Mash output
    """

    ## convert with StringIO and added headers (for development)
    df = pd.read_csv(StringIO(cmd.stdout), sep='\t',
    names=['Ref ID', 'Query ID', 'Mash Dist', 'P-value', 'Kmer'],
    index_col=False)

    ## maybe could use function whichFile?
    dfRef = df['Ref ID'].str.split('/', expand=True) #split the full path
    dfRef = dfRef.iloc[: , -1] #independent of full path, get final column
    dfRef = dfRef.str[:-12] #remove _cleaned.fna; assumes this is there; NOTE should probably exclude this in generating the mash sketch
    tmpDF = pd.DataFrame(columns=['Genus', 'Species', 'GeneBank Identifier', '% Seq Sim'])
    tmpDF[["Genus", "Species", "GeneBank Identifier"]] = dfRef.str.split("_",
    expand=True) #expand into 3 columns
    df = df.join(tmpDF)

    ## split the kmers for sorting because xx/xxxx
    df[['KmersCount','sketchSize']] = df.Kmer.str.split("/", expand=True,)
    df['KmersCount'] = df.KmersCount.astype(int)

    ## add column that is (1 - Mash Distance) * 100, which is % sequence similarity
    df['% Seq Sim'] =  (1 - df['Mash Dist'])*100

    ## now sort and get top species; test for a tie in kmerscount value
    dfSorted = df.sort_values('KmersCount', ascending=False)
    dfSortOut = isTie(dfSorted)
    bestGenusSort = dfSortOut[0]
    bestSpeciesSort = dfSortOut[1]

    ## use column (axis = 1), to create minimal dataframe
    dfSortedDropped = dfSorted.drop(['Ref ID', 'Query ID', 'KmersCount',
    'sketchSize' ], axis=1)

    ## noResult function - confirm mash distance is < than user specified
    ## even if mash distance !< user specified, return the top five hits
    noMash = noResult(dfSortedDropped, inMaxDis, bestGenusSort, bestSpeciesSort)
    bestGenus = noMash[0]
    bestSpecies = noMash[1]

    # change order
    dfSortedDropped = dfSortedDropped[['Genus', 'Species', 'GeneBank Identifier',
    'Mash Dist', '% Seq Sim', 'P-value', 'Kmer']]
    dfTop = dfSortedDropped[:200]

##TO DO - scienfitic notation for P-value

    dfTop.reset_index(drop=True, inplace=True) #make index start at 0
    return bestGenus, bestSpecies, dfTop

def makeTable(dateTime, inMaxDist, results, mFlag):
    """
    Parse results into text output and include relavant variables

    Parameters
    ----------
    dataTime : str
        get current date and time for when analysis is run
    maxDist : float, optional
        optional input to specify the value of maximum mash distance
    results : tuple
        output from running and parsing mash commands

    Returns
    ----------
    txt file
        text file with each isolates results appended that were run through

    """
    with open(f"Results_{dateString}.txt",'a+') as f:
        f.writelines("\n" + "\n" + "Legionella Species ID Tool using Mash" + "\n")
        f.writelines("Date and Time = " + dtString + "\n") #+str(variable)
        #f.write("Input query file(s):" + filesTested + "\n")
        f.write("Maximum mash distance: " + str(inMaxDis) + "\n")
        f.write("Genome size estimate for fastq files: " + mFlag[1] + " " +"(bp)" +"\n") #make into variable
        f.write("Genome coverage estimate for fastq files: " + mFlag[2]  + "\n") #make into variables
        f.write("Minimum kmer copy number to be included in the sketch:"  + "\n" + "\n")
        f.write("Best species match: " + results[0] + " " + results[1] + "\n" + "\n")
        f.write("Top 5 hits:" + "\n")
        f.writelines(u'\u2500' * 100 + "\n")
        f.writelines(tabulate(results[2], headers='keys', tablefmt='pqsl', numalign="center", stralign="center"))

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

now = datetime.now()
dtString = now.strftime("%B %d, %Y %H:%M:%S")
dateString = now.strftime("%Y-%m-%d")

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
logging.info("Great, I was able to concatenate the files...")

logging.info("Determining minimum kmer to use unless specified as input...")
mFlag = cal_kmer()
logging.info("Minimum kmer identified...")

logging.info("Running Mash Dist command with kmer...")
outputFastq2 = get_results(mFlag[0])
logging.info("Completed running mash dist command...")

logging.info("Beginning to parse the output results from mash dist...")
results = parseResults(outputFastq2, inMaxDis)
logging.info("Okay, completed parsing of the results...")

logging.info("Generating table of results as a text file...")
makeTable(dtString, inMaxDis, results, mFlag)
logging.info("Completed analysis for this isolate...")

## TO DO: Move the function to make folder and run log outside of this script
## that should allow things to be appended and not re-wrritten
## TO DO: Check for .gz files, if found then gunzip them before cat_files function
