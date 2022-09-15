#!/usr/bin/env python3.7

import argparse, sys, os
import logging
import shutil
import subprocess
import pandas as pd
from io import StringIO
from datetime import datetime
from tabulate import tabulate

inKSize = os.getenv('kSize')
print("The kmer size is exported from database using mash info: %s" % inKSize)

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

    required.add_argument("--database", "-b", required=True,
                        help="Pre-built Mash Sketch",
                        type=lambda x: parser.is_valid_mash(parser, x))
    required.add_argument("--read1", "-r1", required=True,
                        help="Input Read 1 (forward) file",
                        type=lambda x: parser.is_valid_fastq(parser, x))
    required.add_argument("--read2", "-r2", required=True,
                        help="Input Read 2 (reverse) file",
                        type=lambda x: parser.is_valid_fastq(parser, x))
    optional.add_argument("--max_dist", "-d", default=0.05,
                        help="User specified mash distance (default: 0.05)",
                        type=lambda x: parser.is_valid_distance(parser, x))
    optional.add_argument("--kmer_min", "-m", default=2,
                        help="Minimum copies of kmer count (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    optional.add_argument("--num_threads", "-p", default=2,
                        help="Number of computing threads to use (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    return parser

#inKSize = os.getenv('kSize')
#print("The kmer size is exported from database using mash info: %s" % inKSize)

###############
## FUNCTIONS ##
###############
def fastq_name(inRead1):
    """
    Gets stripped read name for appending to output files

    Parameters
    ----------
    inRead1 : xx

    Returns
    -------
    xxx
        xxx
    """
    if inRead1.endswith("_1.fastq"):
        name = inRead1.split("_1.fastq")[0]
        return(name)
    elif inRead1.endswith("_R1_001.fastq"):
        name = inRead1.split("_R1_001.fastq")[0]
        return(name)
    elif inRead1.endswith("_R1.fastq"):
        name = inRead1.split("_R1.fastq")[0]
        return(name)
    else:
        logging.critical("Please check your file endings, assumes either \
_1.fastq(.gz), _R1_001.fastq(.gz), or _R1.fastq(.gz)")
        sys.exit(1)

def make_output_log(log):
    """
    Makes the output directory and the log file

    Parameters
    ----------
    log : str
        Name of the log file

    Returns
    -------
    None
        Exits the program if unable to make output directory
    """
    logging.basicConfig(filename=log, filemode="a", level=logging.DEBUG,
    format="%(asctime)s - %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p")
    logging.info("New log file created in output directory - %s... " % log)
    logging.info("Starting the tool...")
    sysOutput=(os.uname())
    logging.info("Returning information identifying the current operating system... \n \n \
* System: %s  \n \
* Node Name: %s  \n \
* Release: %s  \n \
* Version: %s  \n \
* Machine %s \n" %
(sysOutput[0], sysOutput[1], sysOutput[2], sysOutput[3], sysOutput[4]))

def get_input(inRead1, inRead2, inMash, inMaxDis, inKmer, inKSize, inThreads):
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
 * Size of Kmer: %s \n \
 * Number of Threads: %s \n " %
 (inRead1, inRead2, inMash, inMaxDis, inKmer, inKSize, inThreads))

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
##TO DO:  - do i want to check if the beginning of the file name is a match between the two files?

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
    ##assumes that program name is lower case
    logging.info("Checking for program %s..." % program_name)
    path = shutil.which(program_name)
    ver = sys.version_info[0:3]
    ver  = ''.join(str(ver))
    ver = ver.replace(",", ".")
    ver = ver.replace('(','').replace(')','')

    if path != None:
        if program_name == 'python' and sys.version_info >= (3,7):
            logging.info("Great, the program %s is loaded ..." % program_name)
            logging.info("The version of python is: %s..." % ver)
        elif program_name != 'python':
            logging.info("Great, the program %s is loaded..." % program_name)
        else:
            logging.info("You do not have an appropriate version of python. \
 Requires Python version >= 3.7. Exiting.")
            sys.exit(1)
    else:
        logging.critical("Program %s not found! Cannot continue; dependency\
not fulfilled. Exiting." % program_name)
        sys.exit(1)


def check_mash(): 
    
    dirpath2 = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'test-data/'))

    filePath1 = os.path.join(dirpath2, 'myCatFile')
    filePath2 = os.path.join(dirpath2, 'myMashDatabase.msh')

    ## open(filepath)  # or whatever you use to open the file
    f1 = open(filePath1, 'r')
    f2 = open(filePath2, 'r')
    mashCheck = ['mash', 'dist', '-k', '25', '-s', '100000', filePath1, filePath2]
    
    result = subprocess.run(mashCheck, capture_output=True,\
        check=True, text=True)
    
    df = pd.read_csv(StringIO(result.stdout), sep='\t',
    names=['Ref ID', 'Query ID', 'Mash Dist', 'P-value', 'Kmer'],
    index_col=False)

    dfDropped = df.drop(['Ref ID', 'P-value', 'Kmer' ], axis=1)
    dfCheckSpecies = df['Query ID'].str.split('/')
    dfCheckSpecies = dfCheckSpecies[0] ##should be feeli
    dfCheckDist = dfDropped['Mash Dist'][0]
    dfCheckSpecies =(''.join(dfCheckSpecies))
    #logging.info("This is the sorted df from test: %s " % dfCheck2)
    if dfCheckSpecies == 'Legionella_fallonii_LLAP-10_GCA_000953135.1.fna' and int(dfCheckDist) == int(0.0185156) :
        logging.info("Great, the test to confirm Mash is running properly return our expected answers...")
        logging.info("This is what dfCheckSpecies is supposed to be: Legionella_fallonii_LLAP-10_GCA_000953135.1.fna")
        logging.info("This is what dfCheckSpecies returned: %s" % dfCheckSpecies)
        logging.info("This is what dfCheckDist is supposed to be: 0.0185156")
        logging.info("This is what dfCheckDist returned: %s" %dfCheckDist) 
    else:
        logging.info("This is what dfCheckSpecies is supposed to be: Legionella_fallonii_LLAP-10_GCA_000953135.1.fna")
        logging.info("This is what dfCheckSpecies returned: %s" % dfCheckSpecies)
        logging.info("This is what dfCheckDist is supposed to be: 0.0185156")
        logging.info("This is what dfCheckDist returned: %s" %dfCheckDist)
        logging.critical("The unit test to confirm species and mash value return a different answer than expected. Exiting")
        sys.exit(1)

def cat_files(inRead1, inRead2):
    """
    XXXX

    Parameters
    ----------
    inRead1 : XX
        XXX
    inRead2 : XX
        XXX

    Returns
    -------
    None
        XXX
    """
    if inRead1 and inRead2 != None and inRead1.endswith('.gz'):
        logging.critical("The files are still gzipped. Exiting")
        sys.exit(1)
    else:
        logging.info("The files have been gunzipped ...")
        with open(inRead1) as readFile:
            read1 = readFile.read()
        with open(inRead2) as readFile:
            read2 = readFile.read()
            read1 += read2
        with open('myCatFile', 'w') as readFile:
            readFile.write(read1)
            readFile.close()

def minKmer(calculatedKmer, inKmer):
    """
    Determine the value of the kmers (-m flag); if less than 2, set as 2

    Parameters
    ----------
    calculatedKmer : int
        Value that is calculated based on genomeCoverage/3
    inKmer : int, optional, default is 2
        Input kmer value specified by user; to be used instead of calucated

    Returns
    ----------
    int
        integer value used for min_kmer (-m flag) with paired-end reads
    """
    if int(inKmer) == 2:
        logging.info("Should kmer value be different than default (2)...")
        logging.info("Min. kmer = genome coverage divided by 3..." )
        #return calculatedKmer
        if (calculatedKmer < 2 or int(inKmer) < 2):
            logging.info("The calucated kmer is less than 2, so will use 2...")
            calculatedKmer = 2
        return calculatedKmer
    elif (int(inKmer) > 2):
        logging.info("User specified a value for minimum kmer: %s ..." % inKmer)
        return int(inKmer)

def run_cmd(command):
    """
    XXXX

    Parameters
    ----------
    XX : XX
        XXXcd

    Returns
    -------
    XXX : XXX
        XXX
    """

    try:
        result = subprocess.run(command, capture_output=True,\
        check=True, text=True)
        new_cmd=(' '.join(command))
        logging.info("This is the command... \n %s " % new_cmd)
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
    fastqCmd1 = ['mash', 'dist', inMash, '-r', 'myCatFile', '-p', inThreads, '-S', '42']

    outputFastq1 = run_cmd(fastqCmd1)

    ## get genome size and coverage; will provide as ouput for user
    gSize = outputFastq1.stderr.splitlines()[0]
    gSize = gSize[23:]
    logging.info("Estimated Genome Size to determine -m flag: %s " % gSize)
    gCoverage = outputFastq1.stderr.splitlines()[1]
    gCoverage = gCoverage[23:]
    logging.info("Estimated Genome coverage to determine -m flag: %s " % gCoverage)

    minKmers = int(float(gCoverage))/3
    minKmers = int(float(minKmers))

    ## this is used the calucate the minimum kmer copies to use (-m flag)
    mFlag = minKmer(minKmers, inKmer) # returned as an integer
    return mFlag, gSize, gCoverage

def get_results(mFlag, inThreads):
    fastqCmd2 = ['mash', 'dist', '-r', '-m', str(mFlag), inMash, 'myCatFile', '-p', inThreads, '-S', '123456']
    outputFastq2 = run_cmd(fastqCmd2)
    
    ## get genome size and coverage; will provide as ouput for user
    gSizeRun2 = outputFastq2.stderr.splitlines()[0]
    gSizeRun2 = gSizeRun2[23:]
    logging.info("Estimated Genome Size with the added -m flag: %s " % gSizeRun2)
    gCoverageRun2 = outputFastq2.stderr.splitlines()[1]
    gCoverageRun2 = gCoverageRun2[23:]
    logging.info("Estimated Genome coverage with the added -m flag: %s "% gCoverageRun2)
    #print(type(outputFastq2))
    #print(outputFastq2)

    return outputFastq2

def get_SC(outputFastq2):
	
	gSizeRun2 = outputFastq2.stderr.splitlines()[0]
	gSizeRun2 = gSizeRun2[23:]
	gCoverageRun2 = outputFastq2.stderr.splitlines()[1]
	gCoverageRun2 = gCoverageRun2[23:]
	return gSizeRun2, gCoverageRun2


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

    dfGenus = df['Ref ID'].str.split('_', 1, expand=True)
    tmpDF = pd.DataFrame(columns=['Genus', 'Species', 'GeneBank Identifier', '% Seq Sim'])
    tmpDF['Genus'] = dfGenus[0]

    dfSpecies = dfGenus.iloc[:, -1].str.split('_', 1, expand=True)
    tmpDF['Species'] = dfSpecies[0]

    dfGB = dfSpecies.iloc[:, -1].str.split('GCA', expand=True)
    dfGB = 'GCA' + dfGB.iloc[:,-1]
    tmpDF['GeneBank Identifier'] = dfGB

    df = df.join(tmpDF)

    ## split the kmers for sorting because xx/xxxx
    df[['KmersCount','sketchSize']] = df.Kmer.str.split("/", expand=True,)	
    df['KmersCount'] = df.KmersCount.astype(int)

    #df['Mash Dist'] = df['Mash Dist'].apply(lambda x: round(x,8)) 

    ## add column that is (1 - Mash Distance) * 100, which is % sequence similarity
    df['% Seq Sim'] =  1 - df['Mash Dist'] # multiplying by 100 gives strange result 
    df['% Seq Sim'] = df['% Seq Sim'] * 100    


    #logging.info("Should give dtype after adding to seq sim: %s" % df.dtypes)
    ## now sort and get top species; test for a tie in kmerscount value
    dfSorted = df.sort_values('KmersCount', ascending=False)
    dfSortOut = isTie(dfSorted)
    bestGenusSort = dfSortOut[0]
    bestSpeciesSort = dfSortOut[1]

    #pd.set_option('display.max_columns', 12)
    #pd.set_option('display.float_format', str)
    #use 7g to have it choose best viz or 7e to print scientific notation
    #pd.options.display.float_format = '{:.7g}'.format
    #print("Line 500", dfSorted)

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
    dfTop = dfSortedDropped[:5]

##TO DO - scienfitic notation for P-value

    dfTop.reset_index(drop=True, inplace=True) #make index start at 0
    return bestGenus, bestSpecies, dfTop

def makeTable(dateTime, name, inRead1, inRead2, inMaxDist, results, mFlag):
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

    with open(f"{name}_results_{dateString}.txt" ,'a+') as f:
        f.writelines("\n" + "Legionella Species ID Tool using Mash" + "\n")
        f.writelines("Date and Time = " + dtString + "\n") #+str(variable)
        f.write("Input query file 1: " + inRead1 + "\n")
        f.write("Input query file 2: " + inRead2 + "\n")
        f.write("Genome size estimate for fastq files with using the -m flag: " + scData[0] + " " +"(bp)" +"\n") #make into variable
        f.write("Genome coverage estimate for fastq files with using the -m flag: " + scData[1]  + "\n") #make into variables
        f.write("Maximum mash distance (-d): " + str(inMaxDis) + "\n")
        f.write("Minimum K-mer copy number (-m) to be included in the sketch: " + str(mFlag[0]) + "\n" )
        f.write("K-mer size used for sketching: " + inKSize + "\n" )
        f.write("Mash Database name: " + inMash + "\n" + "\n")
        f.write("Best species match: " + results[0] + " " + results[1] + "\n" + "\n")
        f.write("Top 5 hits:" + "\n")
        f.writelines(u'\u2500' * 100 + "\n")
        f.writelines(tabulate(results[2], headers='keys', tablefmt='pqsl', numalign="center", stralign="center", floatfmt=(None, None, None, ".5f", ".3f", ".8e"), showindex=False)+ "\n")

if __name__ == '__main__':
    ## parser is created from the function argparser
    ## parse the arguments
    parser = argparser()
    args = parser.parse_args()

inMash = args.database
inMaxDis = args.max_dist
inKmer = args.kmer_min
inThreads = args.num_threads
inRead1 = args.read1
inRead2 = args.read2

now = datetime.now()
dtString = now.strftime("%B %d, %Y %H:%M:%S")
dateString = now.strftime("%Y-%m-%d")

#unique name for log file based on read name
name = fastq_name(inRead1)
log = name + "_run"  + ".log"

req_programs=['mash', 'python']

make_output_log(log)
get_input(inRead1, inRead2, inMash, inMaxDis, inKmer, inKSize, inThreads)

#check_mash()

logging.info("Checking if all the required input files exist...")
check_files(inRead1, inRead2, inMash)
logging.info("Input files are present...")

logging.info("Checking if all the prerequisite programs are installed...")
for program in req_programs:
    check_program(program)
logging.info("All prerequisite programs are accessible...")

logging.info("Peforming internal system checks...")
check_mash()
logging.info("Great, internal system checks passed...")

logging.info("Begin concatenation of the fastq files...")
cat_files(inRead1, inRead2)
logging.info("Great, I was able to concatenate the files...")

logging.info("Calculating estimated genome size and coverage...")
mFlag = cal_kmer()
logging.info("Minimum copies of each kmer required to pass noise filter \
identified ...")

logging.info("Running Mash Dist command with -m flag...")
outputFastq2 = get_results(mFlag[0], inThreads)
logging.info("Completed running mash dist command...")

logging.info("Copying over estimated genome size and coverage to results...")
scData = get_SC(outputFastq2)
logging.info("Successfully, added estimated genome size and coverage to \
output...")

logging.info("Beginning to parse the output results from mash dist...")
results = parseResults(outputFastq2, inMaxDis)
logging.info("Okay, completed parsing of the results...")

logging.info("Generating table of results as a text file...")
makeTable(dtString, name, inRead1, inRead2, inMaxDis, results, mFlag)
logging.info("Completed analysis for the sample: %s..." % name )
