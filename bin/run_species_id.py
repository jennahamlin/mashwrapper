#!/usr/bin/env python3.7

import argparse, sys, os
import logging
import re
import shutil
import subprocess
import pandas as pd
from io import StringIO
from datetime import datetime
from tabulate import tabulate

#inKSize = os.getenv('kSize')
#print("The kmer size is exported from database using mash info: %s" % inKSize)

#############################
## Argument Error Messages ##
#############################

class ParserWithErrors(argparse.ArgumentParser):
    def error(self, message):
        print(f'\nError: {message}\n')
        self.print_help()
        sys.exit(2)

    def is_valid_mash(self, parser, arg):
        base, ext = os.path.splitext(arg)
        if ext != '.msh':
            parser.error(f'File {arg} does not have a .msh extension. Did you generate the mash sketch and specify this file?')
        return arg

    def is_valid_fastq(self, parser, arg):
        base, ext = os.path.splitext(arg)
        if ext not in ('.gz', '.fastq', '.fq', '.fastq.gz'):
            parser.error(f'File {arg} does not have a valid fastq extension (.fq, .gz, .fastq, .fastq.gz).')
        return arg

    def is_valid_distance(self, parser, arg):
        try:
            distance = float(arg)
            if distance < 0:
                raise argparse.ArgumentTypeError(f'Distance must be a positive number. Received: {distance}')
            return arg
        except ValueError:
            parser.error(f'{arg} is not a valid float number (e.g., 0.1).')

    def is_valid_int(self, parser, arg):
        try:
            value = int(arg)
            if value <= 0:
                raise argparse.ArgumentTypeError(f'Integer must be a positive number. Received: {value}')
            return arg
        except ValueError:
            parser.error(f'{arg} is not a valid integer.')

#########################
## ArgParser Arguments ##
#########################

def argparser():
    """
    Returns argument parser for the script with messages for how to use tool.
    """

    description = "A script to run and parse the output from Mash into a table listing the top five matches from the user specified pre-built Mash Database."

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

###############
## FUNCTIONS ##
###############
def make_output_log(log):
    """
    Makes the output directory and the log file which can be appended too.
    Includes printing of operating system (os) information where script is being run.
    That os information is determined positioanlly (e.g., sysOutput[0]). Requires
    logging package.

    Parameters
    ----------
    log : str
        Name of the log file.

    Returns
    -------
    None
        Exits the program if unable to make output directory.
    """
    # Configure logging
    logging.basicConfig(filename=log,
                        filemode="a",
                        level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt="%m/%d/%Y %I:%M:%S %p")

    # Log file creation message
    logging.info("New log file created in output directory - %s" % log)
    logging.info("Starting the tool...")

    # Log system information
    sys_info = os.uname()
    logging.info("System Information:")
    logging.info("   System: %s" % sys_info[0])
    logging.info("   Node Name: %s" % sys_info[1])
    logging.info("   Release: %s" % sys_info[2])
    logging.info("   Version: %s" % sys_info[3])
    logging.info("   Machine: %s" % sys_info[4])

def fastq_name(read1: str):
    """
    Gets stripped read name for appending to output files. 
    When running within NextFlow, allows for .gz files because files are
    gunzipped before running this script.

    Parameters
    ----------
    read1 (str) : Fastq files with the option of three different file endings.

    Returns
    -------
    name (str): File names with one of the three supported endings stripped.
    """
    # Define regex pattern to match any of the supported endings
    pattern = r'(_1|_R1_001|_R1)\.fastq(\.gz)?$'

    match = re.search(pattern, read1)
    if match:
        name = read1[:match.start()]  # Extract everything before the matched pattern
        return name
    else:
        logging.critical(" Please check your file endings, assumes either _1.fastq(.gz), _R1_001.fastq(.gz), or _R1.fastq(.gz)")
        sys.exit(1)

##TODO reduce the number of variables here. Split to two funcitons on for required varaibles and one for optional variables

def get_input(read1: str, read2: str, mash_db: str, max_dis: str, min_kmer: str, k_size: str, threads: str) -> None:
    """
    Logs the command line input parameters.

    Parameters
    ----------
    read1 : str
        Path to the first input read.
    read2 : str
        Path to the second input read.
    mash_db : str
        Path to the Mash database.
    max_dis : str
        Maximum distance parameter.
    min_kmer : str
        Minimum k-mer count.
    k_size : str
        Size of the k-mer.
    threads : str
        Number of threads to use.

    Returns
    -------
    None
    """
    
    logging.info(
        f"The input parameters:\n"
        f" * Read1: {read1}\n"
        f" * Read2: {read2}\n"
        f" * Mash Database: {mash_db}\n"
        f" * Maximum Distance: {max_dis}\n"
        f" * Minimum Kmer Count: {min_kmer}\n"
        f" * Size of Kmer: {k_size}\n"
        f" * Number of Threads: {threads}\n"
    )

##TO DO - do we want to check for corrupt gzip files?
##TO DO - do i want to check if the beginning of the file name is a match between the two files?

def get_k_size(mash_db: str) -> str:
    """
    Retrieves the k-size from the Mash database information.

    Parameters
    ----------
    mash_db : str
        Path to the Mash database.

    Returns
    -------
    str
        The k-size value extracted from the Mash info output.
    """
    try:
        # Run the mash info command and capture its output
        result = subprocess.run(
            ['mash', 'info', mash_db],
            capture_output=True,
            text=True,
            check=True
        )

        # Extract the k-size value using awk from the command's output
        output = result.stdout
        lines = output.splitlines()
        
        if len(lines) >= 3:
            # Assuming the k-size is in the 3rd line and 3rd field
            k_size = lines[2].split()[2]  
            return k_size
        else:
            raise ValueError("Unexpected output format from mash info command.")
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running the command: {e}")
        return ""
    except Exception as e:
        print(f"An error occurred: {e}")
        return ""

def check_files(read1, read2, inMash):
    """
    Checks if all the input files exist; exits if file not found or if file is
    a directory.

    Parameters
    ----------
    read1 : str or None
        Path to input file 1.
    read2 : str or None
        Path to input file 2.
    inMash : str or None
        Path to database file.

    Returns
    -------
    None
        Exits the program if any file doesn't exist.
    """

    def check_file(path, description):
        if path and not os.path.isfile(path):
            logging.critical("%s doesn't exist or is not a file: %s. Exiting." % (description, path))
            sys.exit(1)

    check_file(inMash, "The database file")
    check_file(read1, "Read file 1")
    check_file(read2, "Read file 2")

    if read1 == read2:
        logging.critical("Read1 (%s) and Read2 (%s) are the same file. Exiting." % (read1, read2))
        sys.exit(1)

def check_program(program_name: str) -> None:
    """
    Checks if the supplied program_name exists and if it's an appropriate version of Python.

    Parameters
    ----------
    program_name : str
        Name of the program to check if it exists.

    Returns
    -------
    None
        Exits the program if a dependency doesn't exist.
    """
    logging.info(f"Checking for program {program_name}...")

    # Check if the program exists
    path = shutil.which(program_name)
    
    if path is None:
        logging.critical(f"Program {program_name} not found! Cannot continue; dependency not fulfilled. Exiting.")
        sys.exit(1)

    # If the program is Python, check the version
    if program_name == 'python':
        python_version = '.'.join(map(str, sys.version_info[:3]))
        if sys.version_info >= (3, 7):
            logging.info(f"Great, the program {program_name} is loaded.")
            logging.info(f"The version of python is: {python_version}.")
        else:
            logging.critical("You do not have an appropriate version of Python. Requires Python version >= 3.7. Exiting.")
            sys.exit(1)
    else:
        logging.info(f"Great, the program {program_name} is loaded.")

def check_mash() -> None:
    """
    Checks the output of the Mash command and verifies the results.

    Returns
    -------
    None
        Exits the program if the results do not match the expected values.
    """
    # Define paths
    dirpath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'test-data'))
    file_path1 = os.path.join(dirpath, 'myCatFile')
    file_path2 = os.path.join(dirpath, 'myMashDatabase.msh')

    # Define Mash command
    mash_check = ['mash', 'dist', '-k', '25', '-s', '100000', file_path1, file_path2]

    # Run the Mash command
    result = subprocess.run(mash_check, capture_output=True, check=True, text=True)
    
    # Read and process the output
    df = pd.read_csv(StringIO(result.stdout), sep='\t', names=['Ref ID', 'Query ID', 'Mash Dist', 'P-value', 'Kmer'])
    df_dropped = df.drop(columns=['Ref ID', 'P-value', 'Kmer'])
    
    # Extract relevant information
    df_check_species = df['Query ID'].str.split('/').str[0].iloc[0]
    df_check_dist = df_dropped['Mash Dist'].iloc[0]
    
    expected_species = 'Legionella_fallonii_LLAP-10_GCA_000953135.1.fna'
    expected_dist = 0.0185
    
    # Log and validate the results
    if df_check_species == expected_species and float(round(df_check_dist, 4)) == expected_dist:
        logging.info("Great, the test confirms Mash is running properly and returned expected results.")
        logging.info(f"Expected species: {expected_species}")
        logging.info(f"Returned species: {df_check_species}")
        logging.info(f"Expected distance: {expected_dist}")
        logging.info(f"Returned distance: {round(df_check_dist, 4)}")
    else:
        logging.info(f"Expected species: {expected_species}")
        logging.info(f"Returned species: {df_check_species}")
        logging.info(f"Expected distance: {expected_dist}")
        logging.info(f"Returned distance: {df_check_dist}")
        logging.critical("The unit test to confirm species and mash value did not return expected results. Exiting.")
        sys.exit(1)


def cat_files(read1, read2):
    """
    XXXX

    Parameters
    ----------
    read1 : XX
        XXX
    read2 : XX
        XXX

    Returns
    -------
    None
        XXX
    """
    if read1 and read2 != None and read1.endswith('.gz'):
        logging.critical("The files are still gzipped. Exiting")
        sys.exit(1)
    else:
        logging.info("The files have been gunzipped ...")
        with open(read1) as readFile:
            read1 = readFile.read()
        with open(read2) as readFile:
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
        logging.info("Great, a best species match was found with mash distance \
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

    ## use column (axis = 1), to create minimal dataframe
    dfSortedDropped = dfSorted.drop(['Ref ID', 'Query ID', 'KmersCount',
    'sketchSize' ], axis=1)

    ## noResult function - confirm mash distance is < than user specified
    ## even if mash distance !< user specified, return the top five hits
    noMash = noResult(dfSortedDropped, inMaxDis, bestGenusSort, bestSpeciesSort)
    bestGenus = noMash[0]
    bestSpecies = noMash[1]

    ## change order
    dfSortedDropped = dfSortedDropped[['Genus', 'Species', 'GeneBank Identifier',
    'Mash Dist', '% Seq Sim', 'P-value', 'Kmer']]
    dfTop = dfSortedDropped[:5]

    dfTop.reset_index(drop=True, inplace=True) #make index start at 0
    return bestGenus, bestSpecies, dfTop

def makeTable(dateTime, name, read1, read2, inMaxDist, results, mFlag):
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
        f.write("Input query file 1: " + read1 + "\n")
        f.write("Input query file 2: " + read2 + "\n")
        f.write("Genome size estimate for fastq files with using the -m flag: " + scData[0] + " " +"(bp)" +"\n") #make into variable
        f.write("Genome coverage estimate for fastq files with using the -m flag: " + scData[1]  + "\n") #make into variables
        f.write("Maximum mash distance (-d): " + str(inMaxDis) + "\n")
        f.write("Minimum K-mer copy number (-m) to be included in the sketch: " + str(mFlag[0]) + "\n" )
        f.write("K-mer size used for sketching: " + k_size + "\n" )
        f.write("Mash Database name: " + inMash + "\n" + "\n")
        f.write("Best species match: " + results[0] + " " + results[1] + "\n" + "\n")
        f.write("Top 5 results:" + "\n")
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
read1 = args.read1
read2 = args.read2

now = datetime.now()
dtString = now.strftime("%B %d, %Y %H:%M:%S")
dateString = now.strftime("%Y-%m-%d")

#unique name for log file based on read name
name = fastq_name(read1)
log = name + "_run"  + ".log"

req_programs=['mash', 'python']

make_output_log(log)
k_size = get_k_size(inMash)
get_input(read1, read2, inMash, inMaxDis, inKmer, k_size, inThreads)

#check_mash()

logging.info("Checking if all the required input files exist...")
check_files(read1, read2, inMash)
logging.info("Input files are present...")

logging.info("Checking if all the prerequisite programs are installed...")
for program in req_programs:
    check_program(program)
logging.info("All prerequisite programs are accessible...")

logging.info("Peforming internal system checks...")
check_mash()
logging.info("Great, internal system checks passed...")

logging.info("Begin concatenation of the fastq files...")
cat_files(read1, read2)
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
makeTable(dtString, name, read1, read2, inMaxDis, results, mFlag)
logging.info("Completed analysis for the sample: %s..." % name )
logging.info("Exiting.")
logging.info(" ")
