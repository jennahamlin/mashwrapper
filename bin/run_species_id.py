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

##TO DO - check for corrupt gzip files
##TO DO - check if the beginning of the file name is a match between the two files

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

def cat_files(read1: str, read2: str) -> None:
    """
    Concatenates the contents of two files and writes the result to 'myCatFile'.

    Parameters
    ----------
    read1 : str
        Path to the first file to be concatenated.
    read2 : str
        Path to the second file to be concatenated.

    Returns
    -------
    None
        Exits the program if the files are gzipped.
    """
    # Check if files are gzipped
    if (read1.endswith('.gz') or read2.endswith('.gz')):
        logging.critical("One or both files are gzipped. Exiting.")
        sys.exit(1)

    logging.info("The files are not gzipped. Proceeding with concatenation...")

    # Concatenate the contents of the files
    try:
        with open(read1, 'r') as file1, open(read2, 'r') as file2:
            combined_content = file1.read() + file2.read()

        # Write the concatenated content to 'myCatFile'
        with open('myCatFile', 'w') as output_file:
            output_file.write(combined_content)

        logging.info("Files have been successfully concatenated and written to 'myCatFile'.")

    except FileNotFoundError as e:
        logging.critical(f"Error: {e}")
        sys.exit(1)
    except IOError as e:
        logging.critical(f"IO Error: {e}")
        sys.exit(1)

def minKmer(calculatedKmer, inKmer=2):
    """
    Determine the minimum kmer value. If less than 2, set to 2.

    Parameters
    ----------
    calculatedKmer : int
        Value that is calculated based on genomeCoverage/3
    inKmer : int, optional, default is 2
        Input kmer value specified by the user; used instead of calculatedKmer if greater than 2

    Returns
    ----------
    int
        Integer value used for min_kmer with paired-end reads
    """
    inKmer = int(inKmer)  # Ensure inKmer is an integer

    # Check if user-specified kmer is greater than 2
    if inKmer > 2:
        logging.info(f"User specified a value for minimum kmer: {inKmer} ...")
        return inKmer

    # Log information about default behavior
    logging.info("Should kmer value be different than default (2)...")
    logging.info("Min. kmer = genome coverage divided by 3...")

    # Determine minimum kmer value
    if calculatedKmer < 2:
        logging.info("The calculated kmer is less than 2, so will use 2...")
        return 2
    return calculatedKmer

def run_cmd(command):
    """
    Executes a shell command and logs its output. Exits the program on error.

    Parameters
    ----------
    command : list of str
        The command to be executed, provided as a list of arguments.

    Returns
    -------
    subprocess.CompletedProcess
        The result of the executed command, including stdout, stderr, and return code.
    """

    try:
        # Execute the command
        result = subprocess.run(
            command, 
            capture_output=True,
            check=True,
            text=True
        )
        # Log the command executed
        new_cmd = ' '.join(command)
        logging.info(f"This is the command...\n{new_cmd}")
    except subprocess.CalledProcessError as e:
        # Log critical error and exit if the command fails
        logging.critical(f"CRITICAL ERROR. The following command had an error:\n{e}")
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
    fastqCmd1 = ['mash', 'dist', str(inMash), '-r', 'myCatFile', '-p', str(inThreads), '-S', '42']

    outputFastq1 = run_cmd(fastqCmd1)

    ## get genome size and coverage; will provide as ouput for user
    ## currently gets message at postion 0 but this is some error aboutr lang locale
    ## so it errors out if I can change the two split lines to be on 3 and 4 
    ## than the value can be calcuated
    gSize = outputFastq1.stderr.splitlines()[3]
    gSize = gSize[23:]
    logging.info("Estimated Genome Size to determine -m flag: %s " % gSize)
    gCoverage = outputFastq1.stderr.splitlines()[4]
    gCoverage = gCoverage[23:]
    logging.info("Estimated Genome coverage to determine -m flag: %s " % gCoverage)

    minKmers = int(float(gCoverage))/3
    minKmers = int(float(minKmers))

    ## this is used the calucate the minimum kmer copies to use (-m flag)
    mFlag = minKmer(minKmers, inKmer) # returned as an integer
    return mFlag, gSize, gCoverage

def get_results(mFlag, inThreads):
    fastqCmd2 = ['mash', 'dist', '-r', '-m', str(mFlag), str(inMash), 'myCatFile', '-p', str(inThreads), '-S', '123456']
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

def is_tie(df):
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

def no_result(in_file, in_max_dis, best_g, best_s):
    """
    Determine if the top hit mash distances are >= user-specified mash distance.

    Parameters
    ----------
    in_file : pandas.DataFrame
        Parsed results from running mash.
    in_max_dis : float
        User-specified maximum mash distance as a cut-off.
    best_g : str
        Current best species message.
    best_s : str
        Current best species placeholder.

    Returns
    ----------
    tuple
        Updated best species message and placeholder.
    """
    in_max_dis = float(in_max_dis)
    logging.info("Confirming that best match is less than user-specified distance...")

    if in_file['Mash Dist'].values[0] < in_max_dis:
        logging.info(f"Great, a best species match was found with mash distance less than {in_max_dis}...")
    else:
        best_g = f"No matches found with mash distances < {in_max_dis}..."
        best_s = " "
        logging.info(f"No matches found with mash distances < {in_max_dis}...")

    return best_g, best_s

def parse_results(cmd, in_max_dis):
    """
    Run the initial command and parse the results from mash.

    Parameters
    ----------
    cmd : subprocess.CompletedProcess
        The completed process object from running the mash command.
    in_max_dis : float
        User-specified maximum mash distance for filtering results.

    Returns
    ----------
    tuple
        - best_genus: str
            The most likely genus of the isolate tested.
        - best_species: str
            The most likely species of the isolate tested.
        - df_top: pandas.DataFrame
            The top five results from sorting Mash output.
    """
    # Read the mash output into a DataFrame
    df = pd.read_csv(StringIO(cmd.stdout), sep='\t',
                     names=['Ref ID', 'Query ID', 'Mash Dist', 'P-value', 'Kmer'],
                     index_col=False)

    # Extract Genus, Species, and GeneBank Identifier
    df[['Genus', 'Species']] = df['Ref ID'].str.split('_', 1, expand=True)
    df[['Species', 'GeneBank Identifier']] = df['Species'].str.split('_', 1, expand=True)
    df['GeneBank Identifier'] = 'GCA' + df['GeneBank Identifier'].str.split('GCA', expand=True).iloc[:, -1]

    # Split Kmers and convert to integer
    df[['KmersCount', 'sketchSize']] = df['Kmer'].str.split("/", expand=True)
    df['KmersCount'] = df['KmersCount'].astype(int)

    # Calculate % sequence similarity
    df['% Seq Sim'] = (1 - df['Mash Dist']) * 100

    # Sort by KmersCount and handle ties
    df_sorted = df.sort_values('KmersCount', ascending=False)
    df_sorted_dropped = df_sorted.drop(['Ref ID', 'Query ID', 'KmersCount', 'sketchSize'], axis=1)

    # Determine the best genus and species
    best_genus_sort, best_species_sort = is_tie(df_sorted)
    best_genus, best_species = no_result(df_sorted_dropped, in_max_dis, best_genus_sort, best_species_sort)

    # Prepare the top five results
    df_top = df_sorted_dropped[['Genus', 'Species', 'GeneBank Identifier', 'Mash Dist', '% Seq Sim', 'P-value', 'Kmer']].head(5)
    df_top.reset_index(drop=True, inplace=True)  # Reset index to start at 0

    return best_genus, best_species, df_top

def make_table(date_time, name, read1, read2, max_dist, results, m_flag):
    """
    Parse results into a text output file including relevant variables.

    Parameters
    ----------
    date_time : str
        Current date and time for when analysis is run.
    name : str
        Base name for the output file.
    read1 : str
        Path to the first query file.
    read2 : str
        Path to the second query file.
    max_dist : float
        User-specified maximum mash distance.
    results : tuple
        Output from running and parsing mash commands, where:
        - results[0] is the best genus.
        - results[1] is the best species.
        - results[2] is a pandas DataFrame of the top results.
    m_flag : list
        Contains the minimum k-mer copy number and k-mer size.

    Returns
    ----------
    None
        Writes the results to a text file.
    """
    # Define file name
    file_name = f"{name}_results_{date_time}.txt"

    # Open file in append mode
    with open(file_name, 'a+') as f:
        # Write headers and variable values
        f.write(f"\nLegionella Species ID Tool using Mash\n")
        f.write(f"Date and Time = {date_time}\n")
        f.write(f"Input query file 1: {read1}\n")
        f.write(f"Input query file 2: {read2}\n")
        f.write(f"Maximum mash distance (-d): {max_dist}\n")
        f.write(f"Minimum K-mer copy number (-m) to be included in the sketch: {m_flag[0]}\n")
        f.write(f"K-mer size used for sketching: {m_flag[1]}\n")
        f.write(f"Mash Database name: {m_flag[2]}\n\n")
        f.write(f"Best species match: {results[0]} {results[1]}\n\n")
        
        # Write top results with table formatting
        f.write("Top 5 results:\n")
        f.write(u'\u2500' * 100 + "\n")
        f.write(tabulate(results[2], headers='keys', tablefmt='pqsl', numalign="center",
                         stralign="center", floatfmt=(None, None, None, ".5f", ".3f", ".8e"),
                         showindex=False) + "\n")


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
results = parse_results(outputFastq2, inMaxDis)
logging.info("Okay, completed parsing of the results...")

logging.info("Generating table of results as a text file...")
make_table(dtString, name, read1, read2, inMaxDis, results, mFlag)
logging.info("Completed analysis for the sample: %s..." % name )
logging.info("Exiting.")
logging.info(" ")