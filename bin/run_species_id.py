#!/usr/bin/env python3.7

## Requires python 3.7 due to singularity container
import argparse
import os
import sys
import shutil
import subprocess
import logging
import re
from io import StringIO
from datetime import datetime
from typing import Tuple, Optional, List

import pandas as pd
from tabulate import tabulate

#############################
## Argument Error Messages ##
#############################

class ParserWithErrors(argparse.ArgumentParser):
    """ My own error messages """

    def error(self, message: str) -> None:
        """Override the error method to print the message and show help."""
        print(f'\n{message}\n')
        self.print_help()
        sys.exit(2)

    def is_valid_mash(self, parser: argparse.ArgumentParser, arg: str) -> str:
        """Validate that the argument is a .msh file."""
        _, ext = os.path.splitext(arg)
        if ext != '.msh':
            parser.error(f'ERROR: This is not a file ending with .msh. Did you generate the mash sketch and specify that file to be uploaded?')
        return arg

    def is_valid_fastq(self, parser: argparse.ArgumentParser, arg: str) -> str:
        """Validate that the argument is a .fastq or .fastq.gz file."""
        _, ext = os.path.splitext(arg)
        if ext not in ('.gz', '.fastq') and not arg.endswith('.fastq.gz'):
            parser.error(f'ERROR: This is not a file ending with either .fastq or .fastq.gz. This flag requires the input of a fastq file.')
        return arg

    def is_valid_distance(self, parser: argparse.ArgumentParser, arg: str) -> str:
        """Validate that the argument is a positive float."""
        try:
            value = float(arg)
            if value < 0:
                raise ValueError
        except ValueError:
            parser.error(f'ERROR: {arg} is not a positive float, aka a number with a decimal point.')
        return arg

    def is_valid_int(self, parser: argparse.ArgumentParser, arg: str) -> str:
        """Validate that the argument is a positive integer."""
        if not arg.isdigit() or int(arg) <= 0:
            parser.error(f'ERROR: You input {arg}. This is NOT a positive integer.')
        return arg

#########################
## ArgParser Arguments ##
#########################

def argparser():
    """
    Returns argument parser for the script with messages for how to use tool.
    """

    description = (
        "A script to run and parse the output from Mash into a table listing" 
        " the top five matches from the user specified pre-built Mash Database.")

    ## use class to parse the arguments with formatted error message
    parser = ParserWithErrors(description = description)

    ## Define required and optional groups; uses lambda (anonymous) function
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
                        help="User-specified Mash distance (default: 0.05)",
                        type=lambda x: parser.is_valid_distance(parser, x))
    
    optional.add_argument("--kmer_min", "-m", default=2,
                        help="Min. k-mer copies to pass noise filter  (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    
    optional.add_argument("--num_threads", "-p", default=2,
                        help="Number of computing threads to use (default: 2)",
                        type=lambda x: parser.is_valid_int(parser, x))
    return parser

###############
## FUNCTIONS ##
###############
def make_output_log(log: str) -> None:
    """
    Creates the log file which can be appended to.
    Logs operating system (OS) information where the script is being run.
    Requires the logging package and uses traditional '%' -style formating as 
    is standard with logging module. 

    Parameters
    ----------
    log : str
        Name of the log file.

    Returns
    -------
    None
        Logs an error message if unable to get system information.
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
    try:
        sys_info = os.uname()
        logging.info("System Information:")
        logging.info("   System: %s" % sys_info.sysname)
        logging.info("   Node Name: %s" % sys_info.nodename)
        logging.info("   Release: %s" % sys_info.release)
        logging.info("   Version: %s" % sys_info.version)
        logging.info("   Machine: %s\n" % sys_info.machine)
    except AttributeError:
        # Handle cases where `os.uname` is not available 
        logging.warning("System information is not available on this platform")

def extract_base_name(filename: str) -> str:
    """
    Extract the base name from a filename, stripping common read suffixes.

    Parameters
    ----------
    filename : str
        The filename from which to extract the base name.

    Returns
    -------
    str
        The base name without any read suffixes or extensions.
    """
    
    # Remove any leading directory path
    basename = os.path.basename(filename)
    
    # Remove .gz extension if present
    if basename.endswith('.gz'):
        basename = basename[:-3]
    
    # Remove the .fastq or .fq extension
    if basename.endswith('.fastq'):
        basename = basename[:-6]
    elif basename.endswith('.fq'):
        basename = basename[:-3]

    suffixes = ('_1', '_R1_001', '_R1', '_2', '_R2_001', '_R2')

    # Remove common read suffixes
    for suffix in suffixes:
        if basename.endswith(suffix):
            return basename[:-len(suffix)]
    return basename

def fastq_name(read1: str, read2: str) -> str:
    """
    Gets stripped read name for appending to output files.

    Parameters
    ----------
    read1 : str
        Filename for the first read.
    read2 : str
        Filename for the second read.

    Returns
    -------
    str
        Base name of the file with the suffixes stripped.

    Raises
    ------
    ValueError
        If the base names from read1 and read2 do not match.
    """

    name1 = extract_base_name(read1)
    name2 = extract_base_name(read2)
    
    if name1 != name2:
        raise ValueError(f"Read1 base name ({name1}) and Read2 base name ({name2}) do not match.")
    return name1

def log_required_inputs(read1: str, read2: str, mash_db: str) -> None:
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
    Returns
    -------
    None
    """

    logging.info(
    "The user-specified required parameters:\n"
    " * Read1: %s\n" 
    " * Read2: %s\n" 
    " * Mash Database: %s\n", 
    read1, read2, mash_db
)

def log_optional_inputs(max_dis: str, min_kmer: str, k_size: str, threads: str) -> None:
    """
    Logs the command line input parameters.

    Parameters
    ----------
    max_dis : str
        Maximum distance parameter.
    min_kmer : str
        Minimum k-mer count.
    k_size : str
        Size of the K-mer.
    threads : str
        Number of threads to use.

    Returns
    -------
    None
    """
    
    logging.info(
        "The default or user-specified parameters:\n"
        " * Maximum Distance: %s\n"
        " * Minimum K-mer Count: %s\n"
        " * Size of K-mer: %s\n"
        " * Number of Threads: %s\n",
        max_dis, min_kmer, k_size, threads
    )    

def check_files(read1: str, read2: str, mash_db: str) -> None:
    """
    Checks if all the input files exist; raises an exception if file not found or if file is
    a directory.

    Parameters
    ----------
    read1 : str
        Path to input file 1.
    read2 : str
        Path to input file 2.
    mash_db : str
        Path to database file.

    Raises
    ------
    FileNotFoundError
        If any file doesn't exist or is a directory.
    ValueError
        If read1 and read2 are the same file.
    """
    
    def check_file(path: Optional[str], description: str):
        if path and not os.path.isfile(path):
                raise FileNotFoundError(f"{description} doesn't exist or is not a file: {path}")
    
    check_file(mash_db, "The database file")
    check_file(read1, "Read file 1")
    check_file(read2, "Read file 2")

    # Check if read1 and read2 are the same file
    if read1 == read2:
        raise ValueError(f"Read1 ({read1}) and Read2 ({read2}) are the same file.")

##TODO - check for corrupt gzip files or empty?

def get_k_size(mash_db: str) -> str:

    # Ensure that exception handling is as specific as possible.
    # For example, in get_k_size, catching a general Exception may hide 
    # other issues. Consider handling specific exceptions where possible.
    """
    Retrieves the k-size from the Mash database information.

    Parameters
    ----------
    mash_db : str
        Path to the Mash database.

    Returns
    -------
    Optional[str]
        The k-size value extracted from the Mash info output.
        Returns None if an error occurs or the output format is unexpected
    """
    try:
        # Run the mash info command and capture its output
        result = subprocess.run(
            ['mash', 'info', mash_db],
            capture_output=True,
            text=True,
            check=True
        )

        # Extract the k-size value from the command's output
        lines = result.stdout.splitlines()
        
        # Assumes k-size is in the 3rd line and 3rd field
        if len(lines) >= 3:
            fields = lines[2].split()
            if len(fields) >= 3:
                return fields[2]
            logging.error("Unexpected format or missing k-size information.")
        logging.error("Unexpected output format from mash info command.")
    except subprocess.CalledProcessError as e:
        logging.error("Error occurred while running the command: %s", e)
    except Exception as e:
        logging.error("An unexpected error occurred: %s", e)
        return None
           
def check_program(program_name: str) -> None:
    """
    Checks if the supplied program_name exists and if it's an appropriate version of Python.

    Parameters
    ----------
    program_name : str
        Name of the program to check if it exists.

    Raises
    ------
    SystemExit
        If the program is not found or if the Python version is insufficient.
    """
    logging.info("Checking for program %s..." % program_name)

    # Check if the program exists
    path = shutil.which(program_name)
    
    if path is None:
        logging.critical("Program %s not found! Cannot continue; dependency not fulfilled. Exiting." % program_name)
        sys.exit(1)

    # If the program is Python, check the version
    if program_name == 'python' and sys.version_info < (3, 7):
        logging.critical("Requires Python version >= 3.7. Exiting.")
        sys.exit(1)

    else:
        logging.info("Great, the program %s is loaded." % program_name)

def is_file_empty(mash_db: str) -> bool:
    """
    Checks if the specified file is empty.

    Parameters
    ----------
    mash_db : str
        Path to the file to check.

    Returns
    -------
    bool
        True if the file is empty, otherwise False.
    """
    if os.path.exists(mash_db):
        if os.path.getsize(mash_db) == 0:
            logging.critical("The mash database is empty.")
            sys.exit(1)
            return True
        else:
            return False
    else:
        logging.error("The specified file does not exist: %s", mash_db)
        return False

def check_mash() -> None:
    """
    Checks that Mash is running as expected.

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
    rounded_distance = round(df_check_dist, 4)

    # Log and validate the results
    if df_check_species == expected_species and rounded_distance == expected_dist:
        logging.info(
            "Great, the test confirms Mash is running properly.\n"
            "* Expected species: %s\n"
            "* Returned species: %s\n"
            "* Expected distance: %s\n"
            "* Returned distance: %s\n",
            df_check_species,
            expected_species,
            expected_dist,
            rounded_distance)
    else:
        logging.critical(f"The unit test to confirm species and mash value did not return expected results. Exiting.")
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
        logging.critical("Error: %s", e)
        sys.exit(1)

def update_min_kmer(calculatedKmer: int, min_kmer: int = 2) -> int:
    """
    Determine the minimum kmer value. If less than 2, set to 2.
    
    Parameters
    ----------
    calculatedKmer : int
        Value calculated based on genomeCoverage/3 in cal_kmer function
    min_kmer : int, optional, default is 2
        Input K-mer value specified by the user; used if greater than 2
    
    Returns
    -------
    int
        Integer value used for min_kmer with paired-end reads
    """
    min_kmer = max(int(min_kmer), 2)
    if min_kmer > 2:
        logging.info("User-specified a value for minimum K-mer: %s", min_kmer)
    else:
        logging.info("Min. K-mer = genome coverage divided by 3. Calculated K-mer = %s", calculatedKmer)
    return max(calculatedKmer, 2)

def run_cmd(command: List[str]) -> subprocess.CompletedProcess:
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
        result = subprocess.run(
            command,
            capture_output=True,
            check=True,
            text=True
        )
        logging.info("Executed command: %s", ' '.join(command))
    except subprocess.CalledProcessError as e:
        logging.critical("The following command had an error:\n%s", e)
        sys.exit(1)
    
    return result

def cal_kmer(mash_db: str, threads: int, min_kmer: int) -> Tuple[int, str, str]:
    """
    Calculate the minimum k-mer value and genome statistics from mash dist output.
    
    Parameters
    ----------
    mash_db : str
        The path to the mash database.
    threads : int
        The number of threads to use.
    min_kmer : int
        The minimum K-mer value.
    
    Returns
    -------
    tuple
        A tuple containing:
        - mFlag : int
            The minimum k-mer copies to use (-m flag).
        - gSize : str
            Estimated genome size.
        - gCoverage : str
            Estimated genome coverage.
    """    
    fastqCmd1 = ['mash', 'dist', str(mash_db), '-r', 'myCatFile', '-p', str(threads), '-S', '42']
    outputFastq1 = run_cmd(fastqCmd1)   

    stderr_lines = outputFastq1.stderr.splitlines()
    if len(stderr_lines) < 2:
        raise ValueError("Unexpected output format from the command.")
    
    gSize, gCoverage = stderr_lines[0][23:], stderr_lines[1][23:]
    logging.info("Estimated Genome Size to determine -m flag: %s", gSize)
    logging.info("Estimated Genome Coverage to determine -m flag: %s", gCoverage)
    
    minKmers = int(float(gCoverage) / 3)
    mFlag = update_min_kmer(minKmers, min_kmer)
    
    return mFlag, gSize, gCoverage

def get_results(mFlag: int, threads: int, mash_db: str) -> subprocess.CompletedProcess:
    """
    Runs the mash distance command using value from cal_kmer and extracts genome size 
    and coverage from the command output.
    
    Parameters
    ----------
    mFlag : int
        The -m flag value for the mash command.
    threads : int
        Number of threads to use with the mash command.
    mash_db : str
        Path to the mash database.
    
    Returns
    -------
    subprocess.CompletedProcess
        The result of the mash command.
    """
    fastq_cmd2 = [
        'mash', 'dist', '-r', '-m', str(mFlag),
        str(mash_db), 'myCatFile', '-p', str(threads), '-S', '123456'
    ]
    
    output = run_cmd(fastq_cmd2)
    stderr_lines = output.stderr.splitlines()
    
    if len(stderr_lines) >= 2:
        logging.info("Estimated Genome size calculated with the -m flag:: %s", stderr_lines[0][23:])
        logging.info("Estimated Genome coverage calculated with the -m flag:: %s", stderr_lines[1][23:])
    else:
        logging.warning("Unexpected output format from mash command.")
    
    return output

def parse_results(cmd: subprocess.CompletedProcess, in_max_dis: float) -> Tuple[str, str, pd.DataFrame]:
    """
    Run the initial command and parse the results from mash.
    
    Parameters
    ----------
    cmd : subprocess.CompletedProcess
        The completed process object from running the mash command.
    in_max_dis : float
        User-specified maximum mash distance for filtering results.
    
    Returns
    -------
    tuple
        - best_genus: str
            The most likely genus of the isolate tested.
        - best_species: str
            The most likely species of the isolate tested.
        - df_top: pandas.DataFrame
            The top five results from sorting Mash output.
    """
    df = pd.read_csv(StringIO(cmd.stdout), sep='\t', names=['Ref ID', 'Query ID', 'Mash Dist', 'P-value', 'Kmer'])
    df[['Genus', 'Species']] = df['Ref ID'].str.split('_', 1, expand=True)
    df[['Species', 'GeneBank Identifier']] = df['Species'].str.split('_', 1, expand=True)
    df['GeneBank Identifier'] = 'GCA' + df['GeneBank Identifier'].str.split('GCA', expand=True).iloc[:, -1]
    df[['KmersCount', 'sketchSize']] = df['Kmer'].str.split("/", expand=True)
    df['KmersCount'] = df['KmersCount'].astype(int)
    df['% Seq Sim'] = (1 - df['Mash Dist']) * 100
    
    df_sorted = df.sort_values('KmersCount', ascending=False).drop(['Ref ID', 'Query ID', 'KmersCount', 'sketchSize'], axis=1)
    
    best_genus_sort, best_species_sort = is_tie(df)
    best_genus, best_species = no_result(df_sorted, in_max_dis, best_genus_sort, best_species_sort)
    
    df_top = df_sorted[['Genus', 'Species', 'GeneBank Identifier', 'Mash Dist', '% Seq Sim', 'P-value', 'Kmer']].head(5)
    df_top.reset_index(drop=True, inplace=True)
    
    return best_genus, best_species, df_top

def is_tie(df: pd.DataFrame) -> Tuple[str, str]:
    """
    Determine if the k-mers count value is a tie with the second top isolate.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing k-mers count, genus, and species information.
    
    Returns
    -------
    tuple of str
        - A string indicating the best genus, either a specific genus or a tie message.
        - A string indicating the best species, either a specific species or a blank string if tied.
    """
    logging.info("Checking the number of entries in the DataFrame...")
    
    # Check if there is only one entry
    if len(df) == 1:
        best = df.iloc[0]
        logging.info("Only one entry found. Returning that entry.")
        return best['Genus'], best['Species']
    
    df_sort = df.sort_values('KmersCount', ascending=False)
    logging.info("Checking if k-mers count is tied for top 2 results...")
    
    if df_sort.iloc[0]['KmersCount'] == df_sort.iloc[1]['KmersCount']:
        logging.info("The top two isolates have the same number of matching k-mers, indicating a tie.")
        return "This was a tie, see the top 5 results below", ""
    
    best = df_sort.iloc[0]
    logging.info("There was no tie for k-mers count in the top two species.")
    return best['Genus'], best['Species']

def no_result(in_file: pd.DataFrame, in_max_dis: float, best_g: str, best_s: str) -> Tuple[str, str]:
    """
    Determine if the top hit mash distances are >= user-specified mash distance.
    
    Parameters
    ----------
    in_file : pandas.DataFrame
        DataFrame containing the mash results with a column 'Mash Dist' representing the mash distances.
    in_max_dis : float
        User-specified maximum mash distance as a cut-off.
    best_g : str
        Current best species message.
    best_s : str
        Current best species placeholder.
    
    Returns
    -------
     Tuple[str, str]
        A tuple containing:
        - The updated best species message.
        - The updated best species placeholder.
    """
    logging.info("Confirming that best match is less than user-specified distance...")
    
    if in_file['Mash Dist'].values[0] < float(in_max_dis):
        logging.info("A best species match was found with mash distance less than %s", in_max_dis)
    else:
        best_g = f"No matches found with mash distances < {in_max_dis}..."
        best_s = " "
        logging.info("No matches found with mash distances < %s", in_max_dis)
    
    return best_g, best_s

def make_table(date_time: str, name: str, read1: str, read2: str, max_dist: float, results: Tuple[str, str, pd.DataFrame], m_flag: Tuple[int, int]) -> None:
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
    m_flag : tuple
        Contains the minimum k-mer copy number and k-mer size.
    
    Returns
    -------
    None
        Writes the results to a text file.
    """
    file_name = f"{name}_results_{date_time}.txt"
    
    with open(file_name, 'a+') as f:
        f.write(f"\nLegionella Species ID Tool using Mash\n")
        f.write(f"Date and Time = {date_time}\n")
        f.write(f"Input query file 1: {read1}\n")
        f.write(f"Input query file 2: {read2}\n")
        f.write(f"Maximum Mash distance (-d): {max_dist}\n")
        f.write(f"Minimum K-mer copy number (-m) to be included in the sketch: {m_flag[0]}\n")
        f.write(f"K-mer size used for sketching: {m_flag[1]}\n")
        f.write(f"Mash Database name: {mash_db}\n\n")
        f.write(f"Best species match: {results[0]} {results[1]}\n\n")
        f.write("Top 5 results:\n")
        f.write(u'\u2500' * 100 + "\n")
        f.write(tabulate(results[2], headers='keys', tablefmt='pqsl', numalign="center",
                         stralign="center", floatfmt=(None, None, None, ".5f", ".3f", ".8e"),
                         showindex=False) + "\n")

if __name__ == '__main__':
    # Argument parsing
    parser = argparser()
    args = parser.parse_args()

    mash_db = args.database
    max_dis = args.max_dist
    min_kmer = args.kmer_min
    threads = args.num_threads
    read1 = args.read1
    read2 = args.read2

    now = datetime.now()
    date_time = now.strftime("%Y-%m-%d")

    name = fastq_name(read1, read2)
    log = f"{name}_run.log"

    req_programs = ['mash', 'python']

    make_output_log(log)
    log_required_inputs(read1, read2, mash_db)
    log_optional_inputs(max_dis, min_kmer, get_k_size(mash_db), threads)
    logging.info("Base name for the files: %s", name)

    check_files(read1, read2, mash_db)
    logging.info("Input files are present...")

    for program in req_programs:
        check_program(program)
    logging.info("All prerequisite programs are accessible...\n")

    check_mash()
    logging.info("Internal system checks passed...")

    is_file_empty(mash_db)
    logging.info("Mash database is not empty...")

    cat_files(read1, read2)
    logging.info("Files concatenated successfully...")

    mFlag = cal_kmer(mash_db, threads, min_kmer)
    logging.info("Minimum copies of each K-mer identified...")

    outputFastq2 = get_results(mFlag[0], threads, mash_db)
    logging.info("Mash dist command completed...")

    results = parse_results(outputFastq2, max_dis)
    logging.info("Results parsed successfully...")

    make_table(date_time, name, read1, read2, max_dis, results, mFlag)
    logging.info("Analysis completed for sample: %s", name)
    logging.info("EXITING!")
