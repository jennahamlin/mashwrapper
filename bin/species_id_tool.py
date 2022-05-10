#!/usr/bin/env python3

import os
import pathlib
import os.path
import sys
import subprocess
import argparse
import re
import pandas as pd
from io import StringIO
from datetime import datetime
from tabulate import tabulate

## Author: Jenna Hamlin (jhamlin at cdc.gov, ptx4 at cdc.gov)

## Application Name: Legionella Species Identification Tool

## Description: Wrapper around the tool Mash, which uses a pre-defined database
## of all Legionella species available from NCBI

## Input: fastq or fasta files

## Output: Text file of isolate tested & top 5 closest matches from mash sketch
## if Mash distance is less than user specified distance via -m flag

## Requirements: Python 3, Mash Program, Mash Sketch of organism of interest
## I have generated a pre-built mash sketch of all Legionella species available
## as of Oct. 2021. Ideally, this database should be updated yearly.

## get date and time and format as (Month Day, Year Hours:Minutes:Seconds)
## used in reuslt text file
today = datetime.now()
dateString = today.strftime("%Y-%m-%d")

## dynamically get the folder path where code is run from and assumes
## this is where the reads are located. Will this work on the portal?
folderPath = pathlib.Path().resolve()

## if necessary can join to a path output via posixpath.join
#import posixpath
#print(folderPath)ls
#folderPath2 = posixpath.join(folderPath, "tmp")
#print(folderPath2)

########################
## ArgumentParser Checks
########################
def mashValid(mash):
    """
    check that file extensions match expection of a mash file; if not exit

    Parameters
    ----------
    mash : str
        Pre-assemblied mash sketch file

    Raises
    ----------
    ArgumentTypeError
        If file is not included or file does not end with .msh; then exit
    """

    base, ext = os.path.splitext(mash)
    if ext not in ('.msh'):
        raise argparse.ArgumentTypeError('This is not a file ending with .msh\
 Did you generate the mash database correctly?')
    return mash

def distValid(distIn):
    """
    check for a positive integer; if not exit

    Parameters
    ----------
    distIn : float
        Value to test against calcuated Mash distance, should be pretty small

    Raises
    ----------
    ArgumentTypeError
        If not positive value; then exit
    """

    fdistIn = float(distIn)
    if fdistIn <= 0:
        raise argparse.ArgumentTypeError('%s is not a positive number (e.g., a \
 float aka a number with a decimal point)' % distIn)
    return fdistIn

def valValid(valIn):
    """
    check for a positive integer; if not exit

    Parameters
    ----------
    valIn : int
        Test to see if value entered is positive

    Raises
    ----------
    ArgumentTypeError
        If not positive value; then exit
    """

    ivalIn = int(valIn)
    if ivalIn <= 0:
        raise argparse.ArgumentTypeError('%s is not a positive integer' % valIn)
    return ivalIn

########################
## ArgParser Arguments
########################

## NOTE do I need an argument for sketch size?

parser = argparse.ArgumentParser(description='This is a wrapper around the\
 program Mash and uses a pre-defined database. The tool was built to identify \
 a Legionella isolate of interest to the species level. Thus our pre-defined \
 database is all whole genome isolates of Legionella from NCBI. The required\
 input is a single fasta file (.fasta, .fna, .fas, .fa) or paired-end fastq\
  files (*_R1.fastq and *_R2.fastq), along with your mash database (.msh).',\
 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-mash_sketch', type = mashValid, required = True,\
 help = 'Path to user defined mash database (.msh)') ## should add default to my database
parser.add_argument('-max_dist', type = distValid, default = 0.05,\
 help = 'Maximium Mash distance to return as a species match')
parser.add_argument('-min_kmer', type = valValid,\
 help = 'Minimum copies of each kmer count to use')
parser.add_argument('-nthreads', type = valValid, default = 4,\
 help = 'Number of threads to use')

args = parser.parse_args()

## assign args to variable name, to use in functions and for printing
inMash = args.mash_sketch
inMaxDist = args.max_dist
inKmer = args.min_kmer
inThreads = args.nthreads

########################
## Helper Functions
########################
def whichFile(files1, files2=None):
    """
    Determine if two files or one single file uploaded

    Parameters
    ----------
    files1 : str
        This is either the first fastq file or single fasta files
    files2 : str
        This is the second fastq file

    Returns
    ----------
    str
        The name of file(s) tested with path removed
    """

    if files2 is not None:
        files1 = os.path.basename(files1)
        files2 = os.path.basename(files2)
        f1f2 = [files1, files2]
        combinedSeq = (' and '.join(f1f2))
        return combinedSeq
    else:
        files1 = os.path.basename(files1)
        return files1


def minKmer(calculatedKmer, min_kmer=None):
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

    if args.min_kmer:
        print("Okay, user specified this value for minimum kmer: ", args.min_kmer)
        return min_kmer
    elif (calculatedKmer < 2):
        print("The calucated kmer is less than 2, so will use 2")
        calculatedKmer = 2
        return calculatedKmer
    else:
        print("This is the calcuated kmer: ", calculatedKmer)
        return calculatedKmer


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

## assumes based on position the first and second value are always
## column 11; row 1 and 2
    if dfSort.iloc[0:1, 10:11].equals(dfSort.iloc[1:2, 10:11]):
        bestGenus = "This was a tie, see the top 5 results below"
        bestSpecies = " "
        return bestGenus, bestSpecies
    else:
        best = dfSort.head(1)
        bestGenus = best['Genus']
        bestGenus = bestGenus.str.cat(sep='\n')
        bestSpecies = best['Species']
        bestSpecies = bestSpecies.str.cat(sep='\n')
        return bestGenus, bestSpecies


def noResult(inFile, inDist, bestGenus, bestSpecies):
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

    if (inFile['Mash Dist'].values[0] < inMaxDist):
        return bestGenus, bestSpecies
        print("Okay, a best species match was found with mash distance less \
 than {}, including the best species in the results...". format(inMaxDist))

    else:
        bestGenus = "No matches found with mash distances < {}".format(inMaxDist)
        bestSpecies = " "
        return bestGenus, bestSpecies
        print("No matches found with mash distances < {}".format(inMaxDist))


def generateResults(cmd, genomeSize, genomeCoverage, mFlag):
    """
    run initial command and parse the results from mash

    Parameters
    ----------
    cmd : list
        Initial command to run for either fasta or fastq
    genomeSize : int, optional
        Estimated genome size from testing pair-end fastq files
    genomeCoverage : int, optional
        Estimated genome coverage from testing pair-end fastq files

    Returns
    ----------
    bestGenus
        The most likely genus of the isolate tested
    bestSpecies
        The most likely species of the isolate tested
    dfTop
        The top five results from sorting Mash output
    genomeSize
        Estimated genome size if fastq files were input
    genomeCoverage
        Estimated coverage of genome if fastq files were input; assumes that
        it is calculated like so: Coverage = N x L/G;
        N = number of reads, L = average read legnth, G = genome size
    """
## cmd is created in next function
    cmdRun = subprocess.run(cmd, capture_output=True, check=True, text=True)

    df = pd.read_csv(StringIO(cmdRun.stdout), sep='\t', #convert with StringIO and added headers (for development)
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
    #print("Line 309")
    #print(df)
    #print(type(df))
    #print(type(df['P-value']))

## split the kmers for sorting because xx/xxxx
    df[['KmersCount','sketchSize']] = df.Kmer.str.split("/", expand=True,)
    df['KmersCount'] = df.KmersCount.astype(int)

## add column that is (1 - Mash Distance) * 100, which is % sequence similarity
    df['% Seq Sim'] =  (1 - df['Mash Dist'])*100

## now sort and get top species; test for a tie in kmerscount value
    dfSorted = df.sort_values('KmersCount', ascending=False)
    dfSortOut = isTie(dfSorted)
    bestGenus = dfSortOut[0]
    bestSpecies = dfSortOut[1]

## use column (axis = 1), to create minimal dataframe
    dfSortedDropped = dfSorted.drop(['Ref ID', 'Query ID', 'KmersCount',
    'sketchSize' ], axis=1)

## noResult function - confirm mash distance is < than user specified
## even if mash distance !< user specified, return the top five hits
    noMash = noResult(dfSortedDropped, inMaxDist, bestGenus, bestSpecies)
    bestGenus = noMash[0]
    bestSpecies = noMash[1]

## change order
    dfSortedDropped = dfSortedDropped[['Genus', 'Species', 'GeneBank Identifier',
    'Mash Dist', '% Seq Sim', 'P-value', 'Kmer']]
    dfTop = dfSortedDropped[:5]

##TO DO - scienfitic notation for P-value

    dfTop.reset_index(drop=True, inplace=True) #make index start at 0
    return bestGenus, bestSpecies, dfTop, genomeSize, genomeCoverage, mFlag

########################
## Main function
########################

def wrapperScript(files1, files2=None):
    """
    wrapper script around the mash commands and parsing results for output

    Parameters
    ----------
    files1 : str
        string of file path and name to file 1 (either fastq1 or fasta)

    files2 : str, optional
        string of file path and name to optional file 2 (fastq2)

    Returns
    ----------
    function
        wrapper function to generate & parse results for isolates within folder

    """
    print("Internal Check: this is either a single fasta file or one of the\
 fastq paired-end reads: ", files1)
    print("Internal Check: this is the second paired-end read that will be\
 used for fastq analysis: ", files2)

 ## This script assumes the inclusion of the mash database. As such one cannot
 ## specify the -k (K-mer size) flag because it is inherited from the sketch

    def parseResults(files1, files2=None):
        """
        function to run mash and return results & varibles of interest (e.g.,
        genome size, genome coverage, genus and species of best match)

        Parameters
        ----------
        files1 : str
            string of file path and name to file 1 (either fastq1 or fasta)

        files2 : str, optional
            string of file path and name to optional file 2 (fastq2)

        Returns
        ----------
        function
            determine file type, run script to generate and parse results

        """
        print("Now, testing against mash database...")
        files1Name = os.path.basename(files1)

    ## if len(files1) == 2: #when fastqs this is a tuple, so length is 2
        if files2 is not None:
            files2Name = os.path.basename(files2)
            with open(files1) as readFile:      ##will I be able to do this via geneflow??
                read1 = readFile.read()
            with open(files2) as readFile:
                read2 = readFile.read()
            read1 += read2
            with open('myCatFile', 'w') as readFile:
                readFile.write(read1)
            readFile.close()

            fastqCmd = ['mash', 'dist', '-r', inMash, 'myCatFile']

            outputFastq = subprocess.run(fastqCmd, capture_output=True,
            check=True, text=True)

        ## get genome size and coverage; will provide as ouput for user
            gSize = outputFastq.stderr.splitlines()[0]
            gSize = gSize[23:]
            print("Genome Size: ", gSize)
            gCoverage = outputFastq.stderr.splitlines()[1]
            gCoverage = gCoverage[23:]
            print("Genome coverage: ", gCoverage)

            minKmers = int(float(gCoverage))/3
            minKmers = int(float(minKmers))

        ## this is used the calucate the minimum kmer copies to use (-m flag)
            mFlag = minKmer(minKmers, inKmer) # returned as an integer

            fastqCmd2 = ['mash', 'dist', '-r', '-m',  str(mFlag), inMash,
            'myCatFile']

            print("This is the fastq command: ", fastqCmd2)


            fastqOut = generateResults(fastqCmd2, gSize, gCoverage, str(mFlag))

            os.remove('myCatFile')

            return fastqOut

        elif files1Name: #when single fasta this is a string and the length is greater than 2
            fastaCmd = ['mash', 'dist', inMash, files1Name]
            print("This is the fasta command: ", fastaCmd)

            gSize = 'N/A as the input was fasta'
            gCoverage = 'N/A - as the input was fasta'
            mFlag = 'N/A - as the input was fasta'
            fastaOut = generateResults(fastaCmd, gSize, gCoverage, mFlag)

            return fastaOut

## this is the input variable for when writing the results
    filesTested = whichFile(files1, files2)

## this is the data parsed for makeTable function
    parsedResults = parseResults(files1, files2)

    print("Line 456")
    print(parsedResults)
    print("Writing the results to text file...")

    def makeTable(dateTime, maxDist, results):
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
            f.writelines("\n" + "Legionella Species ID Tool using Mash" + "\n")
            f.writelines("Date and Time = " + dtString + "\n") #+str(variable)
            f.write("Input query file(s):" + filesTested + "\n")
            f.write("Maximum mash distance: " + str(maxDist) + "\n")
            f.write("Genome size estimate for fastq files: " + results[3] + " " +"(bp)" +"\n") #make into variable
            f.write("Genome coverage estimate for fastq files: " + results[4]  + "\n") #make into variables
            f.write("Minimum kmer copy number to be included in the sketch: " + results[5] + "\n" + "\n")
            f.write("Best species match: " + results[0] + " " + results[1] + "\n" + "\n")
            f.write("Top 5 hits:" + "\n")
            f.writelines(u'\u2500' * 105 + "\n")
            f.writelines(tabulate(results[2], headers='keys', tablefmt='pqsl', numalign="center", stralign="center") + "\n")
            f.writelines(u'\u2500' * 105 + "\n")

    makeTable(dtString, inMaxDist, parsedResults)

    print("The analysis for this isolate is complete, checking other files." + "\n")

filesTested = []

"""
This loop checks all file endings in the current folder,
confirms that both fastq files are uploaded and runs the wrapperScript funciton
for any file(s) that pass the checks
"""
print("\n" + 'Starting the tool....' + "\n")

for file in os.listdir(folderPath): #folderPath is defined at beginning
    print('Checking file: ', file)

    f = os.path.join(folderPath, file)
    base, ext = os.path.splitext(f)

## double parathensis makes the file extensions one tuple and thus one argument
## as .endswith only allows 3 arguments

    if file.endswith(("fastq", "fa", "fas", "fna", "fasta")):   #what about .fastq.gz
        print("Great, file matches one of the possible extensions...")
        now = datetime.now()
        dtString = now.strftime("%B %d, %Y %H:%M:%S")

    ## use regular expression to get the file name with _R*.001.fastq removed
    ## add that file name to fileTest list for checking that all files have been
    ## tested
        oneSame = re.split('_R..001.fastq', file)

        if ext in('.fastq') and oneSame[0] in filesTested:
            print("This is the other read from fastq paired-end reads. No need\
 to test as it was already completed. Checking other files..." + "\n")
            continue

        elif ext in('.fastq') and oneSame[0] not in filesTested:

            f1 = f
            f1Check = os.path.exists(f1) ## should return true, if exists

            if f1.endswith("_R1_001.fastq"):
                f2 = f1.replace('_R1_001.fastq', '_R2_001.fastq')
                f2Check = os.path.exists(f2)
            elif f1.endswith("_R2_001.fastq"):
                f2 = f1.replace('_R2_001.fastq', '_R1_001.fastq')
                f2Check = os.path.exists(f2)
            filesTested.append(oneSame[0])

            if f1Check and f2Check:
                print('Okay, both fastq files are available, continuing....')

                try:
                    wrapperScript(f1, f2)
                except:
                    print("Mash could not execute. Is it loaded via module load?" + "\n")

        elif ext not in ('.fastq'):
            print("Looks like a single FASTA file, will continue....")
            fastaIn = f
            try:
                wrapperScript(fastaIn)
            except:
                print("Mash could not execute. Is it loaded via module load?"+ "\n")

    else:
        print('This file does not match one of the required file endings:\
 .fastq, .fasta, .fas, .fna, .fa. The tool will skip this file.' + "\n")

print("All files have been tested, exiting the program.")

## TO DO - scienfitic notation is rounding
## TO DO - what about 'loading' mash tool; see version 1 with a possible solution
## TO DO - any other checks and error messages to add...
