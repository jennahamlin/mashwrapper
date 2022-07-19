#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import argparse
import os
import re
import sys

def parse_args(args=None):
    Description = "Check contents of organism list. Assumes either genus species\
 or species per row."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)

def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip())
    print(error_str)
    sys.exit(1)

def is_file_empty(file_path):
    """ Check if file is not empty by confirming if its size is >0 bytes"""
    # Check if file exist and it is empty
    if os.path.exists(file_path) and os.stat(file_path).st_size != 0:
        print("Okay file exists and is not empty. Continuing...")
        return file_path
    else:
        print_error("No entries to process! Exiting.", \
 "Samplesheet: {}".format(file_path))

def what_is_format(file_ending):
    if file_ending.endswith((".csv", ".tsv", ".msh", ".mash", ".db")):
        print_error("This is not the correct file type, should be .txt! Exiting.",\
 "Samplesheet: {}".format(file_ending))
    else:
        print("Okay, this is a text file. Continuing...")
        return(file_ending)

#filename: str, n: int
def head_file(filename, n):
    try:
        with open(filename) as f:
            head_lines = [next(f).rstrip() for x in range(n)]
    except StopIteration:
        with open(filename) as f:
            head_lines = f.read().splitlines()
    return head_lines

#filename: str, n=2
def detect_delimiter(filename, n=2):
    sample_lines = head_file(filename, n)
    common_delimiters= [',',';','\t',' ','|',':',' ']
    for d in common_delimiters:
        ref = sample_lines[0].count(d)
        if ref > 0:
            if all([ ref == sample_lines[i].count(d) for i in range(1,n)]):
                return d
    return d

def is_space(d):
    delim = ' '
    if d == delim:
        print("okay this is a space")
    else:
        print_error("The delimiter is something other than a space! Exiting.", \
 "Samplesheet: {}".format(d))


def check_organismsheet(file_in, file_out):
    """
    This function checks that the organismsheet follows the following structure:
    Note - there is no header required for this file

    legionella
    fluoribacter
    vibrio natriegens

    """

    file_exists = is_file_empty(file_in)
    correct_ending = what_is_format(file_in)
    d = detect_delimiter(file_in, 1)
    is_space(d)

    with open(file_in, "r") as input:

    # Creating "gfg output file.txt" as output
    # file in write mode
        with open(file_out, "w") as output:

        # Writing each line from input file to
        # output file using loop
            for line in input:
                output.write(line)

def main(args=None):
    args = parse_args(args)
    check_organismsheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
