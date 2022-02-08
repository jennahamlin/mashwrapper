#!/usr/bin/env python

import argparse, sys, os

class ParserWithErrors(argparse.ArgumentParser):
    def error(self, message):
        print('{0}\n\n'.format(message))
        self.print_help()
        sys.exit(2)

    def is_valid_mash(self, parser, arg):
        base, ext = os.path.splitext(arg)
        if ext not in ('.msh'):
            parser.error('This is not a file ending with .msh.\
 Did you generate the mash sketch and specified that file to be uploaded?')
        else:
            return arg

    def is_valid_fastq(self, parser, arg):                                      #should I
        base, ext = os.path.splitext(arg)
        if ext not in ('.gz', '.fastq', '.fastq.gz'):
            parser.error('This is not a file ending with either .fastq or \
 .fastq.gz. This flag requires the input of a fastq file.')
        else:
            return arg

    def is_valid_distance(self, parser, arg):
        #    parser.error('%s is not a positive number (e.g., a float, aka a \
# number with a decimal point)' % arg)
        if float(arg) < 0:
            parser.error("%s is not a positive number (e.g., a float, aka a \
    # number with a decimal point)" % arg)
        else:
            return arg


    def is_valid_int(self, parser, arg):
        isInt = True
        try:
            int(arg)
        except ValueError:
            isInt = False
        if isInt == False :
            parser.error("Input value is NOT an integer")
        elif arg.isnumeric() == False:
            parser.error("This is not a positive number")
        else:
            return arg

def argparser():
    description = """
    A script to parse the output from Mash and give you the most similar \
 species of interest.
    """
    parser = ParserWithErrors(description = description)
    #optional = parser._action_groups.pop()                                     #convert help menu to display as required and optional
    #required = parser.add_argument_group('required arguments')                 #flags. Change word parser in the below add_arguement
    #parser._action_groups.append(optional)                                     #to be either 'required' or 'optional' & uncomment out
                                                                                #these lines to take effect. Not necessary via NF.
    parser.add_argument("--database", "-d", required=True,
                        help="Pre-built Mash Sketch",
                        type=lambda x: parser.is_valid_mash(parser, x))
    parser.add_argument("--read1", "-r1", help="Input Read 1 (forward) file",
                        required=False,
                        type=lambda x: parser.is_valid_fastq(parser, x))
    parser.add_argument("--read2", "-r2", help="Input Read 2 (reverse) file",
                        required=False,
                        type=lambda x: parser.is_valid_fastq(parser, x))
    parser.add_argument("--max_dist", "-m", default=0.05,
                        help="reference fasta file path",
                        type=lambda x: parser.is_valid_distance(parser, x))
    parser.add_argument("--min_kmer", "-k", default=2,
                        help="Minimum copies of each kmer count to use",
                        type=lambda x: parser.is_valid_int(parser, x))
    parser.add_argument("--num_threads", "-t", default=2,
                        help="Number of computing threads to use",
                        type=lambda x: parser.is_valid_int(parser, x))
    return parser




## Add in next function


if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

inMash = args.database
inMaxDist = args.max_dist
inKmer = args.min_kmer
inThreads = args.num_threads
inRead1 = args.read1
inRead2 = args.read2
print(inMash)
print("This is read 1: ", inRead1)
print("This is read 2: ", inRead2)
print("This is the max distance: ", inMaxDist)
print("This is min_kmer :", inKmer)
print("This is default number of threads: ", inThreads)
