#!/usr/bin/env python

import argparse, sys, os

class ParserWithErrors(argparse.ArgumentParser):
    def error(self, message):
        print('{0}\n\n'.format(message))
        self.print_help()
        sys.exit(2)

    def is_valid_file(self, parser, arg):
        base, ext = os.path.splitext(arg)
        if ext not in ('.msh'):
            parser.error('This is not a file ending with .msh.\
 Did you generaete the mash sketch?')
        else:
            return arg

    def is_valid_distance(self, parser, arg):
        fdistIn = float(arg)
        if fdistIn <= 0:
            parser.error('%s is not a positive number (e.g., a float, aka a \
number with a decimal point)' % fdistIn)
        else:
            return arg

    def is_valid_Int(self, parser, arg):                                         #How to get it to print nice message as int?
        ivalIn = int(arg)                                                        #as is the error message is as follows
        if ivalIn <= 0:                                                          #invalid <lambda> value: '0.1' but I need
            parser.error('%s is not a positive number' % ivalIn)                 #but message should tell explicitly that
        else:                                                                    #it requires a whole number (integer)
            return arg

def argparser():
    description = """
    A script to parse the output from Mash and give you the most similar \
 species of interest.
    """
    parser = ParserWithErrors(description = description)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)

    required.add_argument("-d", "--database", required=True,
                        help="Pre-built Mash Sketch",
                        type=lambda x: parser.is_valid_file(parser, x))
    optional.add_argument("-m", "--max_dist", default=0.05,
                        help="reference fasta file path",
                        type=lambda x: parser.is_valid_distance(parser, x))
    optional.add_argument("-k", "--min_kmer",
                        help="Minimum copies of each kmer count to use",
                        type=lambda x: parser.is_valid_Int(parser, x))
    optional.add_argument("-t", "--num_threads", default=4,
                        help="Number of computing threads to use",
                        type=lambda x: parser.is_valid_Int(parser, x))
    return parser




## Add in next function


if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

inMash = args.database
inMaxDist = args.max_dist
inKmer = args.min_kmer
inThreads = args.num_threads
print(inMash)
print("This is the max distance: ", inMaxDist)
print("This is min_kmer :", inKmer)
print("This is default number of threads: ", inThreads)
