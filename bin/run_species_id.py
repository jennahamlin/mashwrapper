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
 number iwth a decimal point)' % fdistIn)
        else:
            return arg



def argparser():
    description = """
    A script to parse the output from Mash and give you the most similar \
 species of interest.
    """
    parser = ParserWithErrors(description = description)
    parser.add_argument("-d", "--database", required=True,
                        help="Pre-built Mash Sketch",
                        type=lambda x: parser.is_valid_file(parser, x))
    parser.add_argument("-m", "--max_dist", default=0.05,
                        help="reference fasta file path",
                        type=lambda x: parser.is_valid_distance(parser, x))
    return parser




## Add in next function


if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

inMash = args.database
inMaxDist = args.max_dist
print(inMash)
print(inMaxDist)
