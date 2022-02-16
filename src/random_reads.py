###############################################################
#
# Name: random_reads.py
# Description: This script is just meant to generate random
#              sequencing reads of a certain length, and a
#              certain number of them.
#
# Author: Omar Ahmed
# Date: February 15, 2022
###############################################################

import argparse
import random 

def generate_reads(args):
    """ Generates the requested number of reads of a certain length """
    for read_num in range(args.num_reads):
        read = "".join(random.choices(['A', 'C', 'G', 'T'], k=args.read_length))
        print(f">read_{read_num}\n{read}")

def parse_arguments():
    """ Parses the command-line arguments, and returns them """
    parser = argparse.ArgumentParser(description="Random Read Generator - reads will be output to stdout.")
    parser.add_argument("-n", dest="num_reads", help="number of reads that will be generated.", type=int, required=True)
    parser.add_argument("-l", dest="read_length", help="length of the sequencing reads.", type=int, required=True)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    generate_reads(args)
