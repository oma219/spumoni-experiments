#!/usr/bin/env python3

# Name: grab_reads.py
# Description: This script is just meant to generate random
#              sequencing reads of a certain length, and a
#              certain number of them.
# Date: June 21, 2022

import argparse
import os
import random

def generate_reads_for_sequence(seq, num_reads, length):
    """ Takes in a single sequences and randomly generates reads from it """
    reads = []
    for i in range(0, num_reads):
        pos = random.randint(0, len(seq) - length)
        reads.append(seq[pos:pos+length])
    return reads

def main(args):
    """ Generates the specified reads, and prints to stdout """
    input_fasta = open(args.input_file, "r")

    # Go through each sequence and grab reads from it
    curr_seq = ""
    full_read_set = []
    for line in input_fasta:
        header_line = (">" in line)
        if header_line and len(curr_seq) > 0:
            full_read_set += generate_reads_for_sequence(curr_seq, args.num_reads, args.read_length)
            curr_seq = ""
        elif not header_line:
            curr_seq += line.strip()
    input_fasta.close()
    
    # Process the last sequence in FASTA file
    full_read_set += generate_reads_for_sequence(curr_seq, args.num_reads, args.read_length)

    # Sample without replacement a set of reads
    final_set = random.sample(full_read_set, args.num_reads)
    
    # Print out each read to stdout
    for i, read in enumerate(final_set):
        if args.fastq_file:
            random_quality = "#" * args.read_length
            print(f">seq_{i}\n{read}\n+\n{random_quality}")
        else:
            print(f">seq_{i}\n{read}")


def parse_arguments():
    """ Parses the command-line arguments, and returns them """
    parser = argparse.ArgumentParser(description="Randomly chooses reads from an input FASTA file of a certain length.")
    parser.add_argument("-i", dest="input_file", help="input file that we will read from", required=True)
    parser.add_argument("-n", dest="num_reads", help="number of reads that will be generated.", type=int, required=True)
    parser.add_argument("-l", dest="read_length", help="length of the sequencing reads.", type=int, required=True)
    parser.add_argument("-q", action="store_true", dest="fastq_file", default=False, help="output a fastq file (default: fasta)")
    args = parser.parse_args()
    return args

def check_args(args):
    """ Checks the command-line arguments and makes sure they are valid """
    if not os.path.isfile(args.input_file):
        print("Error: this input file does not exist.")
        exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    main(args)
