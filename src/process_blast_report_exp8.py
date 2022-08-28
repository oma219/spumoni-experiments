#!/usr/bin/env python3

# Name: process_blast_report_exp8.py
# Description: Takes in the output file from BLAST and summarizes which 
#              contigs appear to have significant contamination
# Date: August 27th, 2022

import os
import sys
import argparse
import random

def main(args):
    """
    Looks through the BLAST output file and determines which 
    contigs have significant alignments, and summarizes them
    in an output file along with average e-score.

    Output Format: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    """
    with open(args.input_file, "r") as input_fd:

        # Collect all the e-values for suspicious contigs
        suspicious_contigs = {}
        for line in input_fd:
            data_items = line.split()
            if float(data_items[10]) < 0.001:
                if data_items[0] not in suspicious_contigs:
                    suspicious_contigs[data_items[0]] = [float(data_items[10])]
                else:
                    suspicious_contigs[data_items[0]].append(float(data_items[10]))
        
        # Summarize each contig with average e-value
        for key in suspicious_contigs:
            avg = sum(suspicious_contigs[key])/len(suspicious_contigs[key])
            suspicious_contigs[key] = [avg]
        sorted_dict = dict(sorted(suspicious_contigs.items()))
        
        # Write out the contigs with suspicious e-values
        with open(args.output_dir + "alignment_report.txt", "w") as out_fd:
            out_fd.write("{:30}{:20}\n".format("name:", "average e-score:"))
            for key in sorted_dict:
                out_fd.write("{:30}{:20}\n".format(key, sorted_dict[key][0]))

def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a BLAST output file, and process it.")
    parser.add_argument("-i", dest="input_file", help="path to BLAST output file", required=True)
    parser.add_argument("-o", dest="output_dir", help="path to output directory", required=True)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verify the correctness of the arguments """
    if not os.path.isfile(args.input_file):
        print("Error: the input file is not valid.")
        exit(1)
    if not os.path.isdir(args.output_dir):
        print("Error: the output directory is not valid")
        exit(1)
    else:
        if args.output_dir[-1] != "/":
            args.output_dir += "/"
    
if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)