#!/usr/bin/env python3

# Name: process_spumoni_report.py
# Description: Process spumoni *.report file and print out the reads
#              that have been successfully classified as FOUND.
# Date: June 26th, 2022

import os
import argparse

def main(args):
    """ process the spumoni report files and determine which reads were found """

    # Remove the last batch file, and initialize the read_stats dictionary
    args.all_batches.remove(args.last_batch)
    read_stats = {}

    with open(args.last_batch, "r") as input_fd:
        all_lines = [x.strip() for x in input_fd.readlines()][1:]
        for curr_read in all_lines:
            line_split = curr_read.split()
            assert len(line_split) == 5, "Assertion Error: report file does not look as expected."

            read_stats[line_split[0]] = [float(line_split[3])]

    # Accumulate the previous KS-statistics from previous batches (only for reads in last batch)    
    for path in args.all_batches:
        with open(path, "r") as input_fd:
            all_lines = [x.strip() for x in input_fd.readlines()][1:]
            for curr_read in all_lines:
                line_split = curr_read.split()
                assert len(line_split) == 5, "Assertion Error: report file does not look as expected."

                if line_split[0] in read_stats:
                    read_stats[line_split[0]].append(int(line_split[3]))
    
    # Go through each read in last batch, and list the read names that are found
    with open(args.output_file, "w") as out_fd:
        for key in read_stats:
            num_above_threshold = sum(read_stats[key])
            if num_above_threshold/(len(read_stats[key]) + 0.0) >= 0.50:
                out_fd.write(f"{key}\n")
    
def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a spumoni *.report and print out the reads that have been found.")
    parser.add_argument("-i", dest="all_batches", help="path to spumoni report file(s)", required=True, nargs='*')
    parser.add_argument("-l", dest="last_batch", help="path to last batch report (tells which reads we should focus on)", required=True, type=str)
    parser.add_argument("-o", dest="output_file", help="path to output file", required=True)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verify the correctness of the arguments """
    assert args.last_batch in args.all_batches, "Assertion Error: last batch needs to be in all batches"
    for path in args.all_batches:
        if not os.path.isfile(path):
            print("Error: at least one of the paths is not valid.")
            exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)

