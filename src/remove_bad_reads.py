#!/usr/bin/env python3

# Name: remove_bad_reads.py
# Description: Removes reads with characters thar are non-ACGT characters
# Date: July 15th, 2022

import os
import sys
import argparse

def main(args):
    """ Parse out the reads, and output all reads where a majority of characters are ACGTs """
    with open(args.full_read_file, "r") as input_fd:
        all_lines = [x.strip() for x in input_fd.readlines()]
        non_common_reads = 0

        # iterate through all reads
        for i in range(0, len(all_lines), 2):
            # verify the header line is there
            assert ">" in all_lines[i], "assertion error: misformed read_id line"

            # determine if a majority of characters are non-common
            non_common = 0
            for ch in all_lines[i+1].upper():
                if ch not in ['A', 'C', 'G', 'T']:
                    non_common += 1
            
            # print out read if it is mostly common characters
            if non_common > 0:
                non_common_reads += 1
            else:
                print(f"{all_lines[i]}\n{all_lines[i+1]}")
        print(f"removed {non_common_reads} due to non-ACGT characters.", file=sys.stderr)
            
def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a FASTA file, remove any reads that contain mostly non-ACGTs, and re-output the read.")
    parser.add_argument("-i", dest="full_read_file", help="path to full input reads.", required=True)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verify the correctness of the arguments """
    if not os.path.isfile(args.full_read_file):
        print("Error: the input file is not valid.")
        exit(1)
    
if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)