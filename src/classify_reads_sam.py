#!/usr/bin/env python3

########################################################################
# Name: classify_reads_sam.py
#
# Description: Take in a SAM file to determine which reads are in each
#              class (truly from bacteria genomes vs. those from the
#              yeast genome).
# 
# Date: June 19th, 2022
#########################################################################

import os
import argparse
import pysam


def main(args):
    """ main method for script """
    sam_file = pysam.AlignmentFile(args.sam_file, "r")

    read_count = 0
    for read in sam_file.fetch():
        read_count += 1

        print(read.reference_name)
    
        if read_count > 10:
            break
    
    print(read_count)


    sam_file.close()

def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a SAM file and separate out the reads from each class.")
    parser.add_argument("-i", dest="sam_file", help="path to SAM file we will be parsing.", required=True)
    parser.add_argument("-p", dest="pos_refs", help="path to txt file with positive reference sequence names", required=True)
    parser.add_argument("-n", dest="null_refs", help="path to txt file with null reference sequence names", required=True)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verifies certain aspects about the input parameters """
    if not os.path.isfile(args.sam_file):
        print("Error: the provided SAM file does not exist.")
        exit(-1)

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)
