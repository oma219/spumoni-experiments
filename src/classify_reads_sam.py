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


def build_reference_name_sets(pos_refs, null_refs):
    """ Builds a list of sets to be used to check for header names """
    ref_list = []

    # add the sets for each class from the positive set of classes (microbial species)
    for pos_ref_file in pos_refs:
        curr_set = set()
        with open(pos_ref_file, "r") as input_fd:
            all_headers = [curr_set.add(x.strip()) for x in input_fd.readlines()]
        ref_list.append(curr_set)
    
    # add the sets for each class from the null set of classes (yeast)
    for null_ref_file in null_refs:
        curr_set = set()
        with open(null_ref_file, "r") as input_fd:
            all_headers = [curr_set.add(x.strip()) for x in input_fd.readlines()]
        ref_list.append(curr_set)
    
    return ref_list

def find_class_of_reference_name(ref_list, ref_name):
    """ Finds the class that a particular reference name occurs in """
    class_ids = []
    for i, name_set in enumerate(ref_list):
        if ref_name in name_set:
            class_ids.append(i+1)
    if len(class_ids) != 1:
        print(f"Error: this read hits {len(class_ids)} which is not expected.")
        exit(1)
    return class_ids[0]

def main(args):
    """ main method for script """
    # Open SAM file and create reference name sets
    sam_file = pysam.AlignmentFile(args.sam_file, "r")
    ref_list = build_reference_name_sets(args.pos_refs, args.null_refs)

    # add all the alignments to dictionary
    read_count = 0
    read_mappings = {}
    print("[log] building a dictionary of the read alignments")

    for read in sam_file.fetch():
        class_num = find_class_of_reference_name(ref_list, read.reference_name)
        if read.query_name in read_mappings:
            read_mappings[read.query_name].append(class_num)
        else:
            read_mappings[read.query_name] = [class_num]
        read_count += 1
        #if read_count > 20:
        #    break
    
    # remove alignments that are ambiguous (multiple hits to different species)
    ambiguous_reads = 0
    reads_to_remove = []
    print(f"[log] found {len(read_mappings)} alignments. removing reads that have ambiguous placement")
    
    for key in read_mappings:
        assert len(read_mappings[key]) > 0
        if len(read_mappings[key]) > 1 and len(set(read_mappings[key])) > 1:
            reads_to_remove.append(key)
            ambiguous_reads += 1
    
    for key in reads_to_remove:
        read_mappings.pop(key)
    
    print(ambiguous_reads)
    print(len(read_mappings))
    sam_file.close()

def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a SAM file and separate out the reads from each class.")
    parser.add_argument("-i", dest="sam_file", help="path to SAM file we will be parsing.", required=True)
    parser.add_argument("-p", dest="pos_refs", help="path to txt file with positive reference sequence names", required=True, nargs='*')
    parser.add_argument("-n", dest="null_refs", help="path to txt file with null reference sequence names", required=True, nargs='*')
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
