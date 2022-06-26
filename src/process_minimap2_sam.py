#!/usr/bin/env python3

# Name: process_spumoni_report.py
# Description: Process minimap2 SAM file and print out the reads that have
#              found in the database. Since MAPQ is not usuable with a pan-genome,
#              we will commit an alignment if all of its alignments are to 
#              same species.
# Date: June 26th, 2022

import os
import argparse
import pysam

def build_reference_name_sets(ref_paths):
    """ Builds a list of sets to be used to check for header names """
    ref_list = []
    for path in ref_paths:
        curr_set = set()
        with open(path, "r") as input_fd:
            all_lines = [curr_set.add(x.strip()) for x in input_fd.readlines()]
        ref_list.append(curr_set)
    return ref_list

def find_class_of_reference_name(ref_list, ref_name):
    """ Finds the class that a particular reference name occurs in """
    if ref_name is None: # handles unmapped reads
        return -1
        
    class_ids = []
    for i, name_set in enumerate(ref_list):
        if ref_name in name_set:
            class_ids.append(i+1)
    if len(class_ids) != 1:
        print(f"Error: this read hits {len(class_ids)} which is not expected.")
        exit(1)
    return class_ids[0]

def main(args):
    """ 
    Process the minimap2 SAM file and identify reads that have been found in
    the database ... meaning that they have alignments and all the alignments are
    to a unique species. Since with a pan-genome, there will multi-mapping so MAPQ
    will not be a good indication of a quality alignment.
    """
    ref_list = build_reference_name_sets(args.ref_names)
    sam_file = pysam.AlignmentFile(args.sam_file, "r")

    read_mappings = {}
    aligment_count = 0

    # Process all the alignments into our dictionary ...
    for read in sam_file.fetch():
        ref_class = find_class_of_reference_name(ref_list, read.reference_name)
        if read.query_name in read_mappings:
            read_mappings[read.query_name].append(ref_class)
        else:
            read_mappings[read.query_name] = [ref_class]
        aligment_count += 1
    
    sam_file.close()
    print(f"[log] found {aligment_count} alignments, and build a dictionary for read mappings.")

    # Determine which reads only have alignments to a single class ...
    output_fd = open(args.output_file, "w")
    committed_count = 0
    for key in read_mappings:
        if -1 not in read_mappings[key] and len(set(read_mappings[key])) == 1:
            committed_count += 1
            output_fd.write(f"{key}\n")

    output_fd.close()
    print(f"[log] found {committed_count} reads in the database, and wrote names to output file.")

def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a SAM file from minimap2, and determine which reads have been found in the database.")
    parser.add_argument("-s", dest="sam_file", help="path to minimap2 SAM file", required=True)
    parser.add_argument("-r", dest="ref_names", help="path to reference name lists", required=True, nargs='*')
    parser.add_argument("-o", dest="output_file", help="path to output file", required=True)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verify the correctness of the arguments """
    if not os.path.isfile(args.sam_file):
        print("Error: the path to SAM file is not valid.")
        exit(1)
    for path in args.ref_names:
        if not os.path.isfile(path):
            print("Error: at least one of the reference name lists is not valid.")
            exit(1)
    
if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)
