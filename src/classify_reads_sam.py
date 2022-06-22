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
    print("\n[log] building a dictionary of the read alignments")

    for read in sam_file.fetch():
        class_num = find_class_of_reference_name(ref_list, read.reference_name)
        if read.query_name in read_mappings:
            read_mappings[read.query_name].append(class_num)
        else:
            read_mappings[read.query_name] = [class_num]
        read_count += 1
    sam_file.close()

    # remove alignments that are ambiguous (multiple hits to different species)
    ambiguous_reads = 0
    reads_to_remove = []
    print(f"[log] found {len(read_mappings)} reads and a total of {read_count} alignments")
    print("\n[log] now will start to remove reads that have ambiguous placement")
    
    for key in read_mappings:
        assert len(read_mappings[key]) > 0
        if len(read_mappings[key]) > 1 and len(set(read_mappings[key])) > 1:
            reads_to_remove.append(key)
            ambiguous_reads += 1
    
    for key in reads_to_remove: # remove those reads
        read_mappings.pop(key)
    
    print(f"[log] removed {ambiguous_reads} reads, now we have {len(read_mappings)} reads left.")

    # quantify the number of alignments in each class
    num_classes = len(args.pos_refs) + len(args.null_refs)
    num_pos_classes = len(args.pos_refs)

    read_distribution = [0 for i in range(num_classes)]
    for key in read_mappings:
        curr_set = set(read_mappings[key])
        assert len(curr_set) == 1

        read_distribution[read_mappings[key][0]-1] += 1
    normalized_read_dist = [round(x/sum(read_distribution), 3) for x in read_distribution]

    print(f"\n[log] reads found for each class: {read_distribution}")
    print(f"[log] percentages found for each class: {normalized_read_dist}")

    # print out reads for each class separately
    sam_file = pysam.AlignmentFile(args.sam_file, "r")
    pos_file = open(args.output_dir + "pos_reads.fa", "w")
    null_file = open(args.output_dir + "null_reads.fa", "w")
    
    num_pos_classes = len(args.pos_refs)
    reads_written = set()

    for read in sam_file.fetch():
        if read.query_name in read_mappings and read.query_name not in reads_written:
            class_num = find_class_of_reference_name(ref_list, read.reference_name)
            if class_num <= num_pos_classes:
                pos_file.write(f">{read.query_name}\n{read.query_sequence}\n")
            else: # a read from null class
                null_file.write(f">{read.query_name}\n{read.query_sequence}\n")
            reads_written.add(read.query_name)
    pos_file.close()
    null_file.close()
    sam_file.close()

    print("\n[log] finished writing the reads to separate files for each class")

def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a SAM file and separate out the reads from each class.")
    parser.add_argument("-i", dest="sam_file", help="path to SAM file we will be parsing.", required=True)
    parser.add_argument("-p", dest="pos_refs", help="path to txt file with positive reference sequence names", required=True, nargs='*')
    parser.add_argument("-n", dest="null_refs", help="path to txt file with null reference sequence names", required=True, nargs='*')
    parser.add_argument("-o", dest="output_dir", help="path to output directory", required=True)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verifies certain aspects about the input parameters """
    if not os.path.isfile(args.sam_file):
        print("Error: the provided SAM file does not exist.")
        exit(-1)
    if not os.path.isdir(args.output_dir):
        print("Error: the provided output directory is not valid.")
        exit(1)
    elif args.output_dir[-1] != "/":
        args.output_dir += "/"
    
if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)
