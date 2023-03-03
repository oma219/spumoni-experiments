#!/usr/bin/env python3

# Name: process_lengths_exp8.py
# Description: Takes in a SPUMONI lengths file and computes report,
#              and csv files in order to plot the data.
# Date: July 22nd, 2022

import os
import sys
import argparse
import random
import numpy as np

def main(args):
    """ Process the lengths file and generate report and csvs """

    ## Step 1: Rank all the contigs by the longest length
    contaminated_bases = 0
    uncontaminated_bases = 0
    top_contigs = []
    full_info_dict = {}

    print("[log] starting to process the lengths file from SPUMONI")
    with open(args.lengths_file, "r") as input_fd:
        curr_seq = ""
        avg_dict = {}
        for i, line in enumerate(input_fd):
            if i % 2 == 0:
                curr_seq = line.strip()
            else:
                lengths = [int(x) for x in line.split()]
                avg = sum(lengths)/len(lengths)
                avg_dict[curr_seq] = avg
                full_info_dict[curr_seq] = [len(lengths), np.percentile(lengths, 25), np.percentile(lengths, 50), np.percentile(lengths, 75), np.percentile(lengths, 90), avg]

                # Add contig if its average is suspiously high
                if avg >= 2:
                    contaminated_bases += len(lengths)
                    top_contigs.append(curr_seq)
                else:
                    uncontaminated_bases += len(lengths)
    print("[log] finished finding avg length for each contig\n")

    ## Step 2: sort them, and assert that we have found some suspicious contigs
    sorted_dict = dict(sorted(avg_dict.items(), key=lambda item: item[1]))

    # top_contigs = []
    # for key in list(sorted_dict)[::-1]:
    #     if sorted_dict[key] >= 2:
    #         top_contigs.append(key)
    
    assert len(top_contigs) > 0, "assertion error: did not find any contigs with large average MS"
    print(f"[log] Out of {len(sorted_dict)} contigs, {len(top_contigs)} contigs appear to have contamination, a total of {contaminated_bases} bp")
    print(f"[log] There is a total of {uncontaminated_bases} 'uncontaminated' bp\n")

    ## Step 3: write out a report containing each contig and max length, as 
    ## well as a csv-file documenting the percentiles for each contig
    with open(args.output_dir + "assembly_length_report.txt", "w") as out_fd:
        out_fd.write("{:30}{:20}\n".format("name:", "avg length:"))
        for name in list(sorted_dict)[::-1]:
            out_fd.write("{:30}{:<20}\n".format(name, sorted_dict[name]))
    print("[log] wrote out report with all average lengths")

    with open(args.output_dir + "assembly_percentile_report.csv", "w") as out_fd:
        out_fd.write("contig,contiglength,25percent,50percent,75percent,90percent,mean\n")
        for key in dict(sorted(full_info_dict.items(), key=lambda item: item[1][-1])).keys():
            out_fd.write(",".join([key] + [str(x) for x in full_info_dict[key]]) + "\n")
    print("[log] wrote out report with percentile info for all contigs")

    ## Step 4: write out a csv-file for each of the top-ten contigs
    with open(args.lengths_file, "r") as input_fd:
        curr_seq = ""
        total_sampled = []
        for i, line in enumerate(input_fd):
            if i % 2 == 0:
                curr_seq = line.strip()
            else:
                if curr_seq in top_contigs:
                    lengths = [int(x) for x in line.split()]
                    index = top_contigs.index(curr_seq)
                    
                    with open(args.output_dir + f"top_contig_{index}_lengths.csv", "w") as out_fd:
                        for i, l in enumerate(lengths):
                            out_fd.write(f"{curr_seq[1:]},{i},{l}\n")
                else:
                    lengths = [int(x) for x in line.split()]
                    for i in range(0, len(lengths), 10000):
                        total_sampled.append(lengths[i])
    print("[log] wrote out the top-contigs to csv (those with avg>=2)\n")

    ## Step 5: write out a csv-file for the "regular" contigs
    with open(args.output_dir + "regular_contigs.csv", "w") as out_fd:
        for length in total_sampled:
            out_fd.write(f"regular,{length}\n")
    print("[log] finished writing sub-sampled csv file for regular contigs")

    ###############################################################################
    # Decided to comment this out, since decided to not include in figure anymore 
    ###############################################################################

    ## Step 6: choose a four random contig that were not suspicious and write out to csv file
    random_normal_contigs = random.sample(list(sorted_dict)[:-len(top_contigs)], k=4)

    # write out a csv-file for each of the top-ten contigs
    with open(args.lengths_file, "r") as input_fd:
        curr_seq = ""
        for i, line in enumerate(input_fd):
            if i % 2 == 0:
                curr_seq = line.strip()
            else:
                if curr_seq in random_normal_contigs:
                    lengths = [int(x) for x in line.split()]
                    index = random_normal_contigs.index(curr_seq)
                    
                    with open(args.output_dir + f"regular_contig_{index}_lengths.csv", "w") as out_fd:
                        for i, l in enumerate(lengths):
                            out_fd.write(f"{curr_seq[1:]},{i},{l}\n")
    print("[log] finished writing csv file for four random regular contigs\n")
    

def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a lengths file from exp8, and process it.")
    parser.add_argument("-i", dest="lengths_file", help="path to lengths file", required=True)
    parser.add_argument("-o", dest="output_dir", help="path to output directory", required=True)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verify the correctness of the arguments """
    if not os.path.isfile(args.lengths_file):
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