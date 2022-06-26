#!/usr/bin/env python3

# Name: extract_batch.py
# Description: Generate a FASTA file that represents the batch of data you would
#              get from an ONT sequencer to classify. It will print out the reads
#              corresponding to the method of classification. So SPUMONI will
#              operate on single batches of data, while minimap2 will operate
#              on successively longer prefixes of the read.
# Date: June 25th, 2022

import os
import argparse

def main(args):
    """ print out the batches to stdout """
    # Load the all the reads
    read_file = open(args.full_read_file, "r")
    all_reads = [x.strip() for x in read_file.readlines()]

    # Build set of reads to remove if given ...
    removed_reads = set()
    if len(args.reads_to_remove) > 0:
        with open(args.reads_to_remove, "r") as input_fd:
            all_lines = [removed_reads.add(x.strip()) for x in input_fd.readlines()]
            
    # Iterate through each read, extract batch and print
    for i in range(0, len(all_reads), 2):
        read_name = all_reads[i].split()[0][1:]
        read_seq = all_reads[i+1]
        read_batches = [read_seq[i:i+args.batch_size] for i in range(0, len(read_seq), args.batch_size)]

        # Output depends on the method
        if read_name not in removed_reads:
            if args.spumoni_method:
                print(f">{read_name}\n{read_batches[args.batch_num-1]}")
            elif args.alignment_method:
                curr_batch = "".join(read_batches[:args.batch_num])
                print(f">{read_name}\n{curr_batch}")
            
    read_file.close()


def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a FASTA file and parse out the batch of data requested.")
    parser.add_argument("-i", dest="full_read_file", help="path to full input reads.", required=True)
    parser.add_argument("-n", dest="batch_num", help="batch number from the ONT sequencer", required=True, default=0, type=int)
    parser.add_argument("-s", dest="batch_size", help="size of the batch of data in bp", required=True, default=180, type=int)
    parser.add_argument("-r", dest="reads_to_remove", help="reads to remove from printing", default="")
    parser.add_argument("--spumoni", dest="spumoni_method", action="store_true", help="using spumoni to classify the batch", default=False)
    parser.add_argument("--alignment", dest="alignment_method", action="store_true", help="using alignment to classify the batch", default=False)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verify the correctness of the arguments """
    if not os.path.isfile(args.full_read_file):
        print("Error: the input file is not valid.")
        exit(1)
    if args.batch_num not in [1, 2, 3, 4]:
        print("Error: the batch number provided is not valid.")
        exit(1)
    if args.batch_size < 0 or args.batch_size > 450:
        print("Error: the batch size in bp is not valid.")
        exit(1)
    if not args.spumoni_method and not args.alignment_method:
        print("Error: you need to specify exactly one method of classification")
        exit(1)
    if args.spumoni_method and args.alignment_method:
        print("Error: you need to specify exactly one method of classification")
        exit(1)
    

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)
