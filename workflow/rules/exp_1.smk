##################################################
# Name: exp_1.smk
# Description: Contains the workflow and methods
#              needed for experiment 1.
#
# Date: February 16, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_group_lists(wildcards):
    """ Retuns a list of genomes given a dataset number """
    file_list = []
    for data_file in os.listdir(f"data/dataset_{wildcards.num}"):
        if data_file.endswith(".fna"):
            file_list.append(f"data/dataset_{wildcards.num}/" + data_file)
    return file_list

def get_files_for_analyzing_document_numbers(wildcards):
    """ Returns the input filelist for all files within a certain read-type/output-type group """
    output_type = wildcards.output
    read_type = wildcards.type
    input_filelist = [f"results/{output_type}/{read_type}/dataset_{num}/dataset_{num}_{read_type}_reads.fa.doc_numbers" for num in range(0, num_datasets+1)]
    return input_filelist

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Generate a list of genomes, and generate
#              references for each individual class for
#              simulating reads. 

rule produce_list_of_genomes:
    input:
        get_group_lists
    output:
        temp("step_1/dataset_{num}/dataset_{num}_list.txt")
    run:
        with open(output[0], "w") as fd:
            for file_name in input:
                fd.write(file_name + "\n")

rule generate_full_file_list:
    input:
        expand("step_1/dataset_{num}/dataset_{num}_list.txt", num=range(1,num_datasets+1))
    output:
        "file_list/genome_file_list.txt"
    run:
        with open(output[0], "w") as fd:
            for i, input_file in enumerate(input):
                with open(input_file, "r") as input_fd:
                    for line in input_fd:
                        fd.write(f"{base_dir}/{line.strip()} {i+1}\n")
            
rule generate_one_seq_per_group:
    input:
        "step_1/dataset_{num}/dataset_{num}_list.txt"
    output:
        temp("step_2/dataset_{num}/dataset_{num}.fa")
    shell:
        "while read line; do cat $line >> {output}; done<{input} "

rule generate_rev_comp_for_group:
    input:
        "step_2/dataset_{num}/dataset_{num}.fa"
    output:
        temp("step_3/dataset_{num}/dataset_{num}_rev_comp.fa")
    shell:
        "seqtk seq -r -U {input} > {output}"

rule generate_combined_seq_for_group:
    input:
        "step_2/dataset_{num}/dataset_{num}.fa",
        "step_3/dataset_{num}/dataset_{num}_rev_comp.fa"
    output:
        "input_data/individual_datasets/dataset_{num}/dataset_{num}_combined.fa"
    run:
        shell("cat {input[0]} > {output}")
        shell("sed 's/^>/\>rev_comp_/' {input[1]} >> {output}")

"""
rule generate_full_reference:
    input:
        expand("input_data/individual_datasets/dataset_{num}/dataset_{num}_combined.fa", num=range(1, num_datasets+1))
    output:
        "input_data/combined_dataset/full_dataset.fa"
    shell:
        "for data_file in {input}; do cat $data_file >> {output[0]}; done"

rule build_fasta_index:
    input:
        "input_data/combined_dataset/full_dataset.fa"
    output:
        "input_data/combined_dataset/full_dataset.fa.fai"
    shell:
        "samtools faidx {input}"
"""

# Section 2.2: Build the SPUMONI indexes for 
#              computing MS and PML based on the 
#              file list produced by earlier rule.

rule build_spumoni_index:
    input:
        "file_list/genome_file_list.txt"
    output:
        "index/spumoni_full_ref.fa.thrbv.ms",
        "index/spumoni_full_ref.fa.thrbv.spumoni"
    run:
        shell("spumoni build -i {input[0]} -b index/ -M -P -f -d")

rule generate_raw_positive_short_reads:
    input:
        "input_data/individual_datasets/dataset_{num}/dataset_{num}_combined.fa"
    output:
        temp("step_4/illumina/dataset_{num}/dataset_{num}_illumina_reads.fq")
    shell:
        """
        art_illumina -ss HS25 -i {input} -na -l 150 -f 2.0 \
        -o step_4/illumina/dataset_{wildcards.num}/dataset_{wildcards.num}_illumina_reads
        """

# Section 2.3: Simulate both long and short 
#              reads, first just more than needed, and 
#              then subset it down to the same number.

rule generate_raw_positive_long_reads:
    input:
        "input_data/individual_datasets/dataset_{num}/dataset_{num}_combined.fa"
    output:
        temp("step_4/ont/dataset_{num}/dataset_{num}_ont_reads.fastq")
    run:
        shell("""
        pbsim --depth 25.0 --prefix step_4/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 {input}
        """)
        shell("cat 'step_4/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads'*.fastq > {output}")
        shell("ls 'step_4/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_'*.fastq | xargs rm")
        shell("ls 'step_4/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_'*.ref | xargs rm")
        shell("ls 'step_4/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_'*.maf | xargs rm")

rule convert_long_reads_to_fasta_and_subset:
    input:
        "step_4/ont/dataset_{num}/dataset_{num}_ont_reads.fastq"
    output:
        "reads/ont/dataset_{num}/dataset_{num}_ont_reads.fa"
    run:
        num_lines = num_reads_per_dataset * 4
        shell("head -n{num_lines} {input} > step_4/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_subset.fq")
        shell("seqtk seq -a step_4/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_subset.fq > {output}")
        shell("if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi")
        shell("rm step_4/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_subset.fq")

rule convert_short_reads_to_fasta_and_subset:
    input:
        "step_4/illumina/dataset_{num}/dataset_{num}_illumina_reads.fq"
    output:
        "reads/illumina/dataset_{num}/dataset_{num}_illumina_reads.fa"
    run:
        num_lines = num_reads_per_dataset * 4
        shell("head -n{num_lines} {input} > step_4/illumina/dataset_{wildcards.num}/dataset_{wildcards.num}_illumina_reads_subset.fq")
        shell("seqtk seq -a step_4/illumina/dataset_{wildcards.num}/dataset_{wildcards.num}_illumina_reads_subset.fq > {output}")
        shell("if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi")
        shell("rm step_4/illumina/dataset_{wildcards.num}/dataset_{wildcards.num}_illumina_reads_subset.fq")

# Section 2.4: Simulate random reads, both short and
#              long reads.

rule generate_random_short_reads:
    output:
        "reads/illumina/dataset_0/dataset_0_illumina_reads.fa"
    shell:
        "python {repo_dir}/src/random_reads.py -n {num_reads_per_dataset} -l 150 > {output}"

rule generate_random_long_reads:
    output:
        "reads/ont/dataset_0/dataset_0_ont_reads.fa"
    shell:
        "python {repo_dir}/src/random_reads.py -n {num_reads_per_dataset} -l 2000 > {output}"

# Section 2.5: Copy the read files over to a results directory
#              and then compute MS and PML for it.

rule copy_reads_to_result_folder:
    input:
        "reads/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa"
    output:
        "results/ms/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa",
        "results/pml/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa"
    run:
        shell("cp {input} {output[0]}")
        shell("cp {input} {output[1]}")

rule run_spumoni_ms_on_reads:
    input:
        "results/ms/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa",
        "index/spumoni_full_ref.fa",
        "index/spumoni_full_ref.fa.thrbv.ms"
    output:
        "results/ms/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa.lengths",
        "results/ms/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa.doc_numbers"
    shell:
        "spumoni run -r {input[1]} -p {input[0]} -M -f -d"

rule run_spumoni_pml_on_reads:
    input:
        "results/pml/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa",
        "index/spumoni_full_ref.fa",
        "index/spumoni_full_ref.fa.thrbv.spumoni"
    output:
        "results/pml/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa.pseudo_lengths",
        "results/pml/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa.doc_numbers"
    shell:
        "spumoni run -r {input[1]} -p {input[0]} -P -f -d"

# Section 2.6: Analyze the output *.doc_number files from
#              SPUMONI

rule analyze_document_numbers:
    input:
        get_files_for_analyzing_document_numbers
    output:
        "analysis/{type}_{output}_doc_analysis.csv"
    run:
        def get_percentages(line):
            counts = [0 for i in range(num_datasets)]
            doc_list = [int(x) for x in line.split()]
            for doc in doc_list:
                counts[doc] += 1
            return [x/sum(counts) for x in counts]

        output_fd = open(output[0], "w")
        headers = ['read_set'] + [f'class_{num}_percent' for num in range(1, num_datasets+1)]
        output_fd.write(",".join(headers) + "\n")

        for input_file in input:
            read_set = input_file.split("/")[3]

            # Calculate the average % across documents
            with open(input_file, "r") as input_fd:
                all_lines = input_fd.readlines()
 
                total_percents = [0.0 for x in range(num_datasets)]
                num_reads = 0
                for line in all_lines:
                    if ">" not in line:
                        num_reads += 1
                        percents = get_percentages(line) 
                        for i, val in enumerate(percents):
                            total_percents[i] += val
                total_percents = [round(x/num_reads, 4) for x in total_percents]       
            output_fd.write(",".join([read_set] + [str(x) for x in total_percents]) + "\n")
        output_fd.close()