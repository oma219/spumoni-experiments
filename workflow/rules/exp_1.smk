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

def get_group_lists_exp1(wildcards):
    """ Returns a list of genomes given a dataset number """
    file_list = []
    for data_file in os.listdir(f"data/dataset_{wildcards.num}"):
        if data_file.endswith(".fna"):
            file_list.append(f"data/dataset_{wildcards.num}/" + data_file)
    return file_list

def get_files_for_analyzing_document_numbers_exp1(wildcards):
    """ Returns the input filelist for all files within a certain read-type/output-type group """
    output_type = wildcards.output
    read_type = wildcards.type
    input_filelist = [f"exp1_results/{output_type}/{read_type}/dataset_{num}/dataset_{num}_{read_type}_reads.fa.doc_numbers" for num in range(0, num_datasets+1)]
    return input_filelist

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Generate a list of genomes, and generate
#              references for each individual class for
#              simulating reads. 

rule produce_list_of_genomes_exp1:
    input:
        get_group_lists_exp1
    output:
        "exp1_intermediate/step_1/dataset_{num}/dataset_{num}_list.txt"
    run:
        with open(output[0], "w") as fd:
            for file_name in input:
                fd.write(file_name + "\n")

rule generate_full_file_list_exp1:
    input:
        expand("exp1_intermediate/step_1/dataset_{num}/dataset_{num}_list.txt", num=range(1,num_datasets+1))
    output:
        "exp1_file_list/genome_file_list.txt"
    run:
        with open(output[0], "w") as fd:
            for i, input_file in enumerate(input):
                with open(input_file, "r") as input_fd:
                    for line in input_fd:
                        fd.write(f"{base_dir}/{line.strip()} {i+1}\n")
            

# Section 2.2: Build the SPUMONI indexes for 
#              computing MS and PML based on the 
#              file list produced by earlier rule.

rule build_spumoni_index_exp1:
    input:
        "exp1_file_list/genome_file_list.txt"
    output:
        "exp1_indexes/spumoni_full_ref.bin.thrbv.ms",
        "exp1_indexes/spumoni_full_ref.bin.thrbv.spumoni"
    shell:
        """
        spumoni build -i {input[0]} -b exp1_indexes/ -M -P -m -K 4 -W 11 -d
        """

# Section 2.3: Simulate both long and short 
#              reads, first just more than needed, and 
#              then subset it down to the same number.

rule generate_raw_positive_short_reads_exp1:
    input:
        get_group_lists_exp1
    output:
        "exp1_raw_reads/illumina/dataset_{num}/dataset_{num}_illumina_reads.fq"
    shell:
        """
        set +o pipefail;
        positive_genome=$(ls data/dataset_{wildcards.num}/*.fna | shuf | head -n1)

        art_illumina -ss HS25 -i $positive_genome -na -l 150 -f 4.0 \
        -o exp1_raw_reads/illumina/dataset_{wildcards.num}/dataset_{wildcards.num}_illumina_reads
        """

rule generate_raw_positive_long_reads_exp1:
    input:
        get_group_lists_exp1
    output:
        "exp1_raw_reads/ont/dataset_{num}/dataset_{num}_ont_reads.fastq"
    shell:
        """
        set +o pipefail;
        positive_genome=$(ls data/dataset_{wildcards.num}/*.fna | shuf | head -n1)

        pbsim --depth 45.0 --prefix exp1_raw_reads/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 $positive_genome
        cat 'exp1_raw_reads/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads'*.fastq > {output}
        ls  'exp1_raw_reads/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_'*.fastq | xargs rm
        ls  'exp1_raw_reads/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads'*.ref | xargs rm
        ls  'exp1_raw_reads/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads'*.maf | xargs rm
        """

rule convert_short_reads_to_fasta_and_subset_exp1:
    input:
        "exp1_raw_reads/illumina/dataset_{num}/dataset_{num}_illumina_reads.fq"
    output:
        "exp1_reads/illumina/dataset_{num}/dataset_{num}_illumina_reads.fa"
    shell:
        """
        num_lines=$(({reads_per_dataset_exp1} * 4))
        head -n $num_lines {input} > exp1_reads/illumina/dataset_{wildcards.num}/dataset_{wildcards.num}_illumina_reads_subset.fq
        seqtk seq -a exp1_reads/illumina/dataset_{wildcards.num}/dataset_{wildcards.num}_illumina_reads_subset.fq > {output}

        if [ $(grep -c '>' {output}) != {reads_per_dataset_exp1} ]; then echo 'number of reads assertion failed.'; exit 1; fi
        rm exp1_reads/illumina/dataset_{wildcards.num}/dataset_{wildcards.num}_illumina_reads_subset.fq
        """

rule convert_long_reads_to_fasta_and_subset_exp1:
    input:
        "exp1_raw_reads/ont/dataset_{num}/dataset_{num}_ont_reads.fastq"
    output:
        "exp1_reads/ont/dataset_{num}/dataset_{num}_ont_reads.fa"
    shell:
        """
        num_lines=$(({reads_per_dataset_exp1} * 4))
        head -n $num_lines {input} > exp1_reads/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_subset.fq
        seqtk seq -a exp1_reads/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_subset.fq > {output}

        if [ $(grep -c '>' {output}) != {reads_per_dataset_exp1} ]; then echo 'number of reads assertion failed.'; exit 1; fi
        rm exp1_reads/ont/dataset_{wildcards.num}/dataset_{wildcards.num}_ont_reads_subset.fq
        """

# Section 2.4: Simulate random reads, both short and
#              long reads.

rule generate_random_short_reads_exp1:
    output:
        "exp1_reads/illumina/dataset_0/dataset_0_illumina_reads.fa"
    shell:
        "python3 {repo_dir}/src/random_reads.py -n {reads_per_dataset_exp1} -l 150 > {output}"

rule generate_random_long_reads_exp1:
    output:
        "exp1_reads/ont/dataset_0/dataset_0_ont_reads.fa"
    shell:
        "python3 {repo_dir}/src/random_reads.py -n {reads_per_dataset_exp1} -l 2000 > {output}"

# Section 2.5: Copy the read files over to a results directory
#              and then compute MS and PML for it.

rule copy_reads_to_result_folder_exp1:
    input:
        "exp1_reads/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa"
    output:
        "exp1_results/ms/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa",
        "exp1_results/pml/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa"
    shell:
        """
        cp {input} {output[0]}
        cp {input} {output[1]}
        """

rule run_spumoni_ms_on_reads_exp1:
    input:
        "exp1_results/ms/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa",
        "exp1_indexes/spumoni_full_ref.bin.thrbv.ms"
    output:
        "exp1_results/ms/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa.lengths",
        "exp1_results/ms/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa.doc_numbers"
    shell:
        "spumoni run -r exp1_indexes/spumoni_full_ref.bin -p {input[0]} -M -m -K 4 -W 11 -d"

rule run_spumoni_pml_on_reads_exp1:
    input:
        "exp1_results/pml/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa",
        "exp1_indexes/spumoni_full_ref.bin.thrbv.spumoni"
    output:
        "exp1_results/pml/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa.pseudo_lengths",
        "exp1_results/pml/{type}/dataset_{num}/dataset_{num}_{type}_reads.fa.doc_numbers"
    shell:
        "spumoni run -r exp1_indexes/spumoni_full_ref.bin -p {input[0]} -P -m -K 4 -W 11 -d"

# Section 2.6: Analyze the output *.doc_number files from
#              SPUMONI

rule analyze_document_numbers_exp1:
    input:
        get_files_for_analyzing_document_numbers_exp1
    output:
        "exp1_analysis/{type}_{output}_doc_analysis.csv"
    run:
        def get_percentages(line):
            counts = [0 for i in range(num_datasets_exp1)]
            doc_list = [int(x) for x in line.split()]
            for doc in doc_list:
                counts[doc] += 1
            return [x/sum(counts) for x in counts]

        output_fd = open(output[0], "w")
        headers = ['read_set'] + [f'class_{num}_percent' for num in range(1, num_datasets_exp1+1)]
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