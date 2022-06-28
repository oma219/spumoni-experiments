##################################################
# Name: exp_6.smk
# Description: Contains the workflow and methods
#              needed for experiment 3.
#
# Date: June 27, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_genomes_in_dataset_exp6(wildcards):
    """ 
    Returns a list of genomes being used in the experiment. In this experiment, 
    we are assuming that we are just given 1 dataset. That is why we have 
    hard-coded the path to the dataset folder.
    """
    file_list = []
    for data_file in os.listdir(f"data/dataset_1"):
        if data_file.endswith(".fna"):
            file_list.append(f"{base_dir}/data/dataset_1/" + data_file)
    return file_list

def get_all_resource_files_exp6(wildcards):
    """
    Returns a list of all report files needed to generate the final 
    results file
    """
    file_list = ["exp6_results/full_index/short_positive_reads.fa.resources"]

    for wind_size in window_sizes_exp6:
        for index_type in ["promotion", "dna"]:
            for read_length in ["short"]:
                for read_class in ["positive"]:
                    file_list.append(f"exp6_results/k4_w{wind_size}/index_{index_type}/{read_length}_{read_class}_reads.fa.resources")
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Build file-list to use with SPUMONI

rule produce_list_of_genomes_exp6:
    input:
        get_genomes_in_dataset_exp6
    output:
        "exp6_file_list/dataset_file_list.txt"
    run:
        with open(output[0], "w") as fd:
            for file_name in input:
                fd.write(file_name + "\n")

# Section 2.2: Builds index using promoted-minimizer alphabet, dna-minimizers, or full text

rule build_index_promotion_minimizer_exp6:
    input:
        "exp6_file_list/dataset_file_list.txt"
    output:
        "exp6_promotion_index_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp6_promotion_index_k{k}_w{w}/spumoni_full_ref.bin.thrbv.ms",
        "exp6_promotion_index_k{k}_w{w}/spumoni_full_ref.bin"
    run:
        shell("spumoni build -i {input[0]} -b exp6_promotion_index_k{wildcards.k}_w{wildcards.w}/ -M -P -m -K {wildcards.k} -W {wildcards.w} &> exp6_promotion_index_k{wildcards.k}_w{wildcards.w}/build.log")

rule build_index_dna_minimizer_exp6:
    input:
        "exp6_file_list/dataset_file_list.txt"
    output:
        "exp6_dna_index_k{k}_w{w}/spumoni_full_ref.fa.thrbv.spumoni",
        "exp6_dna_index_k{k}_w{w}/spumoni_full_ref.fa.thrbv.ms",
        "exp6_dna_index_k{k}_w{w}/spumoni_full_ref.fa"
    run:
        shell("spumoni build -i {input[0]} -b exp6_dna_index_k{wildcards.k}_w{wildcards.w}/ -M -P -t -K {wildcards.k} -W {wildcards.w} &> exp6_dna_index_k{wildcards.k}_w{wildcards.w}/build.log")

rule build_full_index_exp6:
    input:
        "exp6_file_list/dataset_file_list.txt"
    output:
        "exp6_full_index/spumoni_full_ref.fa.thrbv.spumoni",
        "exp6_full_index/spumoni_full_ref.fa.thrbv.ms",
        "exp6_full_index/spumoni_full_ref.fa"
    run:
        shell("spumoni build -i {input[0]} -b exp6_full_index/ -M -P -n  &> exp6_full_index/build.log")

# Section 2.3: Choose random positive genome, simulate short reads

rule simulate_short_positive_reads_exp6:
    input:
        get_genomes_in_dataset_exp6
    output:
        "exp6_reads/short/positive/positive_reads.fa"
    shell:
        """
        set +o pipefail;
        positive_genome=$(ls data/dataset_1/*.fna | shuf | head -n1)
        mason_simulator -ir $positive_genome -n {num_reads_exp6} -v -o {output[0]} --illumina-read-length 150 --illumina-prob-mismatch {illumina_mismatch_prob_exp6}
        """

# Section 2.4: Let SPUMONI classify the short reads using promoted index, dna index, and full index

rule classify_using_promotion_index_exp6:
    input:
        "exp6_reads/short/positive/positive_reads.fa",
        "exp6_promotion_index_k{k}_w{w}/spumoni_full_ref.bin",
        "exp6_promotion_index_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp6_promotion_index_k{k}_w{w}/spumoni_full_ref.bin.thrbv.ms"
    output:
        "exp6_results/k{k}_w{w}/index_promotion/short_positive_reads.fa.report",
        "exp6_results/k{k}_w{w}/index_promotion/index_stats.txt",
        "exp6_results/k{k}_w{w}/index_promotion/short_positive_reads.fa.resources"
    shell:
        """
        short_positive_reads="exp6_results/k{wildcards.k}_w{wildcards.w}/index_promotion/short_positive_reads.fa"
        cp {input[0]} $short_positive_reads

        {time_prog} {time_format} --output={output[2]} spumoni run -r {input[1]} -p $short_positive_reads -P -c -m -K {wildcards.k} -W {wildcards.w}
        ls -l {input[2]} | awk '{{print $5}}' > {output[1]}
        """

rule classify_using_dna_index_exp6:
    input:
        "exp6_reads/short/positive/positive_reads.fa",
        "exp6_dna_index_k{k}_w{w}/spumoni_full_ref.fa",
        "exp6_dna_index_k{k}_w{w}/spumoni_full_ref.fa.thrbv.spumoni",
        "exp6_dna_index_k{k}_w{w}/spumoni_full_ref.fa.thrbv.ms"
    output:
        "exp6_results/k{k}_w{w}/index_dna/short_positive_reads.fa.report",
        "exp6_results/k{k}_w{w}/index_dna/index_stats.txt",
        "exp6_results/k{k}_w{w}/index_dna/short_positive_reads.fa.resources"
    shell:
        """
        short_positive_reads="exp6_results/k{wildcards.k}_w{wildcards.w}/index_dna/short_positive_reads.fa"
        cp {input[0]} $short_positive_reads

        {time_prog} {time_format} --output={output[2]} spumoni run -r {input[1]} -p $short_positive_reads -P -c -a -K {wildcards.k} -W {wildcards.w}
        ls -l {input[2]} | awk '{{print $5}}' > {output[1]}
        """

rule classify_using_full_index_exp6:
    input:
        "exp6_reads/short/positive/positive_reads.fa",
        "exp6_full_index/spumoni_full_ref.fa",
        "exp6_full_index/spumoni_full_ref.fa.thrbv.spumoni",
        "exp6_full_index/spumoni_full_ref.fa.thrbv.ms"
    output:
        "exp6_results/full_index/short_positive_reads.fa.report",
        "exp6_results/full_index/index_stats.txt",
        "exp6_results/full_index/short_positive_reads.fa.resources"
    shell:
        """
        short_positive_reads="exp6_results/full_index/short_positive_reads.fa"
        cp {input[0]} $short_positive_reads

        {time_prog} {time_format} --output={output[2]} spumoni run -r {input[1]} -p $short_positive_reads -P -c -n 
        ls -l {input[2]} | awk '{{print $5}}' > {output[1]}
        """

# Section 2.5: Generate all the needed result files, and summarize them into a handful of files

rule generate_analysis_files_exp6:
    input:
        get_all_resource_files_exp6
    output:
        "exp6_analysis/exp6_total_results.csv",
        "exp6_analysis/exp6_full_index_results.csv",
    run:
        # find the total time needed using full text
        with open(input[0], "r") as input_fd:
            all_lines = input_fd.readlines()
            assert len(all_lines) == 1

            total_time = all_lines[0].split()[5]
            with open(output[1], "w") as out_fd:
                out_fd.write(f"full_index,{total_time}\n")
        
        # find the total time for each minimizer based approach
        with open(output[0], "w") as out_fd:
            out_fd.write("k,w,type,totaltime\n")
            for path in input[1:]:
                k = path.split("/")[1].split("_")[0][1:]
                w = path.split("/")[1].split("_")[1][1:]
                index_type = path.split("/")[2].split("_")[1]

                total_time = ""
                with open(path, "r") as input_fd:
                    all_lines = input_fd.readlines()
                    assert len(all_lines) == 1
                    total_time = all_lines[0].split()[5]
                out_fd.write(f"{k},{w},{index_type},{total_time}\n")

