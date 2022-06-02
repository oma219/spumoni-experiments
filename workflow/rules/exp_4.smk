##################################################
# Name: exp_3.smk
# Description: Contains the workflow and methods
#              needed for experiment 4.
#
# Date: June 2, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_genomes_in_dataset_exp4(wildcards):
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


def get_null_genomes_exp4(wildcards):
    """
    Returns a list of genomes used as null genomes. These are genomes that
    are considered to not be in the database. In this experiment, we 
    are assuming there is ONLY 1 null genome.
    """
    file_list = []
    for data_file in os.listdir("data/dataset_null"):
        if data_file.endswith(".fna") or data_file.endswith(".fa"):
            file_list.append(f"{base_dir}/data/dataset_null/" + data_file)

    assert len(file_list) == 1, "there should only be 1 null genome in experiment 3"
    return file_list[0]

def get_all_report_files_exp4(wildcards):
    """
    Returns a list of all report files needed to generate the final 
    results file
    """
    file_list = []
    for bin_size in bin_sizes_exp4:
        for wind_size in window_sizes_exp4:
            for index_type in ["promotion", "dna"]:
                for read_length in ["short", "long"]:
                    for read_class in ["positive", "null"]:
                        file_list.append(f"exp4_results/bin_{bin_size}/k4_w{wind_size}/index_{index_type}/{read_length}_{read_class}_reads.fa.report")
                file_list.append(f"exp4_results/bin_{bin_size}/k4_w{wind_size}/index_{index_type}/index_stats.txt")
    return file_list


####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Build file-list to use with SPUMONI

rule produce_list_of_genomes_exp4:
    input:
        get_genomes_in_dataset_exp4
    output:
        "exp4_file_list/dataset_file_list.txt"
    run:
        with open(output[0], "w") as fd:
            for file_name in input:
                fd.write(file_name + "\n")

# Section 2.2: Builds index using promoted-minimizer alphabet, or dna based index

rule build_index_promotion_minimizer_exp4:
    input:
        "exp4_file_list/dataset_file_list.txt"
    output:
        "exp4_promotion_index_k{k}_w{w}_bin{b}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp4_promotion_index_k{k}_w{w}_bin{b}/spumoni_full_ref.bin.thrbv.ms",
        "exp4_promotion_index_k{k}_w{w}_bin{b}/spumoni_full_ref.bin"
    run:
        shell("spumoni build -i {input[0]} -b exp4_promotion_index_k{wildcards.k}_w{wildcards.w}_bin{wildcards.b}/ -M -P -m -K {wildcards.k} -W {wildcards.w} -w {wildcards.b} &> exp4_promotion_index_k{wildcards.k}_w{wildcards.w}_bin{wildcards.b}/build.log")

rule build_index_dna_minimizer_exp4:
    input:
        "exp4_file_list/dataset_file_list.txt"
    output:
        "exp4_dna_index_k{k}_w{w}_bin{b}/spumoni_full_ref.fa.thrbv.spumoni",
        "exp4_dna_index_k{k}_w{w}_bin{b}/spumoni_full_ref.fa.thrbv.ms",
        "exp4_dna_index_k{k}_w{w}_bin{b}/spumoni_full_ref.fa"
    run:
        shell("spumoni build -i {input[0]} -b exp4_dna_index_k{wildcards.k}_w{wildcards.w}_bin{wildcards.b}/ -M -P -t -K {wildcards.k} -W {wildcards.w} -w {wildcards.b} &> exp4_dna_index_k{wildcards.k}_w{wildcards.w}_bin{wildcards.b}/build.log")


# Section 2.3: Choose random positive genome, and simulate short or long reads.

rule simulate_short_positive_reads_exp4:
    input:
        get_genomes_in_dataset_exp4
    output:
        "exp4_reads/short/positive/positive_reads.fa"
    shell:
        """
        positive_genome=$(ls data/dataset_1/*.fna | shuf | head -n1)
        mason_simulator -ir $positive_genome -n {num_reads_exp4} -v -o {output[0]} --illumina-read-length 150 --illumina-prob-mismatch {illumina_mismatch_prob_exp4}
        rm "$positive_genome.fai"
        """

rule simulate_long_positive_reads_exp4:
    input:
        get_genomes_in_dataset_exp4
    output:
        "exp4_reads/long/positive/positive_reads.fa"
    shell:
        """
        positive_genome=$(ls data/dataset_1/*.fna | shuf | head -n1)
        pbsim --depth 10.0 --prefix exp4_reads/long/positive/positive_reads --hmm_model {pbsim_model} --accuracy-mean {long_read_acc_exp4} $positive_genome

        cat exp4_reads/long/positive/positive_reads_*.fastq > exp4_reads/long/positive/positive_reads.fastq
        rm exp4_reads/long/positive/positive_reads_*.fastq
        rm exp4_reads/long/positive/positive_reads_*.maf
        rm exp4_reads/long/positive/positive_reads_*.ref

        seqtk seq -a exp4_reads/long/positive/positive_reads.fastq > {output[0]}.full
        rm exp4_reads/long/positive/positive_reads.fastq

        num_lines=$(({num_reads_exp4} * 2))
        head -n $num_lines {output[0]}.full > {output[0]}
        rm {output[0]}.full
        """

# Section 2.3: Copy over the provided null genome, and simulate short and long reads from it.

rule simulate_short_null_reads_exp4:
    input:
        get_null_genomes_exp4
    output:
        "exp4_reads/short/null/null_reads.fa"
    shell:
        """
        null_genome={input}
        mason_simulator -ir $null_genome -n {num_reads_exp4} -v -o {output[0]} --illumina-read-length 150 --illumina-prob-mismatch {illumina_mismatch_prob_exp4}
        rm "$null_genome.fai"
        """

rule simulate_long_null_reads_exp4:
    input:
        get_null_genomes_exp4
    output:
        "exp4_reads/long/null/null_reads.fa"
    shell:
        """
        null_genome={input}
        pbsim --depth 0.10 --prefix exp4_reads/long/null/null_reads --hmm_model {pbsim_model} --accuracy-mean {long_read_acc_exp4} $null_genome

        cat exp4_reads/long/null/null_reads_*.fastq > exp4_reads/long/null/null_reads.fastq
        rm exp4_reads/long/null/null_reads_*.fastq
        rm exp4_reads/long/null/null_reads_*.maf
        rm exp4_reads/long/null/null_reads_*.ref

        seqtk seq -a exp4_reads/long/null/null_reads.fastq > {output[0]}.full
        rm exp4_reads/long/null/null_reads.fastq

        num_lines=$(({num_reads_exp4} * 2))
        head -n $num_lines {output[0]}.full > {output[0]}
        rm {output[0]}.full
        """

# Section 2.4: Let SPUMONI classify both short and long reads using promotion-minimizer index

rule classify_using_promotion_index_exp4:
    input:
        "exp4_reads/short/positive/positive_reads.fa",
        "exp4_reads/short/null/null_reads.fa",
        "exp4_reads/long/positive/positive_reads.fa",
        "exp4_reads/long/null/null_reads.fa",
        "exp4_promotion_index_k{k}_w{w}_bin{b}/spumoni_full_ref.bin",
        "exp4_promotion_index_k{k}_w{w}_bin{b}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp4_promotion_index_k{k}_w{w}_bin{b}/spumoni_full_ref.bin.thrbv.ms"
    output:
        "exp4_results/bin_{b}/k{k}_w{w}/index_promotion/short_positive_reads.fa.report",
        "exp4_results/bin_{b}/k{k}_w{w}/index_promotion/short_null_reads.fa.report",
        "exp4_results/bin_{b}/k{k}_w{w}/index_promotion/long_positive_reads.fa.report",
        "exp4_results/bin_{b}/k{k}_w{w}/index_promotion/long_null_reads.fa.report",
        "exp4_results/bin_{b}/k{k}_w{w}/index_promotion/index_stats.txt"
    shell:
        """
        short_positive_reads="exp4_results/bin_{wildcards.b}/k{wildcards.k}_w{wildcards.w}/index_promotion/short_positive_reads.fa"
        short_null_reads="exp4_results/bin_{wildcards.b}/k{wildcards.k}_w{wildcards.w}/index_promotion/short_null_reads.fa"
        long_positive_reads="exp4_results/bin_{wildcards.b}/k{wildcards.k}_w{wildcards.w}/index_promotion/long_positive_reads.fa"
        long_null_reads="exp4_results/bin_{wildcards.b}/k{wildcards.k}_w{wildcards.w}/index_promotion/long_null_reads.fa"

        cp {input[0]} $short_positive_reads
        cp {input[1]} $short_null_reads
        cp {input[2]} $long_positive_reads
        cp {input[3]} $long_null_reads

        spumoni run -r {input[4]} -p $short_positive_reads -P -c -m -K {wildcards.k} -W {wildcards.w} -w {wildcards.b}
        spumoni run -r {input[4]} -p $short_null_reads -P -c -m -K {wildcards.k} -W {wildcards.w} -w {wildcards.b}
        spumoni run -r {input[4]} -p $long_positive_reads -P -c -m -K {wildcards.k} -W {wildcards.w} -w {wildcards.b}
        spumoni run -r {input[4]} -p $long_null_reads -P -c -m -K {wildcards.k} -W {wildcards.w} -w {wildcards.b}

        ls -l {input[5]} | awk '{{print $5}}' > {output[4]}

        rm -r exp4_promotion_index_k{wildcards.k}_w{wildcards.w}_bin{wildcards.b}/
        """

# Section 2.5: Let SPUMONI classify both short and long reads using DNA-minimizer index

rule classify_using_dna_index_exp4:
    input:
        "exp4_reads/short/positive/positive_reads.fa",
        "exp4_reads/short/null/null_reads.fa",
        "exp4_reads/long/positive/positive_reads.fa",
        "exp4_reads/long/null/null_reads.fa",
        "exp4_dna_index_k{k}_w{w}_bin{b}/spumoni_full_ref.fa",
        "exp4_dna_index_k{k}_w{w}_bin{b}/spumoni_full_ref.fa.thrbv.spumoni",
        "exp4_dna_index_k{k}_w{w}_bin{b}/spumoni_full_ref.fa.thrbv.ms"
    output:
        "exp4_results/bin_{b}/k{k}_w{w}/index_dna/short_positive_reads.fa.report",
        "exp4_results/bin_{b}/k{k}_w{w}/index_dna/short_null_reads.fa.report",
        "exp4_results/bin_{b}/k{k}_w{w}/index_dna/long_positive_reads.fa.report",
        "exp4_results/bin_{b}/k{k}_w{w}/index_dna/long_null_reads.fa.report",
        "exp4_results/bin_{b}/k{k}_w{w}/index_dna/index_stats.txt"
    shell:
        """
        short_positive_reads="exp4_results/bin_{wildcards.b}/k{wildcards.k}_w{wildcards.w}/index_dna/short_positive_reads.fa"
        short_null_reads="exp4_results/bin_{wildcards.b}/k{wildcards.k}_w{wildcards.w}/index_dna/short_null_reads.fa"
        long_positive_reads="exp4_results/bin_{wildcards.b}/k{wildcards.k}_w{wildcards.w}/index_dna/long_positive_reads.fa"
        long_null_reads="exp4_results/bin_{wildcards.b}/k{wildcards.k}_w{wildcards.w}/index_dna/long_null_reads.fa"

        cp {input[0]} $short_positive_reads
        cp {input[1]} $short_null_reads
        cp {input[2]} $long_positive_reads
        cp {input[3]} $long_null_reads

        spumoni run -r {input[4]} -p $short_positive_reads -P -c -a -K {wildcards.k} -W {wildcards.w} -w {wildcards.b}
        spumoni run -r {input[4]} -p $short_null_reads -P -c -a -K {wildcards.k} -W {wildcards.w} -w {wildcards.b}
        spumoni run -r {input[4]} -p $long_positive_reads -P -c -a -K {wildcards.k} -W {wildcards.w} -w {wildcards.b}
        spumoni run -r {input[4]} -p $long_null_reads -P -c -a -K {wildcards.k} -W {wildcards.w} -w {wildcards.b}

        ls -l {input[5]} | awk '{{print $5}}' > {output[4]}

        rm -r exp3_dna_index_k{wildcards.k}_w{wildcards.w}/
        """

# Section 2.6: Generate report files and compute the confusion matrix

rule generate_analysis_file_exp4:
    input:
        get_all_report_files_exp4
    output:
        "exp4_analysis/exp4_total_results.csv"
    run:
        def return_stats_from_file(report_file):
            """
            Returns the number of reads found and not found in the index 
            in the provided report file
            """
            with open(report_file, "r") as input_fd:
                output_stats = [0, 0]
                for line in input_fd:
                    if "FOUND" in line:
                        output_stats[0] += 1
                    elif "NOT_PRESENT" in line:
                        output_stats[1] += 1
            return output_stats
        
        def get_index_size(stats_file):
            """
            Extracts size of SPUMONI index stored in a file
            """
            with open(stats_file, "r") as fd:
                all_lines = fd.readlines()
                size = all_lines[0].strip()
                return size 
        
        with open(output[0], "w") as out_fd:
            
            # Process the report files for minimizer indexes
            out_fd.write("binsize,k,w,indextype,readlength,TP,FN,TN,FP,accuracy,pmlindexsize\n")
            
            for i in range(0, len(input), 5):

                bin_size = input[i].split("/")[1].split("_")[1]
                k = input[i].split("/")[2].split("_")[0][1:]
                w = input[i].split("/")[2].split("_")[1][1:]
                index_type = input[i].split("/")[3].split("_")[1]
                read_length = input[i].split("/")[4].split("_")[0]
                
                # Classify using the short reads, and print stats  
                index_size = get_index_size(input[i+4])

                TP, FN = return_stats_from_file(input[i])
                FP, TN = return_stats_from_file(input[i+1])
                accuracy = (TP+TN)/(TP + FN + TN + FP)
                out_fd.write(f"{bin_size},{k},{w},{index_type},short,{TP},{FN},{TN},{FP},{accuracy:.4f},{index_size}\n")

                # Classify using the long reads, and print stats
                TP, FN = return_stats_from_file(input[i+2])
                FP, TN = return_stats_from_file(input[i+3])
                accuracy = (TP+TN)/(TP + FN + TN + FP)
                out_fd.write(f"{bin_size},{k},{w},{index_type},long,{TP},{FN},{TN},{FP},{accuracy:.4f},{index_size}\n")
            