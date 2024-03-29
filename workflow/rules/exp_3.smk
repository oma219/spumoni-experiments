##################################################
# Name: exp_3.smk
# Description: Contains the workflow and methods
#              needed for experiment 3.
#
# Date: May 7, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_genomes_in_dataset_exp3(wildcards):
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


def get_null_genomes_exp3(wildcards):
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

def get_all_report_files(wildcards):
    """
    Returns a list of all report files needed to generate the final 
    results file
    """
    file_list = ["exp3_results/full_index/short_positive_reads.fa.report",
                "exp3_results/full_index/short_null_reads.fa.report",
                "exp3_results/full_index/long_positive_reads.fa.report",
                "exp3_results/full_index/long_null_reads.fa.report",
                "exp3_results/full_index/index_stats.txt"]

    for wind_size in window_sizes_exp3:
        for index_type in ["promotion", "dna"]:
            for read_length in ["short", "long"]:
                for read_class in ["positive", "null"]:
                    file_list.append(f"exp3_results/k4_w{wind_size}/index_{index_type}/{read_length}_{read_class}_reads.fa.report")
            file_list.append(f"exp3_results/k4_w{wind_size}/index_{index_type}/index_stats.txt")
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Build file-list to use with SPUMONI

rule produce_list_of_genomes_exp3:
    input:
        get_genomes_in_dataset_exp3
    output:
        "exp3_file_list/dataset_file_list.txt"
    run:
        with open(output[0], "w") as fd:
            for file_name in input:
                fd.write(file_name + "\n")

# Section 2.2: Builds index using promoted-minimizer alphabet, or dna based index

rule build_index_promotion_minimizer_exp3:
    input:
        "exp3_file_list/dataset_file_list.txt"
    output:
        "exp3_promotion_index_k{k}_w{w}/output.bin.thrbv.spumoni",
        "exp3_promotion_index_k{k}_w{w}/output.bin.thrbv.ms",
        "exp3_promotion_index_k{k}_w{w}/output.bin"
    shell:
        """
        spumoni build -i {input[0]} \
                      -o exp3_promotion_index_k{wildcards.k}_w{wildcards.w}/output \
                      -M \
                      -P \
                      -m \
                      -K {wildcards.k} \
                      -W {wildcards.w} &> exp3_promotion_index_k{wildcards.k}_w{wildcards.w}/build.log
        """

rule build_index_dna_minimizer_exp3:
    input:
        "exp3_file_list/dataset_file_list.txt"
    output:
        "exp3_dna_index_k{k}_w{w}/output.fa.thrbv.spumoni",
        "exp3_dna_index_k{k}_w{w}/output.fa.thrbv.ms",
        "exp3_dna_index_k{k}_w{w}/output.fa"
    shell:
        """
        spumoni build -i {input[0]} \
                      -o exp3_dna_index_k{wildcards.k}_w{wildcards.w}/output \
                      -M \
                      -P \
                      -t \
                      -K {wildcards.k} \
                      -W {wildcards.w} &> exp3_dna_index_k{wildcards.k}_w{wildcards.w}/build.log
        """

rule build_full_index_exp3:
    input:
        "exp3_file_list/dataset_file_list.txt"
    output:
        "exp3_full_index/output.fa.thrbv.spumoni",
        "exp3_full_index/output.fa.thrbv.ms",
        "exp3_full_index/output.fa"
    shell:
        """
        spumoni build -i {input[0]} \
                      -o exp3_full_index/output \
                      -M \
                      -P \
                      -n  &> exp3_full_index/build.log
        """

# Section 2.3: Choose random positive genome, and simulate short or long reads.

rule simulate_short_positive_reads_exp3:
    input:
        get_genomes_in_dataset_exp3
    output:
        "exp3_reads/short/positive/positive_reads.fa"
    shell:
        """
        set +o pipefail;
        positive_genome=$(ls data/dataset_1/*.fna | shuf | head -n1)

        mason_simulator -ir $positive_genome \
                        -n {num_reads_exp3} \
                        -v \
                        -o {output[0]}.before_filter.fa \
                        --illumina-read-length 150 \
                        --illumina-prob-mismatch {illumina_mismatch_prob_exp3}

        python3 {repo_dir}/src/remove_bad_reads.py -i {output[0]}.before_filter.fa > {output[0]}
        rm {output[0]}.before_filter.fa
        """

rule simulate_long_positive_reads_exp3:
    input:
        get_genomes_in_dataset_exp3
    output:
        "exp3_reads/long/positive/positive_reads.fa"
    shell:
        """
        set +o pipefail;
        positive_genome=$(ls data/dataset_1/*.fna | shuf | head -n1)
        
        pbsim --depth 50.0 \
              --prefix exp3_reads/long/positive/positive_reads \
              --hmm_model {pbsim_model} \
              --accuracy-mean {long_read_acc_exp3} $positive_genome

        cat exp3_reads/long/positive/positive_reads_*.fastq > exp3_reads/long/positive/positive_reads.fastq
        rm exp3_reads/long/positive/positive_reads_*.fastq
        rm exp3_reads/long/positive/positive_reads_*.maf
        rm exp3_reads/long/positive/positive_reads_*.ref

        seqtk seq -a exp3_reads/long/positive/positive_reads.fastq > {output[0]}.full
        rm exp3_reads/long/positive/positive_reads.fastq
        num_lines=$(({num_reads_exp3} * 2))

        total_lines=$(wc -l {output[0]}.full | awk '{{print $1}}')
        if [ $total_lines -lt $num_lines ]; then
            echo "Error: not enough coverage to reach the reads needed.";
            exit;
        fi

        head -n $num_lines {output[0]}.full > {output[0]}.before_filter
        python3 {repo_dir}/src/remove_bad_reads.py -i {output[0]}.before_filter > {output[0]}
        rm {output[0]}.full {output[0]}.before_filter
        """

# Section 2.3: Copy over the provided null genome, and simulate short and long reads from it.

rule simulate_short_null_reads_exp3:
    input:
        get_null_genomes_exp3
    output:
        "exp3_reads/short/null/null_reads.fa"
    shell:
        """
        null_genome={input}
        
        mason_simulator -ir $null_genome \
                        -n {num_reads_exp3} \
                        -v \
                        -o {output[0]}.before_filter.fa \
                        --illumina-read-length 150 \
                        --illumina-prob-mismatch {illumina_mismatch_prob_exp3}
        
        python3 {repo_dir}/src/remove_bad_reads.py -i {output[0]}.before_filter.fa > {output[0]}
        rm {output[0]}.before_filter.fa
        """

rule simulate_long_null_reads_exp3:
    input:
        get_null_genomes_exp3
    output:
        "exp3_reads/long/null/null_reads.fa"
    shell:
        """
        null_genome={input}

        pbsim --depth 0.20 \
              --prefix exp3_reads/long/null/null_reads \
              --hmm_model {pbsim_model} \
              --accuracy-mean {long_read_acc_exp3} $null_genome

        cat exp3_reads/long/null/null_reads_*.fastq > exp3_reads/long/null/null_reads.fastq
        rm exp3_reads/long/null/null_reads_*.fastq
        rm exp3_reads/long/null/null_reads_*.maf
        rm exp3_reads/long/null/null_reads_*.ref

        seqtk seq -a exp3_reads/long/null/null_reads.fastq > {output[0]}.full
        rm exp3_reads/long/null/null_reads.fastq
        num_lines=$(({num_reads_exp3} * 2))

        total_lines=$(wc -l {output[0]}.full | awk '{{print $1}}')
        if [ $total_lines -lt $num_lines ]; then
            echo "Error: not enough coverage to reach the reads needed.";
            exit;
        fi

        head -n $num_lines {output[0]}.full > {output[0]}.before_filter
        python3 {repo_dir}/src/remove_bad_reads.py -i {output[0]}.before_filter > {output[0]}
        rm {output[0]}.full {output[0]}.before_filter
        """
    
# Section 2.4: Let SPUMONI classify both short and long reads using promotion-minimizer index

rule classify_using_promotion_index_exp3:
    input:
        "exp3_reads/short/positive/positive_reads.fa",
        "exp3_reads/short/null/null_reads.fa",
        "exp3_reads/long/positive/positive_reads.fa",
        "exp3_reads/long/null/null_reads.fa",
        "exp3_promotion_index_k{k}_w{w}/output.bin",
        "exp3_promotion_index_k{k}_w{w}/output.bin.thrbv.spumoni",
        "exp3_promotion_index_k{k}_w{w}/output.bin.thrbv.ms"
    output:
        "exp3_results/k{k}_w{w}/index_promotion/short_positive_reads.fa.report",
        "exp3_results/k{k}_w{w}/index_promotion/short_null_reads.fa.report",
        "exp3_results/k{k}_w{w}/index_promotion/long_positive_reads.fa.report",
        "exp3_results/k{k}_w{w}/index_promotion/long_null_reads.fa.report",
        "exp3_results/k{k}_w{w}/index_promotion/index_stats.txt"
    shell:
        """
        short_positive_reads="exp3_results/k{wildcards.k}_w{wildcards.w}/index_promotion/short_positive_reads.fa"
        short_null_reads="exp3_results/k{wildcards.k}_w{wildcards.w}/index_promotion/short_null_reads.fa"
        long_positive_reads="exp3_results/k{wildcards.k}_w{wildcards.w}/index_promotion/long_positive_reads.fa"
        long_null_reads="exp3_results/k{wildcards.k}_w{wildcards.w}/index_promotion/long_null_reads.fa"

        cp {input[0]} $short_positive_reads
        cp {input[1]} $short_null_reads
        cp {input[2]} $long_positive_reads
        cp {input[3]} $long_null_reads

        spumoni run -r exp3_promotion_index_k{wildcards.k}_w{wildcards.w}/output -p $short_positive_reads -P -c -m -K {wildcards.k} -W {wildcards.w}
        spumoni run -r exp3_promotion_index_k{wildcards.k}_w{wildcards.w}/output -p $short_null_reads -P -c -m -K {wildcards.k} -W {wildcards.w}
        spumoni run -r exp3_promotion_index_k{wildcards.k}_w{wildcards.w}/output -p $long_positive_reads -P -c -m -K {wildcards.k} -W {wildcards.w}
        spumoni run -r exp3_promotion_index_k{wildcards.k}_w{wildcards.w}/output -p $long_null_reads -P -c -m -K {wildcards.k} -W {wildcards.w}

        ls -l {input[5]} | awk '{{print $5}}' > {output[4]}

        #rm -r exp3_promotion_index_k{wildcards.k}_w{wildcards.w}/
        """

# Section 2.5: Let SPUMONI classify both short and long reads using DNA-minimizer index

rule classify_using_dna_index_exp3:
    input:
        "exp3_reads/short/positive/positive_reads.fa",
        "exp3_reads/short/null/null_reads.fa",
        "exp3_reads/long/positive/positive_reads.fa",
        "exp3_reads/long/null/null_reads.fa",
        "exp3_dna_index_k{k}_w{w}/output.fa",
        "exp3_dna_index_k{k}_w{w}/output.fa.thrbv.spumoni",
        "exp3_dna_index_k{k}_w{w}/output.fa.thrbv.ms"
    output:
        "exp3_results/k{k}_w{w}/index_dna/short_positive_reads.fa.report",
        "exp3_results/k{k}_w{w}/index_dna/short_null_reads.fa.report",
        "exp3_results/k{k}_w{w}/index_dna/long_positive_reads.fa.report",
        "exp3_results/k{k}_w{w}/index_dna/long_null_reads.fa.report",
        "exp3_results/k{k}_w{w}/index_dna/index_stats.txt"
    shell:
        """
        short_positive_reads="exp3_results/k{wildcards.k}_w{wildcards.w}/index_dna/short_positive_reads.fa"
        short_null_reads="exp3_results/k{wildcards.k}_w{wildcards.w}/index_dna/short_null_reads.fa"
        long_positive_reads="exp3_results/k{wildcards.k}_w{wildcards.w}/index_dna/long_positive_reads.fa"
        long_null_reads="exp3_results/k{wildcards.k}_w{wildcards.w}/index_dna/long_null_reads.fa"

        cp {input[0]} $short_positive_reads
        cp {input[1]} $short_null_reads
        cp {input[2]} $long_positive_reads
        cp {input[3]} $long_null_reads

        spumoni run -r exp3_dna_index_k{wildcards.k}_w{wildcards.w}/output -p $short_positive_reads -P -c -a -K {wildcards.k} -W {wildcards.w}
        spumoni run -r exp3_dna_index_k{wildcards.k}_w{wildcards.w}/output -p $short_null_reads -P -c -a -K {wildcards.k} -W {wildcards.w}
        spumoni run -r exp3_dna_index_k{wildcards.k}_w{wildcards.w}/output -p $long_positive_reads -P -c -a -K {wildcards.k} -W {wildcards.w}
        spumoni run -r exp3_dna_index_k{wildcards.k}_w{wildcards.w}/output -p $long_null_reads -P -c -a -K {wildcards.k} -W {wildcards.w}

        ls -l {input[5]} | awk '{{print $5}}' > {output[4]}

        #rm -r exp3_dna_index_k{wildcards.k}_w{wildcards.w}/
        """

# Section 2.6: Let SPUMONI classify both short and long reads using full index

rule classify_using_full_index_exp3:
    input:
        "exp3_reads/short/positive/positive_reads.fa",
        "exp3_reads/short/null/null_reads.fa",
        "exp3_reads/long/positive/positive_reads.fa",
        "exp3_reads/long/null/null_reads.fa",
        "exp3_full_index/output.fa",
        "exp3_full_index/output.fa.thrbv.spumoni",
        "exp3_full_index/output.fa.thrbv.ms"
    output:
        "exp3_results/full_index/short_positive_reads.fa.report",
        "exp3_results/full_index/short_null_reads.fa.report",
        "exp3_results/full_index/long_positive_reads.fa.report",
        "exp3_results/full_index/long_null_reads.fa.report",
        "exp3_results/full_index/index_stats.txt"
    shell:
        """
        short_positive_reads="exp3_results/full_index/short_positive_reads.fa"
        short_null_reads="exp3_results/full_index/short_null_reads.fa"
        long_positive_reads="exp3_results/full_index/long_positive_reads.fa"
        long_null_reads="exp3_results/full_index/long_null_reads.fa"

        cp {input[0]} $short_positive_reads
        cp {input[1]} $short_null_reads
        cp {input[2]} $long_positive_reads
        cp {input[3]} $long_null_reads

        spumoni run -r exp3_full_index/output -p $short_positive_reads -P -c -n 
        spumoni run -r exp3_full_index/output -p $short_null_reads -P -c -n 
        spumoni run -r exp3_full_index/output -p $long_positive_reads -P -c -n 
        spumoni run -r exp3_full_index/output -p $long_null_reads -P -c -n 

        ls -l {input[5]} | awk '{{print $5}}' > {output[4]}
        """

# Section 2.7: Generate report files and compute the confusion matrix

rule generate_analysis_file_exp3:
    input:
        get_all_report_files
    output:
        "exp3_analysis/exp3_total_results.csv",
        "exp3_analysis/exp3_full_index_results.csv",
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
            
            full_out_fd = open("exp3_analysis/exp3_full_index_results.csv", "w")
            full_out_fd.write("readlength,TP,FN,TN,FP,accuracy,pmlindexsize\n")

            # Process the report files for the full index 
            index_size = get_index_size(input[4])

            TP, FN = return_stats_from_file(input[0])
            FP, TN = return_stats_from_file(input[1])
            accuracy = (TP+TN)/(TP + FN + TN + FP)
            full_out_fd.write(f"short,{TP},{FN},{TN},{FP},{accuracy:.4f},{index_size}\n")

            TP, FN = return_stats_from_file(input[2])
            FP, TN = return_stats_from_file(input[3])
            accuracy = (TP+TN)/(TP + FN + TN + FP)
            full_out_fd.write(f"long,{TP},{FN},{TN},{FP},{accuracy:.4f},{index_size}\n")
            full_out_fd.close()

            # Process the report files for minimizer indexes
            out_fd.write("k,w,indextype,readlength,TP,FN,TN,FP,accuracy,pmlindexsize\n")
            input_list = input[5:]
            for i in range(0, len(input_list), 5):

                k = input_list[i].split("/")[1].split("_")[0][1:]
                w = input_list[i].split("/")[1].split("_")[1][1:]
                index_type = input_list[i].split("/")[2].split("_")[1]
                read_length = input_list[i].split("/")[3].split("_")[0]
                
                # Classify using the short reads, and print stats  
                index_size = get_index_size(input_list[i+4])

                TP, FN = return_stats_from_file(input_list[i])
                FP, TN = return_stats_from_file(input_list[i+1])
                accuracy = (TP+TN)/(TP + FN + TN + FP)
                out_fd.write(f"{k},{w},{index_type},short,{TP},{FN},{TN},{FP},{accuracy:.4f},{index_size}\n")

                # Classify using the long reads, and print stats
                TP, FN = return_stats_from_file(input_list[i+2])
                FP, TN = return_stats_from_file(input_list[i+3])
                accuracy = (TP+TN)/(TP + FN + TN + FP)
                out_fd.write(f"{k},{w},{index_type},long,{TP},{FN},{TN},{FP},{accuracy:.4f},{index_size}\n")


