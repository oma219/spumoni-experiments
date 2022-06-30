##################################################
# Name: exp_7.smk
# Description: Contains the workflow and methods
#              needed for experiment 7.
#
# Date: June 30, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def list_of_human_genomes_to_index_exp7(wildcards):
    """
    This method assumes that all the human genomes
    are stored in the folder called data/dataset_1
    """
    input_files = []
    for data_file in os.listdir(f"data/dataset_1"):
        if data_file.endswith(".fna"):
            input_files.append(f"{base_dir}/data/dataset_1/{data_file}")
    return input_files

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Sub-sample the microbiome reads, and generate the human reads through simulation

rule subsample_microbiome_reads_exp7:
    input:
        "data/microbiome_reads/SRR9847854_5000bp.fastq"
    output:
        "exp7_intermediate/step_1/null_reads.fa",
        "exp7_read_sets/null_reads.fa"
    shell:
        """
        seqtk seq -A {input} > {output[0]}
        num_lines=$(({num_pos_reads_exp7}*2))
        head -n $num_lines {output[0]} > {output[1]}
        """

rule simulate_long_human_reads_exp7:
    input:
        "data/human_genome/chm13.draft_v1.0.fasta"
    output:
        "exp7_intermediate/step_2/pos_reads.fastq",
        "exp7_intermediate/step_2/pos_reads.fa",
        "exp7_read_sets/pos_reads.fa"
    shell:
        """
        pbsim --depth 2.0 --prefix exp7_intermediate/step_2/pos_reads --hmm_model {pbsim_model} --accuracy-mean {long_read_acc_exp7} {input}

        # remove reads simulated from chrM
        rm exp7_intermediate/step_2/pos_reads_0024.fastq

        cat exp7_intermediate/step_2/pos_reads_*.fastq > {output[0]}
        rm exp7_intermediate/step_2/pos_reads_*.fastq
        rm exp7_intermediate/step_2/pos_reads_*.maf
        rm exp7_intermediate/step_2/pos_reads_*.ref

        seqtk seq -a {output[0]} > {output[1]}

        num_lines=$(({num_null_reads_exp7} * 2))
        head -n $num_lines {output[1]} > {output[2]}
        """

# Section 2.2: Build the complete reference, and then build the indexes needed for classification

rule build_human_genome_database_exp7:
    input:
        list_of_human_genomes_to_index_exp7
    output:
        "exp7_full_ref/human_database.fa"
    shell:
        """
        cat {input} > {output}
        """

rule build_spumoni_promoted_index_exp7:
    input:
        "exp7_full_ref/human_database.fa"
    output:
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.pmlnulldb"
    shell:
        """
        # Copy over reference file to specific folder, and build the index
        curr_ref_file="exp7_indexes/spumoni_promoted_k{wildcards.k}_w{wildcards.w}/full_ref.fa"
        log_file="exp7_indexes/spumoni_promoted_k{wildcards.k}_w{wildcards.w}/full_ref.fa.log"
        cp {input[0]} $curr_ref_file
        spumoni build -r $curr_ref_file -M -P -m -K {wildcards.k} -W {wildcards.w} &> $log_file
        """



