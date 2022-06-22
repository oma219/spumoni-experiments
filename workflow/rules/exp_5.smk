##################################################
# Name: exp_5.smk
# Description: Contains the workflow and methods
#              needed for experiment 5.
#
# Date: June 18, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_readset_for_exp5(wildcards):
    """
    Return a path to the read-set with real ONT reads
    from a Mock Community. This experiment assumes we
    only have 1 FASTQ file to process. 
    """
    file_list = []
    for data_file in os.listdir(f"data/full_read_set"):
        if data_file.endswith(".fastq"):
            file_list.append(f"{base_dir}/data/full_read_set/" + data_file)
    
    if len(file_list) > 1:
        print("Error: Expected only one input Fastq file.")
        exit(-1)

    return file_list

def get_all_zymo_genomes_exp5(wildcards):
    """
    Return paths to all eight of the genomes that were truly
    sequenced in the Zymo Mock Community experiment.
    """
    file_list = []
    for data_file in os.listdir(f"data/true_set"):
        if data_file.endswith(".fasta"):
            file_list.append(f"{base_dir}/data/true_set/" + data_file)

    if len(file_list) != 8:
        print("Error: Expected incorrect number of genomes found.")
        exit(-1)
    
    # Move the yeast to be last in list since it important for later rules
    yeast_file = f"{base_dir}/data/true_set/yeast_complete_genome.fasta"
    assert yeast_file in file_list

    file_list.remove(yeast_file)
    file_list.append(yeast_file)

    return file_list

def list_of_genomes_to_index_exp5(wildcards):
    """
    Returns all the paths to files to be included
    in the index.
    """
    input_files = []
    for num in range(1, num_datasets+1):
        for data_file in os.listdir(f"data/dataset_{num}"):
            if data_file.endswith(".fna"):
                input_files.append(f"{base_dir}/data/dataset_{num}/{data_file}")
    return input_files

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Sub-sample ONT reads that are long enough to classify confidently (>5000 bp) and
#              trim the beginning of the ONT to exclude the bases that will be used for experiment

rule sample_from_full_ONT_reads_exp5:
    input:
        get_readset_for_exp5
    output:
        "exp5_intermediate/step_1/sampled_reads.fastq"
    shell:
        """
        seqtk seq -L 5000 {input[0]} > {output}
        """

rule trim_ONT_reads_for_gold_standard_exp5:
    input:
        "exp5_intermediate/step_1/sampled_reads.fastq"
    output:
        "exp5_intermediate/step_2/trim_sampled_reads.fastq"
    shell:
        """
        seqtk trimfq -b 720 {input} > {output}
        """

# Section 2.2: Build a multi-FASTA file for the 8 genomes in the Zymo Mock Community and build
#              a minimap2 index for it.

rule build_ZymoMC_ref_exp5:
    input:
        get_all_zymo_genomes_exp5
    output:
        "exp5_intermediate/step_3/zymo_mc.fa"
    shell:
        """
        cat {input} > {output}
        """

rule build_minimap2_index_over_zymo_mc_exp5:
    input:
        "exp5_intermediate/step_3/zymo_mc.fa"
    output:
        "exp5_intermediate/step_4/zymo_mc.mmi"
    shell:
        """
        minimap2 -x map-ont -d {output} {input}
        """

# Section 2.3: Align all of the reads to develop the gold standard classifications, 
#              and filter them based on MAPQ. Generate files that contain the 
#              reference headers to separate out the microbial reads versus yeast 
#              reads.

rule align_all_reads_to_zymo_mc_exp5:
    input:
        "exp5_intermediate/step_2/trim_sampled_reads.fastq",
        "exp5_intermediate/step_4/zymo_mc.mmi"
    output:
        "exp5_intermediate/step_5/zymo_mc_read_aln.sam"
    shell:
        """
        minimap2 -p 0.6 -a {input[1]} {input[0]} > {output}
        """

rule filter_all_reads_by_mapq_exp5:
    input:
        "exp5_intermediate/step_5/zymo_mc_read_aln.sam"
    output:
        "exp5_intermediate/step_6/zymo_mc_read_aln.filtered.sam",
        "exp5_intermediate/step_6/zymo_mc_read_aln.filtered.sorted.sam"
    shell:
        """
        samtools view -u -q 30 {input} > {output[0]}
        samtools sort {output[0]} -o {output[1]}
        """

rule generate_reference_seq_name_list_exp5:
    input:
        get_all_zymo_genomes_exp5
    output:
        expand("exp5_ref_name_lists/class_{n}.txt", n=range(1, 9)),
        "exp5_ref_name_lists/README"
    shell:
        """
        i=1
        for file in {input}; do
            grep '^>' $file | awk '{{print substr($1,2)}}' > "exp5_ref_name_lists/class_$i.txt"
            echo $file >> "exp5_ref_name_lists/README"
            i=$((i+1))
        done
        """

# Section 2.4: Build indexes using SPUMONI (either promoted minimizers or
#              DNA minimizers) or minimap2 indexes

rule build_full_reference_to_index_exp5:
    input:
        list_of_genomes_to_index_exp5
    output:
        "exp5_full_ref/full_ref.fa"
    shell:
        """
        cat {input} > {output}
        """

rule build_spumoni_promoted_index_exp5:
    input:
        "exp5_full_ref/full_ref.fa"
    output:
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin"
    shell:
        """
        # Copy over reference file to specific folder, and build the index
        curr_ref_file="exp5_indexes/spumoni_promoted_k{wildcards.k}_w{wildcards.w}/full_ref.fa"
        log_file="exp5_indexes/spumoni_promoted_k{wildcards.k}_w{wildcards.w}/full_ref.fa.log"
        cp {input[0]} $curr_ref_file
        spumoni build -r $curr_ref_file -M -P -m -K {wildcards.k} -W {wildcards.w} &> $log_file
        """

rule build_spumoni_dna_index_exp5:
    input:
        "exp5_full_ref/full_ref.fa"
    output:
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa"
    shell:
        """
        # Copy over reference file to specific folder, and build the index
        curr_ref_file="exp5_indexes/spumoni_dna_k{wildcards.k}_w{wildcards.w}/full_ref.fa"
        log_file="exp5_indexes/spumoni_dna_k{wildcards.k}_w{wildcards.w}/full_ref.fa.log"
        cp {input[0]} $curr_ref_file
        spumoni build -r $curr_ref_file -M -P -t -K {wildcards.k} -W {wildcards.w} &> $log_file
        """

rule build_minimap2_index_exp5:
    input:
        "exp5_full_ref/full_ref.fa"
    output:
        "exp5_indexes/minimap2_index/full_ref.mmi"
    shell:
        """
        # Copy over reference file to specific folder, and build the index
        curr_ref_file="exp5_indexes/minimap2_index/full_ref.fa"
        log_file="exp5_indexes/minimap2_index/full_ref.fa.log"
        cp {input[0]} $curr_ref_file 
        minimap2 -x map-ont -d {output} {input} &> $log_file
        """





 