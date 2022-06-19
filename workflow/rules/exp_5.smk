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

    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Sub-sample ONT reads that are long enough to classify confidently (>5000 bp)

rule sample_from_full_ONT_reads_exp5:
    input:
        get_readset_for_exp5
    output:
        "exp5_step_1/sampled_reads.fastq"
    shell:
        """
        seqtk seq -L 5000 {input[0]} > {output}
        """

# Section 2.2: Trim the beginning of ONT reads (720 bp) to exclude the bases that will be used for experiment

rule trim_ONT_reads_for_gold_standard_exp5:
    input:
        "exp5_step_1/sampled_reads.fastq"
    output:
        "exp5_step_2/trim_sampled_reads.fastq"
    shell:
        """
        seqtk trimfq -b 720 {input} > {output}
        """

# Section 2.3: Build a multi-FASTA file for the 8 genomes in the Zymo Mock Community

rule build_ZymoMC_ref_exp5:
    input:
        get_all_zymo_genomes_exp5
    output:
        "exp5_step_3/zymo_mc.fa"
    shell:
        """
        cat {input} > {output}
        """

# Section 2.4: Build minimap2 index over the Zymo MC reference

rule build_minimap2_index_over_zymo_mc_exp5:
    input:
        "exp5_step_3/zymo_mc.fa"
    output:
        "exp5_step_4/zymo_mc.mmi"
    shell:
        """
        minimap2 -x map-ont -d exp5_step_4/zymo_mc.mmi exp5_step_3/zymo_mc.fa
        """

# Section 2.5: Align all of the reads to develop the gold standard classifications

rule align_all_reads_to_zymo_mc_exp5:
    input:
        "exp5_step_4/zymo_mc.mmi",
        "exp5_step_2/trim_sampled_reads.fastq"
    output:
        "exp5_step_5/zymo_mc_read_aln.sam"
    shell:
        """
        minimap2 -p 0.6 -a {input[0]} {input[1]} > {output}
        """

# Section 2.6: Filter out reads that have a high enough MAPQ

rule filter_all_reads_by_mapq_exp5:
    input:
        "exp5_step_5/zymo_mc_read_aln.sam"
    output:
        "exp5_step_6/zymo_mc_read_aln.filtered.sam",
        "exp5_step_6/zymo_mc_read_aln.filtered.sorted.sam"
    shell:
        """
        samtools view -u -q 30 {input} > {output[0]}
        samtools sort {output[0]} -o {output[1]}
        """
 