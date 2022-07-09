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

def get_all_genomes_to_index_for_specific_species(wildcards):
    """
    Returns all the paths to genomes included in index for
    a particular species.
    """
    input_files = []
    for data_file in os.listdir(f"data/dataset_{wildcards.num}"):
        if data_file.endswith(".fna"):
            input_files.append(f"{base_dir}/data/dataset_{wildcards.num}/{data_file}")
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

rule generate_separate_read_files_for_each_class_exp5:
    input:
        expand("exp5_ref_name_lists/class_{n}.txt", n=range(1, 9)),
        "exp5_intermediate/step_6/zymo_mc_read_aln.filtered.sorted.sam"
    output:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa"
    shell:
        """
        python3 {repo_dir}/src/classify_reads_sam.py \
        -i exp5_intermediate/step_6/zymo_mc_read_aln.filtered.sorted.sam \
        -p exp5_ref_name_lists/class_{{1..7}}.txt \
        -n exp5_ref_name_lists/class_8.txt \
        -o exp5_read_sets/ \
        -r exp5_intermediate/step_1/sampled_reads.fastq
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
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.pmlnulldb"
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
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.thrbv.spumoni",
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.pmlnulldb"
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
        minimap2 -x map-ont -d {output} --split-prefix="exp5_indexes/minimap2_index/split_" {input} &> $log_file
        """

# Section 2.5: Extract a batch of data that will be processed by SPUMONI or minimap2.
#              For SPUMONI, it will just be the current batch while for minimap2 it will
#              a consistently longer prefix of the read.

rule extract_first_batch_of_data_for_spumoni_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa"
    output:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp5_read_sets/{wildcards.class}_reads.fa \
        -n 1 -s 180 -a --spumoni > {output}
        """

rule extract_first_batch_of_data_for_minimap2_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa"
    output:
        "exp5_results/minimap2/{class}_reads/batch_1/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp5_read_sets/{wildcards.class}_reads.fa \
        -n 1 -s 180 -a --alignment > {output}
        """

# Section 2.6: Classify the first batch of data using SPUMONI and minimap2

rule classify_first_batch_using_spumoni_promoted_exp5:
    input:
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa"
    output:
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.report",
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} -P -m -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_first_batch_using_spumoni_dna_exp5:
    input:
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa"
    output:
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.report",
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} -P -a -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_first_batch_using_minimap2_exp5:
    input:
        "exp5_indexes/minimap2_index/full_ref.mmi",
        "exp5_results/minimap2/{class}_reads/batch_1/curr_batch.fa"
    output:
        "exp5_results/minimap2/{class}_reads/batch_1/curr_batch.sam",
        "exp5_results/minimap2/{class}_reads/batch_1/curr_batch.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} minimap2 -t 0 -a {input[0]} --split-prefix="exp5_indexes/minimap2_index/split_" {input[1]} -o {output[0]}
        """

# Section 2.7: Parse out the reference names in each class, and store them to be used when
#              we try to identify what species each read is mapping to.

rule parse_ref_names_from_full_dataset_exp5:
    input:
        get_all_genomes_to_index_for_specific_species
    output:
        "exp5_index_ref_name_lists/class_{num}.txt"
    shell:
        """
        for file in {input}; do
            grep '^>' $file | awk '{{print substr($1,2)}}' >> "exp5_index_ref_name_lists/class_{wildcards.num}.txt"
        done
        """
    
# Section 2.8: Determine which reads each method did not find in the database in order
#              to keep them to classify in next batch

rule determine_reads_to_keep_for_second_batch_spumoni_exp5:
    input:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.report"
    output:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_second_batch_minimap2_exp5:
    input:
        "exp5_results/minimap2/{class}_reads/batch_1/curr_batch.sam",
        expand("exp5_index_ref_name_lists/class_{num}.txt", num=range(1, 8))
    output:
        "exp5_results/minimap2/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_minimap2_sam.py -s {input[0]} -o {output} -r exp5_index_ref_name_lists/class_{{1..7}}.txt
        """

# Section 2.9: Extract the second batch of data for minimap2 and SPUMONI. This
#              time we will take into account the reads that need to be excluded
#              since they have already been classified.

rule extract_second_batch_of_data_for_spumoni_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    output:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp5_read_sets/{wildcards.class}_reads.fa \
        -n 2 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_second_batch_of_data_for_minimap2_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/minimap2/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    output:
        "exp5_results/minimap2/{class}_reads/batch_2/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp5_read_sets/{wildcards.class}_reads.fa \
        -n 2 -s 180 -k {input[2]} --alignment > {output}
        """

# Section 2.10: Classify the second batch of data given using SPUMONI and minimap2

rule classify_second_batch_using_spumoni_promoted_exp5:
    input:
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa"
    output:
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.report",
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} -P -m -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_second_batch_using_spumoni_dna_exp5:
    input:
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa"
    output:
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.report",
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} -P -a -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_second_batch_using_minimap2_exp5:
    input:
        "exp5_indexes/minimap2_index/full_ref.mmi",
        "exp5_results/minimap2/{class}_reads/batch_2/curr_batch.fa"
    output:
        "exp5_results/minimap2/{class}_reads/batch_2/curr_batch.sam",
        "exp5_results/minimap2/{class}_reads/batch_2/curr_batch.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} minimap2 -t 0 -a {input[0]} --split-prefix="exp5_indexes/minimap2_index/split_" {input[1]} -o {output[0]}
        """

# Section 2.11: Determine which reads each method did not find yet in order to 
#               keep them for classifying the next batch

rule determine_reads_to_keep_for_third_batch_spumoni_exp5:
    input:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.report"
    output:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_third_batch_minimap2_exp5:
    input:
        "exp5_results/minimap2/{class}_reads/batch_2/curr_batch.sam",
        expand("exp5_index_ref_name_lists/class_{num}.txt", num=range(1, 8))
    output:
        "exp5_results/minimap2/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_minimap2_sam.py -s {input[0]} -o {output} -r exp5_index_ref_name_lists/class_{{1..7}}.txt
        """

# Section 2.12: Extract the third batch of data for minimap2 and SPUMONI. This
#               time we will take into account the reads that need to be excluded
#               since they have already been classified.

rule extract_third_batch_of_data_for_spumoni_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    output:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp5_read_sets/{wildcards.class}_reads.fa\
        -n 3 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_third_batch_of_data_for_minimap2_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/minimap2/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    output:
        "exp5_results/minimap2/{class}_reads/batch_3/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp5_read_sets/{wildcards.class}_reads.fa \
        -n 3 -s 180 -k {input[2]} --alignment > {output}
        """

# Section 2.13: Classify the third batch of data given using SPUMONI and minimap2

rule classify_third_batch_using_spumoni_promoted_exp5:
    input:
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa"
    output:
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.report",
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} -P -m -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_third_batch_using_spumoni_dna_exp5:
    input:
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa"
    output:
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.report",
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} -P -a -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_third_batch_using_minimap2_exp5:
    input:
        "exp5_indexes/minimap2_index/full_ref.mmi",
        "exp5_results/minimap2/{class}_reads/batch_3/curr_batch.fa"
    output:
        "exp5_results/minimap2/{class}_reads/batch_3/curr_batch.sam",
        "exp5_results/minimap2/{class}_reads/batch_3/curr_batch.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} minimap2 -t 0 -a {input[0]} --split-prefix="exp5_indexes/minimap2_index/split_" {input[1]} -o {output[0]}
        """

# Section 2.14: Determine which reads each method did not find yet in order to 
#               keep those reads for the next batch.

rule determine_reads_to_keep_for_fourth_batch_spumoni_exp5:
    input:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.report"
    output:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_fourth_batch_minimap2_exp5:
    input:
        "exp5_results/minimap2/{class}_reads/batch_3/curr_batch.sam",
        expand("exp5_index_ref_name_lists/class_{num}.txt", num=range(1, 8))
    output:
        "exp5_results/minimap2/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_minimap2_sam.py -s {input[0]} -o {output} -r exp5_index_ref_name_lists/class_{{1..7}}.txt
        """

# Section 2.15: Extract the fouth batch of data for minimap2 and SPUMONI. This
#               time we will take into account the reads that need to be excluded
#               since they have already been classified.

rule extract_fourth_batch_of_data_for_spumoni_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    output:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp5_read_sets/{wildcards.class}_reads.fa\
        -n 4 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_fourth_batch_of_data_for_minimap2_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/minimap2/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    output:
        "exp5_results/minimap2/{class}_reads/batch_4/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp5_read_sets/{wildcards.class}_reads.fa \
        -n 4 -s 180 -k {input[2]} --alignment > {output}
        """

# Section 2.16: Classify the fourth batch of data given using SPUMONI and minimap2

rule classify_fourth_batch_using_spumoni_promoted_exp5:
    input:
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa"
    output:
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.report",
        "exp5_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} -P -m -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_fourth_batch_using_spumoni_dna_exp5:
    input:
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa"
    output:
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.report",
        "exp5_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} -P -a -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_fourth_batch_using_minimap2_exp5:
    input:
        "exp5_indexes/minimap2_index/full_ref.mmi",
        "exp5_results/minimap2/{class}_reads/batch_4/curr_batch.fa"
    output:
        "exp5_results/minimap2/{class}_reads/batch_4/curr_batch.sam",
        "exp5_results/minimap2/{class}_reads/batch_4/curr_batch.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} minimap2 -t 0 -a {input[0]} --split-prefix="exp5_indexes/minimap2_index/split_" {input[1]} -o {output[0]}
        """

# Section 2.17: Determine which reads each method did not find yet in order to 
#               keep those reads for the next batch.

rule determine_reads_to_keep_for_fifth_batch_spumoni_exp5:
    input:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.report"
    output:
        "exp5_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_fifth_batch_minimap2_exp5:
    input:
        "exp5_results/minimap2/{class}_reads/batch_4/curr_batch.sam",
        expand("exp5_index_ref_name_lists/class_{num}.txt", num=range(1, 8))
    output:
        "exp5_results/minimap2/{class}_reads/batch_4/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_minimap2_sam.py -s {input[0]} -o {output} -r exp5_index_ref_name_lists/class_{{1..7}}.txt
        """

# Section 2.18: Analyze each classification approach individually in terms of 
#               TP, FN, FP, TN, time, peak memory, and index size

rule analyze_spumoni_promoted_results_into_csv_line_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp5_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp5_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_1/curr_batch.fa.resources",
        "exp5_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_2/curr_batch.fa.resources",
        "exp5_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_3/curr_batch.fa.resources",
        "exp5_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_4/curr_batch.fa.resources",
        "exp5_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_1/curr_batch.fa.resources",
        "exp5_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_2/curr_batch.fa.resources",
        "exp5_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_3/curr_batch.fa.resources",
        "exp5_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_4/curr_batch.fa.resources",
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp5_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.pmlnulldb"
    output:
        "exp5_results/spumoni_promoted_k{k}_w{w}/analysis.csv"
    shell:
        """
        # Compute the confusion matrix stats, and accuracy
        num_pos_reads=$(grep -c ">" {input[0]})
        num_null_reads=$(grep -c ">" {input[1]})

        FN=$(wc -l {input[2]} | awk '{{print $1}}')
        TN=$(wc -l {input[3]} | awk '{{print $1}}')
        TP=$(($num_pos_reads-$FN))
        FP=$(($num_null_reads-$TN))

        accuracy=$(echo "scale=4; ($TP+$TN)/($TP+$FN+$FP+$TN)" | bc)
        sensitivity=$(echo "scale=4; ($TP)/($TP+$FN)" | bc)
        specificity=$(echo "scale=4; ($TN)/($FP+$TN)" | bc)

        # Compute time used in each stage
        time_1=$(cat {input[4]} | awk '{{print $6}}')
        time_2=$(cat {input[5]} | awk '{{print $6}}')
        time_3=$(cat {input[6]} | awk '{{print $6}}')
        time_4=$(cat {input[7]} | awk '{{print $6}}')
        time_5=$(cat {input[8]} | awk '{{print $6}}')
        time_6=$(cat {input[9]} | awk '{{print $6}}')
        time_7=$(cat {input[10]} | awk '{{print $6}}')
        time_8=$(cat {input[11]} | awk '{{print $6}}')

        total_time=$(echo "$time_1 + $time_2 + $time_3 + $time_4 + $time_5 + $time_6 + $time_7 + $time_8" | bc)

        # Find peak memory
        peak_mem=$(cat {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]} | sort -k10 -n | tail -n1 | awk '{{print $10}}')

        # Find index size
        component_1=$(ls -l {input[12]} | awk '{{print $5}}')
        component_2=$(ls -l {input[13]} | awk '{{print $5}}')
        total_index_size=$(($component_1+$component_2))

        # Print to output file
        echo "approach,TP,FN,FP,TN,accuracy,sensitivity,specificity,totaltime,peakmem,indexsize" >> {output}
        echo "spumoni_promoted_k{wildcards.k}_w{wildcards.w},$TP,$FN,$FP,$TN,$accuracy,$sensitivity,$specificity,$total_time,$peak_mem,$total_index_size" >> {output}
        """

rule analyze_spumoni_dna_results_into_csv_line_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp5_results/spumoni_dna_k{k}_w{w}/null_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp5_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_1/curr_batch.fa.resources",
        "exp5_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_2/curr_batch.fa.resources",
        "exp5_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_3/curr_batch.fa.resources",
        "exp5_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_4/curr_batch.fa.resources",
        "exp5_results/spumoni_dna_k{k}_w{w}/null_reads/batch_1/curr_batch.fa.resources",
        "exp5_results/spumoni_dna_k{k}_w{w}/null_reads/batch_2/curr_batch.fa.resources",
        "exp5_results/spumoni_dna_k{k}_w{w}/null_reads/batch_3/curr_batch.fa.resources",
        "exp5_results/spumoni_dna_k{k}_w{w}/null_reads/batch_4/curr_batch.fa.resources",
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.thrbv.spumoni",
        "exp5_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.pmlnulldb"
    output:
        "exp5_results/spumoni_dna_k{k}_w{w}/analysis.csv"
    shell:
        """
        # Compute the confusion matrix stats, and accuracy
        num_pos_reads=$(grep -c ">" {input[0]})
        num_null_reads=$(grep -c ">" {input[1]})

        FN=$(wc -l {input[2]} | awk '{{print $1}}')
        TN=$(wc -l {input[3]} | awk '{{print $1}}')
        TP=$(($num_pos_reads-$FN))
        FP=$(($num_null_reads-$TN))

        accuracy=$(echo "scale=4; ($TP+$TN)/($TP+$FN+$FP+$TN)" | bc)
        sensitivity=$(echo "scale=4; ($TP)/($TP+$FN)" | bc)
        specificity=$(echo "scale=4; ($TN)/($FP+$TN)" | bc)

        # Compute time used in each stage
        time_1=$(cat {input[4]} | awk '{{print $6}}')
        time_2=$(cat {input[5]} | awk '{{print $6}}')
        time_3=$(cat {input[6]} | awk '{{print $6}}')
        time_4=$(cat {input[7]} | awk '{{print $6}}')
        time_5=$(cat {input[8]} | awk '{{print $6}}')
        time_6=$(cat {input[9]} | awk '{{print $6}}')
        time_7=$(cat {input[10]} | awk '{{print $6}}')
        time_8=$(cat {input[11]} | awk '{{print $6}}')

        total_time=$(echo "$time_1 + $time_2 + $time_3 + $time_4 + $time_5 + $time_6 + $time_7 + $time_8" | bc)

        # Find peak memory
        peak_mem=$(cat {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]} | sort -k10 -n | tail -n1 | awk '{{print $10}}')

        # Find index size
        component_1=$(ls -l {input[12]} | awk '{{print $5}}')
        component_2=$(ls -l {input[13]} | awk '{{print $5}}')
        total_index_size=$(($component_1+$component_2))

        # Print to output file
        echo "approach,TP,FN,FP,TN,accuracy,sensitivity,specificity,totaltime,peakmem,indexsize" >> {output}
        echo "spumoni_dna_k{wildcards.k}_w{wildcards.w},$TP,$FN,$FP,$TN,$accuracy,$sensitivity,$specificity,$total_time,$peak_mem,$total_index_size" >> {output}
        """

rule analyze_minimap2_results_into_csv_line_exp5:
    input:
        "exp5_read_sets/pos_reads.fa",
        "exp5_read_sets/null_reads.fa",
        "exp5_results/minimap2/pos_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp5_results/minimap2/null_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp5_results/minimap2/pos_reads/batch_1/curr_batch.resources",
        "exp5_results/minimap2/pos_reads/batch_2/curr_batch.resources",
        "exp5_results/minimap2/pos_reads/batch_3/curr_batch.resources",
        "exp5_results/minimap2/pos_reads/batch_4/curr_batch.resources",
        "exp5_results/minimap2/null_reads/batch_1/curr_batch.resources",
        "exp5_results/minimap2/null_reads/batch_2/curr_batch.resources",
        "exp5_results/minimap2/null_reads/batch_3/curr_batch.resources",
        "exp5_results/minimap2/null_reads/batch_4/curr_batch.resources",
        "exp5_indexes/minimap2_index/full_ref.mmi"
    output:
        "exp5_results/minimap2/analysis.csv"
    shell:
        """
        # Compute the confusion matrix stats, and accuracy
        num_pos_reads=$(grep -c ">" {input[0]})
        num_null_reads=$(grep -c ">" {input[1]})

        FN=$(wc -l {input[2]} | awk '{{print $1}}')
        TN=$(wc -l {input[3]} | awk '{{print $1}}')
        TP=$(($num_pos_reads-$FN))
        FP=$(($num_null_reads-$TN))

        accuracy=$(echo "scale=4; ($TP+$TN)/($TP+$FN+$FP+$TN)" | bc)
        sensitivity=$(echo "scale=4; ($TP)/($TP+$FN)" | bc)
        specificity=$(echo "scale=4; ($TN)/($FP+$TN)" | bc)

        # Compute time used in each stage
        time_1=$(cat {input[4]} | awk '{{print $6}}')
        time_2=$(cat {input[5]} | awk '{{print $6}}')
        time_3=$(cat {input[6]} | awk '{{print $6}}')
        time_4=$(cat {input[7]} | awk '{{print $6}}')
        time_5=$(cat {input[8]} | awk '{{print $6}}')
        time_6=$(cat {input[9]} | awk '{{print $6}}')
        time_7=$(cat {input[10]} | awk '{{print $6}}')
        time_8=$(cat {input[11]} | awk '{{print $6}}')

        total_time=$(echo "$time_1 + $time_2 + $time_3 + $time_4 + $time_5 + $time_6 + $time_7 + $time_8" | bc)

        # Find peak memory
        peak_mem=$(cat {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]} | sort -k10 -n | tail -n1 | awk '{{print $10}}')

        # Find index size
        total_index_size=$(ls -l {input[12]} | awk '{{print $5}}')

        # Print to output file
        echo "approach,TP,FN,FP,TN,accuracy,sensitivity,specificity,totaltime,peakmem,indexsize" >> {output}
        echo "minimap2,$TP,$FN,$FP,$TN,$accuracy,$sensitivity,$specificity,$total_time,$peak_mem,$total_index_size" >> {output}
        """

# Section 2.19: Collect all the individual analysis files into one csv

rule collect_all_analysis_files_exp5:
    input:
        expand("exp5_results/spumoni_{type}_k4_w{w}/analysis.csv", type=['promoted','dna'], w=[10, 11]),
        "exp5_results/minimap2/analysis.csv"
    output:
        "exp5_final_output/exp5_analysis.csv"
    shell:
        """
        echo "approach,TP,FN,FP,TN,accuracy,sensitivity,specificity,totaltime,peakmem,indexsize" >> {output}
        for file in {input}; do
            tail -n1 $file >> {output}
        done
        """

