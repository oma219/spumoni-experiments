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
        "exp7_intermediate/step_1/null_reads_2000bp.fa",
        "exp7_read_sets/null_reads.fa"
    shell:
        """
        seqtk seq -A {input} > {output[0]}
        seqtk seq -L 2000 {output[0]} > {output[1]}
        num_lines=$(({num_pos_reads_exp7}*2))
        head -n $num_lines {output[1]} > {output[2]}.before_filter
        python3 {repo_dir}/src/remove_bad_reads.py -i {output[2]}.before_filter > {output[2]}
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
        pbsim --depth 2.0 --prefix exp7_intermediate/step_2/pos_reads --hmm_model {pbsim_model} --accuracy-mean {long_read_acc_exp7} --length-min 2000 {input}

        # remove reads simulated from chrM
        rm exp7_intermediate/step_2/pos_reads_0024.fastq

        cat exp7_intermediate/step_2/pos_reads_*.fastq > {output[0]}
        rm exp7_intermediate/step_2/pos_reads_*.fastq
        rm exp7_intermediate/step_2/pos_reads_*.maf
        rm exp7_intermediate/step_2/pos_reads_*.ref

        seqtk seq -a {output[0]} > {output[1]}

        num_lines=$(({num_null_reads_exp7} * 2))
        head -n $num_lines {output[1]} > {output[2]}.before_filter
        python3 {repo_dir}/src/remove_bad_reads.py -i {output[2]}.before_filter > {output[2]}
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

rule build_spumoni_full_index_exp7:
    input:
        "exp7_full_ref/human_database.fa"
    output:
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa",
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa.thrbv.spumoni",
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa.thrbv.ms",
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa.pmlnulldb",
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa.msnulldb"
    shell:
        """
        # Copy over reference file to specific folder, and build the index
        curr_ref_file="exp7_indexes/spumoni_full_index/full_ref.fa"
        log_file="exp7_indexes/spumoni_full_index/full_ref.fa.log"
        cp {input[0]} $curr_ref_file
        spumoni build -r $curr_ref_file -M -P -k -n &> $log_file
        """

rule build_spumoni_promoted_index_exp7:
    input:
        "exp7_full_ref/human_database.fa"
    output:
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.thrbv.ms",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.pmlnulldb",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.msnulldb"
    shell:
        """
        # Copy over reference file to specific folder, and build the index
        curr_ref_file="exp7_indexes/spumoni_promoted_k{wildcards.k}_w{wildcards.w}/full_ref.fa"
        log_file="exp7_indexes/spumoni_promoted_k{wildcards.k}_w{wildcards.w}/full_ref.fa.log"
        cp {input[0]} $curr_ref_file
        spumoni build -r $curr_ref_file -M -P -k -m -K {wildcards.k} -W {wildcards.w} &> $log_file
        """

rule build_spumoni_dna_index_exp7:
    input:
        "exp7_full_ref/human_database.fa"
    output:
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.thrbv.spumoni",
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.thrbv.ms",
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.pmlnulldb",
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.msnulldb"
    shell:
        """
        # Copy over reference file to specific folder, and build the index
        curr_ref_file="exp7_indexes/spumoni_dna_k{wildcards.k}_w{wildcards.w}/full_ref.fa"
        log_file="exp7_indexes/spumoni_dna_k{wildcards.k}_w{wildcards.w}/full_ref.fa.log"
        cp {input[0]} $curr_ref_file
        spumoni build -r $curr_ref_file -M -P -k -t -K {wildcards.k} -W {wildcards.w} &> $log_file
        """

rule build_minimap2_index_exp7:
    input:
        "exp7_full_ref/human_database.fa"
    output:
        "exp7_indexes/minimap2_index/full_ref.mmi"
    shell:
        """
        # Copy over reference file to specific folder, and build the index
        curr_ref_file="exp7_indexes/minimap2_index/full_ref.fa"
        log_file="exp7_indexes/minimap2_index/full_ref.fa.log"
        cp {input[0]} $curr_ref_file 
        minimap2 -x map-ont -d {output} --split-prefix="exp7_indexes/minimap2_index/split_" {input} &> $log_file
        """

# Section 2.3: Extract a batch of data that will be processed by SPUMONI or minimap2.
#              For SPUMONI, it will just be the current batch while for minimap2 it will
#              a consistently longer prefix of the read.

rule extract_first_batch_of_data_for_spumoni_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa"
    output:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa \
        -n 1 -s 180 -a --spumoni > {output}
        """

rule extract_first_batch_of_data_for_spumoni_full_index_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_1/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa \
        -n 1 -s 180 -a --spumoni > {output}
        """

rule extract_first_batch_of_data_for_minimap2_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa"
    output:
        "exp7_results/minimap2/{class}_reads/batch_1/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa \
        -n 1 -s 180 -a --alignment > {output}
        """

# Section 2.4: Classify the first batch of data using SPUMONI and minimap2

rule classify_first_batch_using_spumoni_promoted_exp7:
    input:
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa"
    output:
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.report",
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -m -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_first_batch_using_spumoni_dna_exp7:
    input:
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa"
    output:
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.report",
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -a -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_first_batch_using_spumoni_full_index_exp7:
    input:
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa",
        "exp7_results/spumoni_full_index/{class}_reads/batch_1/curr_batch.fa"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_1/curr_batch.fa.report",
        "exp7_results/spumoni_full_index/{class}_reads/batch_1/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -n -c
        """

rule classify_first_batch_using_minimap2_exp7:
    input:
        "exp7_indexes/minimap2_index/full_ref.mmi",
        "exp7_results/minimap2/{class}_reads/batch_1/curr_batch.fa"
    output:
        "exp7_results/minimap2/{class}_reads/batch_1/curr_batch.sam",
        "exp7_results/minimap2/{class}_reads/batch_1/curr_batch.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} minimap2 -t 0 -a {input[0]} --split-prefix="exp7_indexes/minimap2_index/split_" {input[1]} > {output[0]}
        """

# Section 2.5: Parse out the reference names in all human genomes put
#              in the index

rule parse_ref_names_from_full_dataset_exp7:
    input:
        list_of_human_genomes_to_index_exp7
    output:
        "exp7_index_ref_name_lists/class_all.txt"
    shell:
        """
        for file in {input}; do
            grep '^>' $file | awk '{{print substr($1,2)}}' >> "exp7_index_ref_name_lists/class_all.txt"
        done
        """

# Section 2.6: Determine which reads each method did not find in the database in order
#              to keep them to classify in next batch

rule determine_reads_to_keep_for_second_batch_spumoni_exp7:
    input:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.report"
    output:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_second_batch_spumoni_full_index_exp7:
    input:
        "exp7_results/spumoni_full_index/{class}_reads/batch_1/curr_batch.fa.report"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_second_batch_minimap2_exp7:
    input:
        "exp7_results/minimap2/{class}_reads/batch_1/curr_batch.sam",
        "exp7_index_ref_name_lists/class_all.txt"
    output:
        "exp7_results/minimap2/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_minimap2_sam.py -s {input[0]} -o {output} -r exp7_index_ref_name_lists/class_all.txt
        """

# Section 2.7: Extract the second batch of data for minimap2 and SPUMONI. This
#              time we will take into account the reads that need to be excluded
#              since they have already been classified.

rule extract_second_batch_of_data_for_spumoni_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa \
        -n 2 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_second_batch_of_data_for_spumoni_full_index_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_full_index/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_2/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa \
        -n 2 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_second_batch_of_data_for_minimap2_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/minimap2/{class}_reads/batch_1/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/minimap2/{class}_reads/batch_2/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa \
        -n 2 -s 180 -k {input[2]} --alignment > {output}
        """

# Section 2.8: Classify the second batch of data given using SPUMONI and minimap2

rule classify_second_batch_using_spumoni_promoted_exp7:
    input:
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa"
    output:
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.report",
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -m -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_second_batch_using_spumoni_dna_exp7:
    input:
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa"
    output:
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.report",
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -a -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_second_batch_using_spumoni_full_index_exp7:
    input:
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa",
        "exp7_results/spumoni_full_index/{class}_reads/batch_2/curr_batch.fa"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_2/curr_batch.fa.report",
        "exp7_results/spumoni_full_index/{class}_reads/batch_2/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -n -c
        """

rule classify_second_batch_using_minimap2_exp7:
    input:
        "exp7_indexes/minimap2_index/full_ref.mmi",
        "exp7_results/minimap2/{class}_reads/batch_2/curr_batch.fa"
    output:
        "exp7_results/minimap2/{class}_reads/batch_2/curr_batch.sam",
        "exp7_results/minimap2/{class}_reads/batch_2/curr_batch.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} minimap2 -t 0 -a {input[0]} --split-prefix="exp7_indexes/minimap2_index/split_" {input[1]} > {output[0]}
        """

# Section 2.9: Determine which reads each method did not find yet in order to 
#              keep them for classifying the next batch

rule determine_reads_to_keep_for_third_batch_spumoni_exp7:
    input:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.report"
    output:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_third_batch_spumoni_full_index_exp7:
    input:
        "exp7_results/spumoni_full_index/{class}_reads/batch_2/curr_batch.fa.report"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_third_batch_minimap2_exp7:
    input:
        "exp7_results/minimap2/{class}_reads/batch_2/curr_batch.sam",
        "exp7_index_ref_name_lists/class_all.txt"
    output:
        "exp7_results/minimap2/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_minimap2_sam.py -s {input[0]} -o {output} -r exp7_index_ref_name_lists/class_all.txt
        """


# Section 2.10: Extract the third batch of data for minimap2 and SPUMONI. This
#               time we will take into account the reads that need to be excluded
#               since they have already been classified.

rule extract_third_batch_of_data_for_spumoni_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa\
        -n 3 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_third_batch_of_data_for_spumoni_full_index_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_full_index/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_3/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa\
        -n 3 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_third_batch_of_data_for_minimap2_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/minimap2/{class}_reads/batch_2/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/minimap2/{class}_reads/batch_3/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa \
        -n 3 -s 180 -k {input[2]} --alignment > {output}
        """

# Section 2.11: Classify the third batch of data given using SPUMONI and minimap2

rule classify_third_batch_using_spumoni_promoted_exp7:
    input:
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa"
    output:
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.report",
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -m -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_third_batch_using_spumoni_dna_exp7:
    input:
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa"
    output:
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.report",
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -a -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_third_batch_using_spumoni_full_index_exp7:
    input:
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa",
        "exp7_results/spumoni_full_index/{class}_reads/batch_3/curr_batch.fa"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_3/curr_batch.fa.report",
        "exp7_results/spumoni_full_index/{class}_reads/batch_3/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -n -c
        """

rule classify_third_batch_using_minimap2_exp7:
    input:
        "exp7_indexes/minimap2_index/full_ref.mmi",
        "exp7_results/minimap2/{class}_reads/batch_3/curr_batch.fa"
    output:
        "exp7_results/minimap2/{class}_reads/batch_3/curr_batch.sam",
        "exp7_results/minimap2/{class}_reads/batch_3/curr_batch.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} minimap2 -t 0 -a {input[0]} --split-prefix="exp7_indexes/minimap2_index/split_" {input[1]} > {output[0]}
        """

# Section 2.12: Determine which reads each method did not find yet in order to 
#               keep those reads for the next batch.

rule determine_reads_to_keep_for_fourth_batch_spumoni_exp7:
    input:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.report"
    output:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_fourth_batch_spumoni_full_index_exp7:
    input:
        "exp7_results/spumoni_full_index/{class}_reads/batch_3/curr_batch.fa.report"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_fourth_batch_minimap2_exp7:
    input:
        "exp7_results/minimap2/{class}_reads/batch_3/curr_batch.sam",
        "exp7_index_ref_name_lists/class_all.txt"
    output:
        "exp7_results/minimap2/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_minimap2_sam.py -s {input[0]} -o {output} -r exp7_index_ref_name_lists/class_all.txt
        """

# Section 2.13: Extract the fouth batch of data for minimap2 and SPUMONI. This
#               time we will take into account the reads that need to be excluded
#               since they have already been classified.

rule extract_fourth_batch_of_data_for_spumoni_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa\
        -n 4 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_fourth_batch_of_data_for_spumoni_full_index_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_full_index/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_4/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa\
        -n 4 -s 180 -k {input[2]}  --spumoni > {output}
        """

rule extract_fourth_batch_of_data_for_minimap2_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/minimap2/{class}_reads/batch_3/curr_batch.fa.reads_to_keep"
    output:
        "exp7_results/minimap2/{class}_reads/batch_4/curr_batch.fa"
    shell:
        """
        python3 {repo_dir}/src/extract_batch.py -i exp7_read_sets/{wildcards.class}_reads.fa \
        -n 4 -s 180 -k {input[2]} --alignment > {output}
        """

# Section 2.14: Classify the fourth batch of data given using SPUMONI and minimap2

rule classify_fourth_batch_using_spumoni_promoted_exp7:
    input:
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin",
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa"
    output:
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.report",
        "exp7_results/spumoni_promoted_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -m -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_fourth_batch_using_spumoni_dna_exp7:
    input:
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa",
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa"
    output:
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.report",
        "exp7_results/spumoni_dna_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -a -K {wildcards.k} -W {wildcards.w} -c
        """

rule classify_fourth_batch_using_spumoni_full_index_exp7:
    input:
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa",
        "exp7_results/spumoni_full_index/{class}_reads/batch_4/curr_batch.fa"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_4/curr_batch.fa.report",
        "exp7_results/spumoni_full_index/{class}_reads/batch_4/curr_batch.fa.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} spumoni run -r {input[0]} -p {input[1]} {output_type_exp7} -n -c
        """

rule classify_fourth_batch_using_minimap2_exp7:
    input:
        "exp7_indexes/minimap2_index/full_ref.mmi",
        "exp7_results/minimap2/{class}_reads/batch_4/curr_batch.fa"
    output:
        "exp7_results/minimap2/{class}_reads/batch_4/curr_batch.sam",
        "exp7_results/minimap2/{class}_reads/batch_4/curr_batch.resources"
    shell:
        """
        {time_prog} {time_format} --output={output[1]} minimap2 -t 0 -a {input[0]} --split-prefix="exp7_indexes/minimap2_index/split_" {input[1]} > {output[0]}
        """

# Section 2.15: Determine which reads each method did not find yet in order to 
#               keep those reads for the next batch.

rule determine_reads_to_keep_for_fifth_batch_spumoni_exp7:
    input:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.report"
    output:
        "exp7_results/spumoni_{type}_k{k}_w{w}/{class}_reads/batch_4/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_fifth_batch_spumoni_full_index_exp7:
    input:
        "exp7_results/spumoni_full_index/{class}_reads/batch_4/curr_batch.fa.report"
    output:
        "exp7_results/spumoni_full_index/{class}_reads/batch_4/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_spumoni_report.py -i {input[0]} -l {input[0]} -o {output}
        """

rule determine_reads_to_keep_for_fifth_batch_minimap2_exp7:
    input:
        "exp7_results/minimap2/{class}_reads/batch_4/curr_batch.sam",
        "exp7_index_ref_name_lists/class_all.txt"
    output:
        "exp7_results/minimap2/{class}_reads/batch_4/curr_batch.fa.reads_to_keep"
    shell:
        """
        python3 {repo_dir}/src/process_minimap2_sam.py -s {input[0]} -o {output} -r exp7_index_ref_name_lists/class_all.txt
        """

# Section 2.16: Analyze each classification approach individually in terms of 
#               TP, FN, FP, TN, time, peak memory, and index size

rule analyze_spumoni_promoted_results_into_csv_line_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp7_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp7_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_1/curr_batch.fa.resources",
        "exp7_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_2/curr_batch.fa.resources",
        "exp7_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_3/curr_batch.fa.resources",
        "exp7_results/spumoni_promoted_k{k}_w{w}/pos_reads/batch_4/curr_batch.fa.resources",
        "exp7_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_1/curr_batch.fa.resources",
        "exp7_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_2/curr_batch.fa.resources",
        "exp7_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_3/curr_batch.fa.resources",
        "exp7_results/spumoni_promoted_k{k}_w{w}/null_reads/batch_4/curr_batch.fa.resources",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.pmlnulldb",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.thrbv.ms",
        "exp7_indexes/spumoni_promoted_k{k}_w{w}/spumoni_full_ref.bin.msnulldb"
    output:
        "exp7_results/spumoni_promoted_k{k}_w{w}/analysis.csv"
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

        # Change index size measurements when using MS-index (default is PML)
        if [[ "{output_type_exp5}" == *"-M"* ]]; then
            component_1=$(ls -l {input[14]} | awk '{{print $5}}')
            component_2=$(ls -l {input[15]} | awk '{{print $5}}')
            total_index_size=$(($component_1+$component_2))
        fi

        # Print to output file
        echo "approach,TP,FN,FP,TN,accuracy,sensitivity,specificity,totaltime,peakmem,indexsize" >> {output}
        echo "spumoni_promoted_k{wildcards.k}_w{wildcards.w},$TP,$FN,$FP,$TN,$accuracy,$sensitivity,$specificity,$total_time,$peak_mem,$total_index_size" >> {output}
        """

rule analyze_spumoni_dna_results_into_csv_line_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp7_results/spumoni_dna_k{k}_w{w}/null_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp7_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_1/curr_batch.fa.resources",
        "exp7_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_2/curr_batch.fa.resources",
        "exp7_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_3/curr_batch.fa.resources",
        "exp7_results/spumoni_dna_k{k}_w{w}/pos_reads/batch_4/curr_batch.fa.resources",
        "exp7_results/spumoni_dna_k{k}_w{w}/null_reads/batch_1/curr_batch.fa.resources",
        "exp7_results/spumoni_dna_k{k}_w{w}/null_reads/batch_2/curr_batch.fa.resources",
        "exp7_results/spumoni_dna_k{k}_w{w}/null_reads/batch_3/curr_batch.fa.resources",
        "exp7_results/spumoni_dna_k{k}_w{w}/null_reads/batch_4/curr_batch.fa.resources",
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.thrbv.spumoni",
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.pmlnulldb",
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.thrbv.ms",
        "exp7_indexes/spumoni_dna_k{k}_w{w}/spumoni_full_ref.fa.msnulldb"
    output:
        "exp7_results/spumoni_dna_k{k}_w{w}/analysis.csv"
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

        # Change index size measurements when using MS-index (default is PML)
        if [[ "{output_type_exp5}" == *"-M"* ]]; then
            component_1=$(ls -l {input[14]} | awk '{{print $5}}')
            component_2=$(ls -l {input[15]} | awk '{{print $5}}')
            total_index_size=$(($component_1+$component_2))
        fi

        # Print to output file
        echo "approach,TP,FN,FP,TN,accuracy,sensitivity,specificity,totaltime,peakmem,indexsize" >> {output}
        echo "spumoni_dna_k{wildcards.k}_w{wildcards.w},$TP,$FN,$FP,$TN,$accuracy,$sensitivity,$specificity,$total_time,$peak_mem,$total_index_size" >> {output}
        """

rule analyze_spumoni_full_index_results_into_csv_line_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/spumoni_full_index/pos_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp7_results/spumoni_full_index/null_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp7_results/spumoni_full_index/pos_reads/batch_1/curr_batch.fa.resources",
        "exp7_results/spumoni_full_index/pos_reads/batch_2/curr_batch.fa.resources",
        "exp7_results/spumoni_full_index/pos_reads/batch_3/curr_batch.fa.resources",
        "exp7_results/spumoni_full_index/pos_reads/batch_4/curr_batch.fa.resources",
        "exp7_results/spumoni_full_index/null_reads/batch_1/curr_batch.fa.resources",
        "exp7_results/spumoni_full_index/null_reads/batch_2/curr_batch.fa.resources",
        "exp7_results/spumoni_full_index/null_reads/batch_3/curr_batch.fa.resources",
        "exp7_results/spumoni_full_index/null_reads/batch_4/curr_batch.fa.resources",
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa.thrbv.spumoni",
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa.pmlnulldb",
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa.thrbv.ms",
        "exp7_indexes/spumoni_full_index/spumoni_full_ref.fa.msnulldb"
    output:
        "exp7_results/spumoni_full_index/analysis.csv"
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

        # Change index size measurements when using MS-index (default is PML)
        if [[ "{output_type_exp5}" == *"-M"* ]]; then
            component_1=$(ls -l {input[14]} | awk '{{print $5}}')
            component_2=$(ls -l {input[15]} | awk '{{print $5}}')
            total_index_size=$(($component_1+$component_2))
        fi

        # Print to output file
        echo "approach,TP,FN,FP,TN,accuracy,sensitivity,specificity,totaltime,peakmem,indexsize" >> {output}
        echo "spumoni_full_index,$TP,$FN,$FP,$TN,$accuracy,$sensitivity,$specificity,$total_time,$peak_mem,$total_index_size" >> {output}
        """

rule analyze_minimap2_results_into_csv_line_exp7:
    input:
        "exp7_read_sets/pos_reads.fa",
        "exp7_read_sets/null_reads.fa",
        "exp7_results/minimap2/pos_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp7_results/minimap2/null_reads/batch_4/curr_batch.fa.reads_to_keep",
        "exp7_results/minimap2/pos_reads/batch_1/curr_batch.resources",
        "exp7_results/minimap2/pos_reads/batch_2/curr_batch.resources",
        "exp7_results/minimap2/pos_reads/batch_3/curr_batch.resources",
        "exp7_results/minimap2/pos_reads/batch_4/curr_batch.resources",
        "exp7_results/minimap2/null_reads/batch_1/curr_batch.resources",
        "exp7_results/minimap2/null_reads/batch_2/curr_batch.resources",
        "exp7_results/minimap2/null_reads/batch_3/curr_batch.resources",
        "exp7_results/minimap2/null_reads/batch_4/curr_batch.resources",
        "exp7_indexes/minimap2_index/full_ref.mmi"
    output:
        "exp7_results/minimap2/analysis.csv"
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

# Section 2.17: Collect all the individual analysis files into one csv

rule collect_all_analysis_files_exp7:
    input:
        expand("exp7_results/spumoni_{type}_k4_w{w}/analysis.csv", type=['promoted','dna'], w=[10, 11]),
        "exp7_results/minimap2/analysis.csv",
        "exp7_results/spumoni_full_index/analysis.csv"
    output:
        "exp7_final_output/exp7_analysis.csv"
    shell:
        """
        echo "approach,TP,FN,FP,TN,accuracy,sensitivity,specificity,totaltime,peakmem,indexsize" >> {output}
        for file in {input}; do
            tail -n1 $file >> {output}
        done
        """


