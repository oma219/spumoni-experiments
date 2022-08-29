##################################################
# Name: exp_8.smk
# Description: Contains the workflow and methods
#              needed for experiment 8.
#
# Date: July 20, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def list_of_all_genomes_to_index_exp8(wildcards):
    """
    This method grabs all the files that we need for
    our index to use for contamination detection.
    """
    input_files = []
    for num in range(1, num_datasets_exp8):
        for data_file in os.listdir(f"data/dataset_{num}"):
            if data_file.endswith(".fna"):
                input_files.append(f"{base_dir}/data/dataset_{num}/{data_file}")
    return input_files

def subset_database_for_minimap2_index_exp8(wildcards):
    """
    This method does a similar thing to previous method, however,
    I only want a subset of the genomes from each class since we
    want to visualize it in a dot plot.
    """
    input_files = []
    for num in range(1, num_datasets_exp8):
        for data_file in os.listdir(f"data/dataset_{num}")[:5]:
            if data_file.endswith(".fna"):
                input_files.append(f"{base_dir}/data/dataset_{num}/{data_file}")
    return input_files

def get_current_assembly_exp8(wildcards):
    """
    Returns the path to the current assembly that we want
    to compute MS for against the database
    """
    input_files = []
    for data_file in os.listdir(f"data/assembly_{wildcards.asm_num}"):
        if data_file.endswith(".fna") or data_file.endswith(".fa"):
            input_files.append(f"{base_dir}/data/assembly_{wildcards.asm_num}/{data_file}")
    assert len(input_files) == 1, "assertion error: there should only be one assembly"
    return input_files[0]

def find_report_files_for_each_approach_exp8(wildcards):
    """
    Returns a list of paths for the *report files that were
    outputed by the time prog that include the elapsed time
    """
    input_files = []
    method_list = ["spumoni_1", "spumoni_2", "blast"]
    for method in method_list:
        for query_num in range(1, num_assemblies_exp8+1):
            next_data_file = f"exp8_results/{method}/assembly_{query_num}/input_assembly.fa.report"
            input_files.append(next_data_file)
            assert os.path.exists(next_data_file), "assertion error: a data file path is not valid"
    return input_files

####################################################
# Section 2: Rules needed for this experiment type
####################################################

#   Section 2.1: Build the various indexes that will be used for
#   finding contamination such as SPUMONI 1 (no digestion), 
#   SPUMONI 2 (minimizer-alphabet k=4/w=11), and BLAST+ index

rule build_full_reference_file_exp8:
    input:
        list_of_all_genomes_to_index_exp8
    output:
        "exp8_full_ref/full_ref.fa"
    shell:
        """
        cat {input} > {output}
        """

rule build_spumoni1_index_for_contamination_exp8:
    input:
        "exp8_full_ref/full_ref.fa"
    output:
        "exp8_indexes/spumoni1_index/spumoni_full_ref.fa",
        "exp8_indexes/spumoni1_index/spumoni_full_ref.fa.thrbv.ms",
        "exp8_indexes/spumoni1_index/spumoni_full_ref.fa.msnulldb"
    shell:
        """
        input_file="exp8_indexes/spumoni1_index/full_ref.fa"
        cp {input} $input_file
        spumoni build -r $input_file -n -M -P &> exp8_indexes/spumoni1_index/full_ref.fa.log
        """

rule build_spumoni2_index_for_contamination_exp8:
    input:
        "exp8_full_ref/full_ref.fa"
    output:
        "exp8_indexes/spumoni2_index/spumoni_full_ref.bin",
        "exp8_indexes/spumoni2_index/spumoni_full_ref.bin.thrbv.ms",
        "exp8_indexes/spumoni2_index/spumoni_full_ref.bin.msnulldb",
        "exp8_indexes/spumoni2_index/spumoni_full_ref.bin.thrbv.spumoni",
        "exp8_indexes/spumoni2_index/spumoni_full_ref.bin.pmlnulldb"
    shell:
        """
        input_file="exp8_indexes/spumoni2_index/full_ref.fa"
        cp {input} $input_file
        spumoni build -r $input_file -m -M -P &> exp8_indexes/spumoni2_index/full_ref.fa.log
        """

rule build_blast_index_for_contamination_exp8:
    input:
        "exp8_full_ref/full_ref.fa"
    output:
        "exp8_indexes/blast_index/contaminant_db.ndb"
    shell:
        """
        input_file="exp8_indexes/blast_index/full_ref.fa"
        out_db="exp8_indexes/blast_index/contaminant_db"
        cp {input} $input_file
        makeblastdb -in $input_file -dbtype 'nucl' -parse_seqids -out $out_db
        """

#   Section 2.2: For each method, "query" the assembly against the
#   the contamination databaes. For SPUMONI 1 and 2, compute PMLs w.r.t
#   to the their databases, and then for BLAST compute the 

rule run_spumoni1_against_database_exp8:
    input:
        get_current_assembly_exp8,
        "exp8_indexes/spumoni1_index/spumoni_full_ref.fa"
    output:
        "exp8_results/spumoni_1/assembly_{asm_num}/input_assembly.fa",
        "exp8_results/spumoni_1/assembly_{asm_num}/input_assembly.fa.pseudo_lengths",
        "exp8_results/spumoni_1/assembly_{asm_num}/input_assembly.fa.report"
    shell:
        """
        cp {input[0]} {output[0]}
        {time_prog} {time_format} --output={output[2]} spumoni run -t 16 -r {input[1]} -p {output[0]} -P -n
        """

rule run_spumoni2_against_database_exp8:
    input:
        get_current_assembly_exp8,
        "exp8_indexes/spumoni2_index/spumoni_full_ref.bin"
    output:
        "exp8_results/spumoni_2/assembly_{asm_num}/input_assembly.fa",
        "exp8_results/spumoni_2/assembly_{asm_num}/input_assembly.fa.pseudo_lengths",
        "exp8_results/spumoni_2/assembly_{asm_num}/input_assembly.fa.report"
    shell:
        """
        cp {input[0]} {output[0]}
        {time_prog} {time_format} --output={output[2]} spumoni run -t 16 -r {input[1]} -p {output[0]} -P -m
        """

rule run_blast_against_database_exp8:
    input:
        get_current_assembly_exp8,
        "exp8_indexes/blast_index/contaminant_db.ndb"
    output:
        "exp8_results/blast/assembly_{asm_num}/input_assembly.fa",
        "exp8_results/blast/assembly_{asm_num}/output_results.txt",
        "exp8_results/blast/assembly_{asm_num}/input_assembly.fa.report"
    shell:
        """
        input_db="exp8_indexes/blast_index/contaminant_db"
        cp {input[0]} {output[0]}
        {time_prog} {time_format} --output={output[2]} blastn -db $input_db -query {output[0]} -out {output[1]} -max_target_seqs 10 -num_threads 16 -mt_mode 1 -outfmt 6
        """

#   Section 2.3: Analyze the time taken for querying for each method
#   for each assembly and save it to csv file

rule analyze_query_time_results_exp8:
    input:
        find_report_files_for_each_approach_exp8
    output:
        "exp8_analysis/time/query_time_report.csv"
    run:
        def extract_elapsed_time_exp8(data_path):
            # return the elapsed time based on time program output format
            with open(data_path, "r") as in_fd:
                all_lines = [x.strip() for x in in_fd.readlines()]
                assert len(all_lines) == 1
                return float(all_lines[0].split()[5])

        # Extract the data from the *.report files and write to csv
        with open(output[0], "w") as out_fd:
            header = ",".join(["method"] + [f"assembly_{x}" for x in range(1, num_assemblies_exp8+1)]) 
            out_fd.write(header + "\n")    
            pos = 0
            for method in ["spumoni_1", "spumoni_2", "blast"]:
                out_list = [method]
                for query_num in range(1, num_assemblies_exp8+1):
                    out_list.append(extract_elapsed_time_exp8(input[pos]))
                    pos += 1
                out_fd.write(",".join([str(x) for x in out_list]) + "\n")

#   Section 2.4: Analyze the contamination detection by examining either
#   lengths files for SPUMONI 1/2 or the local alignments file from 
#   BLAST+

rule analyze_spumoni_results_against_database_exp8:
    input:
        "exp8_results/spumoni_{type}/assembly_{asm_num}/input_assembly.fa.pseudo_lengths"
    output:
        "exp8_analysis/detection/spumoni_{type}/assembly_{asm_num}/assembly_length_report.txt"
    shell:
        """
        python3 {repo_dir}/src/process_lengths_exp8.py -i {input} \
        -o exp8_analysis/detection/spumoni_{wildcards.type}/assembly_{wildcards.asm_num}/
        """

rule analyze_blast_results_against_database_exp8:
    input:
        "exp8_results/blast/assembly_{asm_num}/output_results.txt"
    output:
        "exp8_analysis/detection/blast/assembly_{asm_num}/alignment_report.txt"
    shell:
        """
        python3 {repo_dir}/src/process_blast_report_exp8.py -i {input[0]} \
         -o exp8_analysis/detection/blast/assembly_{wildcards.asm_num}/
        """

#   Section 2.5: Generate PAF files with assembly alignments to the 
#   contaminant database in order to create dot plots

rule build_minimap2_index_for_paf_creation_exp8:
    input:
        subset_database_for_minimap2_index_exp8
    output:
        "exp8_dotplot_data/minimap2_index/full_ref.mmi"
    shell:
        """
        ref_file="exp8_dotplot_data/minimap2_index/full_ref.fa"
        cat {input} > $ref_file
        minimap2 -d {output} $ref_file 
        """

rule run_minimap2_for_paf_files_exp8:
    input:
        get_current_assembly_exp8,
        "exp8_dotplot_data/minimap2_index/full_ref.mmi"
    output:
        "exp8_dotplot_data/minimap2_results/assembly_{asm_num}/aln_results.paf"
    shell:
        """
        index="exp8_dotplot_data/minimap2_index/full_ref.fa"
        input_file="exp8_dotplot_data/minimap2_results/assembly_{wildcards.asm_num}/input_assembly.fa"
        cp {input[0]} $input_file
        minimap2 -c -x asm5 $index $input_file > {output}
        """



