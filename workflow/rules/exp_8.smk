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
    our SPUMONI index to use for contamination detection.
    """
    input_files = []
    for num in range(1, num_datasets_exp8):
        for data_file in os.listdir(f"data/dataset_{num}"):
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

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Build the SPUMONI index over genomes that are common sources of contamination

rule build_full_reference_file_exp8:
    input:
        list_of_all_genomes_to_index_exp8
    output:
        "exp8_full_reference/full_ref.fa"
    shell:
        """
        cat {input} > {output}
        """

rule build_spumoni_index_of_contamination_exp8:
    input:
        "exp8_full_reference/full_ref.fa"
    output:
        "exp8_index/spumoni_full_ref.fa",
        "exp8_index/spumoni_full_ref.fa.thrbv.ms",
        "exp8_index/spumoni_full_ref.fa.msnulldb"
    shell:
        """
        cp {input} exp8_index/full_ref.fa
        spumoni build -r exp8_index/full_ref.fa -n -M -P &> exp8_index/full_ref.fa.log
        """


# Section 2.2: Compute MS for each assembly against the contamination database

rule run_spumoni_against_database_exp8:
    input:
        get_current_assembly_exp8,
        "exp8_index/spumoni_full_ref.fa"
    output:
        "exp8_results/assembly_{asm_num}/input_assembly.fa",
        "exp8_results/assembly_{asm_num}/input_assembly.fa.lengths"
    shell:
        """
        cp {input[0]} {output[0]}
        spumoni run -t 16 -r {input[1]} -p {output[0]} -M -n
        """

# Section 2.3: Generate analysis/csv files to use for plotting

rule analyze_spumoni_results_against_database_exp8:
    input:
        "exp8_results/assembly_{asm_num}/input_assembly.fa.lengths"
    output:
        "exp8_analysis/assembly_{asm_num}/assembly_length_report.txt"
    shell:
        """
        python3 {repo_dir}/src/process_lengths_exp8.py -i {input} \
        -o exp8_analysis/assembly_{wildcards.asm_num}/
        """

