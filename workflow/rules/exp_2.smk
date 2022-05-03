##################################################
# Name: exp_2.smk
# Description: Contains the workflow and methods
#              needed for experiment 2.
#
# Date: May 3, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_genomes_in_dataset(wildcards):
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

####################################################
# Section 2: Rules needed for this experiment type
####################################################

rule produce_list_of_genomes_exp2:
    input:
        get_genomes_in_dataset
    output:
        "file_list/dataset_file_list.txt"
    run:
        with open(output[0], "w") as fd:
            for file_name in input:
                fd.write(file_name + "\n")

rule build_index_promotion_minimizer_exp2:
    input:
        "file_list/dataset_file_list.txt"
    output:
        "current_index_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "current_index_k{k}_w{w}/spumoni_full_ref.bin.thrbv.ms"
    run:
        shell("spumoni build -i {input[0]} -b current_index_k{wildcards.k}_w{wildcards.w}/ -M -P -m -K {wildcards.k} -W {wildcards.w} &> current_index_k{wildcards.k}_w{wildcards.w}/build.log")

rule gather_index_stats_promotion_minimizer_exp2:
    input:
        "current_index_k{k}_w{w}/spumoni_full_ref.bin.thrbv.spumoni",
        "current_index_k{k}_w{w}/spumoni_full_ref.bin.thrbv.ms"
    output:
        "results/promotion/index_results_k{k}_w{w}.txt"
    run:
        shell("ls -lh {input[0]} > {output[0]}")
        #shell("rm -r current_index_k{wildcards.k}_w{wildcards.w}/")



