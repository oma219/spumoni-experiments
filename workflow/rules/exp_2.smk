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

# Section 2.1: Build file-list to use with SPUMONI

rule produce_list_of_genomes_exp2:
    input:
        get_genomes_in_dataset
    output:
        "file_list/dataset_file_list.txt"
    run:
        with open(output[0], "w") as fd:
            for file_name in input:
                fd.write(file_name + "\n")

# Section 2.2: Builds index using promoted-minimizer alphabet and extracts stats

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
    shell:
        """
        input_file="current_index_k{wildcards.k}_w{wildcards.w}/build.log"
        ms_index_file="current_index_k{wildcards.k}_w{wildcards.w}/spumoni_full_ref.bin.thrbv.ms"
        slp_file="current_index_k{wildcards.k}_w{wildcards.w}/spumoni_full_ref.bin.slp"
        pml_index_file="current_index_k{wildcards.k}_w{wildcards.w}/spumoni_full_ref.bin.thrbv.spumoni"
        type="promotion"
        k={wildcards.k}
        w={wildcards.w}
        r=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $6}}' | sed 's/,//')
        n=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $9}}' | sed 's/,//')
        n_over_r=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $12}}')
        ms_size=$(ls -l $ms_index_file | awk '{{print $5}}')
        slp_size=$(ls -l $slp_file | awk '{{print $5}}')
        pml_size=$(ls -l $pml_index_file | awk '{{print $5}}')
        printf '%s,%s,%d,%d,%s,%s,%s,%s,%s,%s\n' {dataset_name_exp2} $type $k $w $n $r $n_over_r $ms_size $slp_size $pml_size > {output[0]}
        rm -r current_index_k{wildcards.k}_w{wildcards.w}/
        """

# Section 2.3: Builds index using DNA-minimizer alpabet and extracts stats

rule build_index_dna_minimizer_exp2:
    input:
        "file_list/dataset_file_list.txt"
    output:
        "current_index_k{k}_w{w}_dna/spumoni_full_ref.fa.thrbv.spumoni",
        "current_index_k{k}_w{w}_dna/spumoni_full_ref.fa.thrbv.ms"
    run:
        shell("spumoni build -i {input[0]} -b current_index_k{wildcards.k}_w{wildcards.w}_dna/ -M -P -t -K {wildcards.k} -W {wildcards.w} &> current_index_k{wildcards.k}_w{wildcards.w}_dna/build.log")

rule gather_index_stats_dna_minimizer_exp2:
    input:
        "current_index_k{k}_w{w}_dna/spumoni_full_ref.fa.thrbv.spumoni",
        "current_index_k{k}_w{w}_dna/spumoni_full_ref.fa.thrbv.ms"
    output:
        "results/dna/index_results_k{k}_w{w}.txt"
    shell:
        """
        input_file="current_index_k{wildcards.k}_w{wildcards.w}_dna/build.log"
        ms_index_file="current_index_k{wildcards.k}_w{wildcards.w}_dna/spumoni_full_ref.fa.thrbv.ms"
        slp_file="current_index_k{wildcards.k}_w{wildcards.w}_dna/spumoni_full_ref.fa.slp"
        pml_index_file="current_index_k{wildcards.k}_w{wildcards.w}_dna/spumoni_full_ref.fa.thrbv.spumoni"
        type="dna"
        k={wildcards.k}
        w={wildcards.w}
        r=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $6}}' | sed 's/,//')
        n=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $9}}' | sed 's/,//')
        n_over_r=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $12}}')
        ms_size=$(ls -l $ms_index_file | awk '{{print $5}}')
        slp_size=$(ls -l $slp_file | awk '{{print $5}}')
        pml_size=$(ls -l $pml_index_file | awk '{{print $5}}')
        printf '%s,%s,%d,%d,%s,%s,%s,%s,%s,%s\n' {dataset_name_exp2} $type $k $w $n $r $n_over_r $ms_size $slp_size $pml_size > {output[0]}
        rm -r current_index_k{wildcards.k}_w{wildcards.w}_dna/
        """

# Section 2.4: Builds index over total input, and records stats

rule build_total_index_exp2:
    input:
        "file_list/dataset_file_list.txt"
    output:
        "full_index/spumoni_full_ref.fa.thrbv.spumoni",
        "full_index/spumoni_full_ref.fa.thrbv.ms"
    run:
        shell("spumoni build -i {input[0]} -b full_index/ -M -P -n &> full_index/build.log")

rule gather_index_stats_for_full_index_exp2:
    input:
        "full_index/spumoni_full_ref.fa.thrbv.spumoni",
        "full_index/spumoni_full_ref.fa.thrbv.ms"
    output:
        "total_results/exp2_full_index_stats.csv"
    shell:
        """
        input_file="full_index/build.log"
        ms_index_file="full_index/spumoni_full_ref.fa.thrbv.ms"
        slp_file="full_index/spumoni_full_ref.fa.slp"
        pml_index_file="full_index/spumoni_full_ref.fa.thrbv.spumoni"
        r=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $6}}' | sed 's/,//')
        n=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $9}}' | sed 's/,//')
        n_over_r=$(grep 'build_ms' $input_file | grep 'bwt statistics' | awk '{{print $12}}')
        ms_size=$(ls -l $ms_index_file | awk '{{print $5}}')
        slp_size=$(ls -l $slp_file | awk '{{print $5}}')
        pml_size=$(ls -l $pml_index_file | awk '{{print $5}}')
        printf 'name,n,r,n_over_r,ms_size,slp_size,pml_size' > {output[0]}
        printf '%s,%s,%s,%s,%s,%s,%s\n' {dataset_name_exp2} $n $r $n_over_r $ms_size $slp_size $pml_size > {output[0]}
        rm -r full_index/
        """    

# Section 2.5: Gathers all statistics and stores in a single csv file

rule gather_all_index_stats_promotion_exp2:
    input:
        expand("results/{type}/index_results_k{k}_w{w}.txt", type=["promotion", "dna"], k="4", w=large_window_sizes)
    output:
        "total_results/exp2_index_stats.csv"
    run:
        out_fd = open(output[0], "w")
        out_fd.write("name,type,k,w,n,r,n_over_r,ms_size,slp_size,pml_size\n")
        out_fd.close()

        # Put all data into the csv-file
        for file in input:
            shell("cat {file} >> {output[0]}")