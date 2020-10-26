# General function
from parse_gene_presence_absence import read_gene_presence_absence
from gff_parser import segment_genome_content
from merge_dicts import merge_dicts_counts, merge_dicts_lists
import concurrent.futures
from numpy import mean, std
import sys

# For plots
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#def get_arguments():
    # TODO - Set up arguments parser
    # ARGUMENTS:
    #


def main():
    # TODO - Get arguments

    # TODO Check if all gff files are present in input folder

    # TODO - Open gene_presence_absence file and return dict with a key for each core gene cluster and all locus_tags as the value for each key.
    core_dict, low_freq_dict = read_gene_presence_absence("/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/Pangenome/gene_presence_absence_roary.csv",
                                                          0.99, 0.05)

    gff_files = ["/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/Recombination_detection_181020/all_aligned_to_single_reference/GCA_900475985.gff",
                 "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/annotations_gff/GCA_004135875.gff"]
    # TODO for-loop over each gff - Try to multiprocess
    # TODO Parse gff and extract core and low frequency genes from gffs
    # TODO - Multi processing

    # Initialise dictionaries to contain results from all gff files
    core_neighbour_pairs = {}
    core_neighbour_distance = {}
    core_neighbour_accessory_count = {}
    core_neighbour_low_freq = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        results = [executor.submit(segment_genome_content, gff, core_dict, low_freq_dict) for gff in gff_files]
        # core_pairs, distance, acc_count, low_freq

        for output in concurrent.futures.as_completed(results):
            # Split the outputs from
            core_pairs, distance, acc_count, low_freq, master_info_dict = output.result()

            # Merge results into single dictionaries
            core_neighbour_pairs = merge_dicts_counts(core_neighbour_pairs, core_pairs)
            core_neighbour_distance = merge_dicts_lists(core_neighbour_distance, distance)
            core_neighbour_accessory_count = merge_dicts_counts(core_neighbour_accessory_count, acc_count)
            core_neighbour_low_freq = merge_dicts_lists(core_neighbour_low_freq, low_freq)

    print(core_neighbour_pairs)
    print(core_neighbour_distance)
    print(core_neighbour_accessory_count)
    print(core_neighbour_low_freq)
    print(master_info_dict)
    print(sys.getsizeof(core_neighbour_pairs))
    print(sys.getsizeof(core_neighbour_distance))
    print(sys.getsizeof(core_neighbour_accessory_count))
    print(sys.getsizeof(core_neighbour_low_freq))
    print(sys.getsizeof(master_info_dict))

    ### FUNCTION ###
    # TODO Get the synteny of genes if genome is complete with score 1-n_core_genes
    ################

    ### DETERMINE CONSENSUS ###
    # TODO Determine start cluster from possible consensus from complete genomes - else determine relative consensus from connections
    # TODO Do a greedy walk through a graph where all nodes are a core cluster and the edge weight is the number of times two core clusters are observed to be neighbours
    # TODO for-loop - Determine following clusters from connections to first cluster
    # TODO Determine all alternative connections between core genes and their frequency.
    ###########################

    ### DO CALCULATIONS ###
    # TODO mean number length between core genes
    # for neighbours in core_neighbour_distance:
    #     print(mean(core_neighbour_distance[neighbours]))
    #     print(std(core_neighbour_distance[neighbours]))

    # TODO SD of length between core genes
    # TODO plot -
    # TODO kruwalski plot
    #######################

    ### WRITE OUTPUTS ###
    # TODO write raw distance outputs in long format
    # TODO possibly construct pseudo core with core-core distances

    #####################

    ### WRITE PLOTS ###

    ###################


if __name__ == "__main__":
    main()
