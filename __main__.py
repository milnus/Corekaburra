# General function
from parse_gene_presence_absence import read_gene_presence_absence
from gff_parser import segment_genome_content
from merge_dicts import merge_dicts_counts, merge_dicts_lists
from output_writer_functions import master_info_writer
from time_calculator import time_calculator
from commandline_interface import get_commandline_arguments
import concurrent.futures
import time
import csv
from numpy import mean, std

import sys

# For plots
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    # TODO - Get arguments
    #args = get_commandline_arguments(sys.argv[1:])

    # TODO Check if all gff files are present in input folder

    local_pres_abs = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/Pangenome/gene_presence_absence_roary.csv"
    # TODO - Open gene_presence_absence file and return dict with a key for each core gene cluster and all locus_tags as the value for each key.
    time_start = time.time()
    core_dict, low_freq_dict = read_gene_presence_absence(local_pres_abs,
                                                          0.99, 0.05)
    time_calculator(time_start, time.time(), "reading in gene presence/absence file")

    gff_files = ["/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/Recombination_detection_181020/all_aligned_to_single_reference/GCA_900475985.gff"]
    # gff_files = ["/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/annotations_gff/GCA_004135875.gff"]
    # gff_files = args.input_gffs

    # TODO for-loop over each gff - Try to multiprocess
    # TODO Parse gff and extract core and low frequency genes from gffs
    # TODO - Multi processing

    time_start = time.time()
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

    time_calculator(time_start, time.time(), "Searching gff files for core genomes")

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
    # TODO write comprehensice output format for everything from master_info except low frequency gene present (But take counts)

    # TODO write comprehensive output of low frequency genes present in core regions
    # TODO give argument of output folder!
    # Write master information to output file
    time_start = time.time()
    master_info_writer(master_info_dict, verbose=True)
    time_calculator(time_start, time.time(), "write comprehensive output")

    #####################

    ### WRITE PLOTS ###

    ###################


if __name__ == "__main__":
    main()
