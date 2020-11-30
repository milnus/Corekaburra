# General function
from parse_gene_presence_absence import read_gene_presence_absence
from gff_parser import segment_genome_content
from merge_dicts import merge_dicts_counts, merge_dicts_lists, merge_first_genes
from output_writer_functions import master_info_writer, \
    write_consensus_core_gene_synteny, \
    write_alternative_core_gene_counts
from consesus_core_genome import determine_core_gene_consesus, identify_rearrangements, characterise_rearrangements
from time_calculator import time_calculator
from commandline_interface import get_commandline_arguments
import concurrent.futures
import time
import csv
import numpy as np

import sys

# For plots
import numpy as np
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt


def main():
    total_time_start = time.time()
    # TODO - Make commandline input take the Panaroo output folder to easier get access to gene_presence_absence and gene_date file
    args = get_commandline_arguments(sys.argv[1:])

    # TODO Check if all gff files are present in input folder

    # local_pres_abs = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/Pangenome/gene_presence_absence_roary.csv"
    # TODO - Open gene_presence_absence file and return dict with a key for each core gene cluster and all locus_tags as the value for each key.
    time_start = time.time()
    core_dict, low_freq_dict = read_gene_presence_absence(args.input_pres_abs,
                                                          0.99, 0.05)
    time_calculator(time_start, time.time(), "reading in gene presence/absence file")

    # TODO - Add in the refound genes into the gff files and print the corrected GFF files.
    # gff_files = ["/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/Recombination_detection_181020/all_aligned_to_single_reference/GCA_900475985.gff",
    #              "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Accessory_exploration/Micro_evolution/Emm75/annotations_gff/GCA_004135875.gff"]
    gff_files = args.input_gffs

    # Loop over all gffs and extract info from each of them.
    time_start = time.time()
    # Initialise dictionaries to contain results from all gff files
    core_neighbour_pairs = {}
    core_neighbour_distance = {}
    core_neighbour_accessory_count = {}
    core_neighbour_low_freq = {}
    master_info_total = {}
    non_core_contig_info = {}
    merged_start_gene_clusters = []
    merged_second_gene_clusters = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        print(f'{len(gff_files)} GFF files to process')
        results = [executor.submit(segment_genome_content, gff, core_dict, low_freq_dict, i)
                   for i, gff in enumerate(gff_files)]

        for output in concurrent.futures.as_completed(results):
            # Split the outputs
            core_pairs, distance, acc_count, \
            low_freq, master_info_return, \
            core_less_contigs_return, start_gene_cluster = output.result()

            # Merge results into single/master dictionaries
            core_neighbour_pairs = merge_dicts_counts(core_neighbour_pairs, core_pairs)
            core_neighbour_distance = merge_dicts_lists(core_neighbour_distance, distance)
            core_neighbour_accessory_count = merge_dicts_lists(core_neighbour_accessory_count, acc_count)
            core_neighbour_low_freq = merge_dicts_lists(core_neighbour_low_freq, low_freq)
            master_info_total.update(master_info_return)
            non_core_contig_info.update(core_less_contigs_return)
            merged_start_gene_clusters, merged_second_gene_clusters = merge_first_genes(start_gene_cluster,
                                                                                        merged_start_gene_clusters,
                                                                                        merged_second_gene_clusters,
                                                                                        core_pairs[0])

    time_calculator(time_start, time.time(), "searching gff files for core genomes")


    ### FUNCTION ###
    # Determine the most common core gene synteny.
    time_start = time.time()
    # Find the core gene synteny and possible core genes with alternative neighbours
    consensus_core_genome, possible_rearrangement_genes, core_path_coverage = determine_core_gene_consesus(core_neighbour_pairs, merged_start_gene_clusters, merged_second_gene_clusters)

    # Identify alternative connections and their occurrence
    alt_core_pairs, alt_core_pair_count = identify_rearrangements(consensus_core_genome, possible_rearrangement_genes, master_info_total)

    # TODO look at number of alternative core-genome pairs. if one give missing info, if two try to find size, if three or more too complex.
    rearrangement_predictions = characterise_rearrangements(alt_core_pairs)

    time_calculator(time_start, time.time(), "determining best core gene synteny")
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
    print("Printing outputs")
    # TODO Write out coreless contigs and their accessory information
    # TODO write raw distance outputs in long format
    # TODO possibly construct pseudo core with core-core distances
    # TODO write comprehensice output format for everything from master_info except low frequency gene present (But take counts)

    # TODO write comprehensive output of low frequency genes present in core regions
    # TODO give argument of output folder!
    # Write master information to output file
    time_start = time.time()
    master_info_writer(master_info_total, verbose=True)
    # Write outputs related to core gene synteny
    print(consensus_core_genome)
    print(core_path_coverage)
    write_consensus_core_gene_synteny(consensus_core_genome, core_path_coverage)
    write_alternative_core_gene_counts(alt_core_pair_count)
    rearrangement_predictions #TODO - write output for rearrangement prediciton - First possibly make predictions for two alternative core pairs

    time_calculator(time_start, time.time(), "writing output files")
    #####################

    ### WRITE PLOTS ###

    ###################
    time_calculator(total_time_start, time.time(), "running the entire program")

if __name__ == "__main__":
    main()
