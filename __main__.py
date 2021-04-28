# General function
from parse_gene_presence_absence import read_gene_presence_absence
from gff_parser import segment_genome_content
from correct_gffs import correct_gffs
from merge_dicts import merge_dicts_counts, merge_dicts_lists, merge_first_genes
from output_writer_functions import master_info_writer, \
    write_consensus_core_gene_synteny, \
    write_core_gene_coverage, \
    write_alternative_core_gene_counts, \
    write_core_gene_types
from consesus_core_genome import determine_core_gene_consesus, identify_rearrangements, \
    characterise_rearrangements, core_pair_matrix
from check_inputs import define_input_source, check_gene_data, check_gff_files, check_gene_alignments
from time_calculator import time_calculator
from commandline_interface import get_commandline_arguments
from construct_multi_fasta_genome import construct_consensus_alignment
# from plots import consesus_genome_coverage
import concurrent.futures
from os import listdir, mkdir
from os.path import join
import time

import sys

# TODO - POSSIBLE NAME - Corkaburra or Corekaburra

def main():
    total_time_start = time.time()
    # get arguments from the commandline
    # TODO - make the command line print the help if no input is given.
    args = get_commandline_arguments(sys.argv[1:])

    ## Check input
    if not args.quiet:
        print("\n----Checking presence of input files in pan genome folder----\n")

    # Check if Panaroo or Roary input folder is given
    source_program, input_pres_abs_file_path = define_input_source(args.input_pan)

    # Check if gene_data file is present if Panaroo input is given an gffs should be annotated
    if args.annotate and source_program is not 'Rorary':
        gene_data_path = check_gene_data(args.input_pan)

    if not args.quiet:
        print(f"Pan genome determined to come from {source_program}")
        print("All files found, let's move on!\n")
        print("--------------------------------------------------------------\n")

    ## Read in gene presence absence file
    time_start = time.time()
    # TODO - check if the genes merged with a ; are neighbours in given genome. If, then keep the group as nothing is in between.
    # TODO - Add the user specified thresholds for core and low frequency genes.
    core_dict, low_freq_dict, acc_gene_dict, attribute_dict = read_gene_presence_absence(input_pres_abs_file_path,
                                                          1, 0.05, source_program)
    time_calculator(time_start, time.time(), "reading in gene presence/absence file")
    # TODO - Check if all gff files are present in input folder and gene presence absence file - report missing files or mismatches.

    # If source program is Panaroo, see if alignments folder is available, and core genes are contained in it.
    if source_program == 'Panaroo':
        alignment_folder = check_gene_alignments(args.input_pan, core_dict)

    # Check presence of gff files.
    if check_gff_files(args.input_gffs):
        print("All .gff files were found!")

    # Construct output folder
    try:
        mkdir(args.output_path)
        print("Output folder constructed")  # TODO - make verbose
    except FileExistsError:
        print("Output folder exists") # TODO - make verbose

    # Add in the refound genes into the gff files and print the corrected GFF files.
    if source_program == "Panaroo" and args.annotate:
        # TODO - Add a timer and a message about the process starting
        corrected_folder = correct_gffs(args.input_gffs, gene_data_path, args.output_path, attribute_dict)

        args.input_gffs = [join(corrected_folder, file) for file in listdir(corrected_folder) if '.gff' in file]


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

    # TODO: Implement that the core genes strand in saved.
    #  Implement that a summary of strand placement is given as output
    #  Implement that the strand is used in the consensus alignment to reverse complement gene alignments.
    with concurrent.futures.ProcessPoolExecutor(max_workers=15) as executor:
        print(f'{len(args.input_gffs)} GFF files to process')
        results = [executor.submit(segment_genome_content, gff, core_dict, low_freq_dict, acc_gene_dict, i)
                   for i, gff in enumerate(args.input_gffs)]

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
    time_calculator(time_start, time.time(), "searching gff files for core genes")


    ### FUNCTION ###
    # Determine the most common core gene synteny.
    # time_start = time.time()
    # Find the core gene synteny and possible core genes with alternative neighbours
    # TODO - FIX how the first and second gene in the genome is determined!
    # consensus_core_genome, \
    #     possible_rearrangement_genes, \
    #     core_path_coverage = determine_core_gene_consesus(core_neighbour_pairs,
    #                                                       merged_start_gene_clusters,
    #                                                       merged_second_gene_clusters, args.output_path)

    # Assign core-gene synteny types:
    # TODO: insert function for core gene synteny types

    # Identify alternative connections and their occurrence
    # alt_core_pairs, alt_core_pair_count, \
    #     core_genome_types, alt_core_comp_types = identify_rearrangements(consensus_core_genome,
    #                                                                      possible_rearrangement_genes,
    #                                                                      master_info_total,
    #                                                                      args.input_gffs)

    # TODO look at number of alternative core-genome pairs. if one give missing info, if two try to find size, if three or more too complex.
    # TODO Determine all alternative connections between core genes and their frequency. - SEE next line
    '''Use consensus core synteny to work out how rearrangements may have occured. If A and B are connected in the concesus and C and D
    are connected, then if A and C are connected, and B and D are found next to sequence breaks then it may be a possible recombination'''
    # rearrangement_predictions = characterise_rearrangements(alt_core_pairs, consensus_core_genome)

    # time_calculator(time_start, time.time(), "determining best core gene synteny")
    ################

    ### DETERMINE CONSENSUS ###
    # TODO Determine start cluster from possible consensus from complete genomes - else determine relative consensus from connections
    # TODO Do a greedy walk through a graph where all nodes are a core cluster and the edge weight is the number of times two core clusters are observed to be neighbours
    # TODO for-loop - Determine following clusters from connections to first cluster
    ###########################

    # TODO - constuct a multiple sequence alignment fasta from a possible input of alignments.
    #  If panaroo then look for: aligned_gene_sequences
    # if source_program == 'Panaroo':
    #     if alignment_folder:
    #         construct_consensus_alignment(consensus_core_genome.copy(), alignment_folder, core_neighbour_distance)

    ### DO CALCULATIONS ###
    # TODO mean number length between core genes
    # for neighbours in core_neighbour_distance:
    #     print(mean(core_neighbour_distance[neighbours]))
    #     print(std(core_neighbour_distance[neighbours]))

    # TODO SD of length between core genes
    # TODO plot -
    # TODO kruwalski plot
    #######################

    # TODO - Insert timer and printer to indicate the process and seperate printed messages nicely.
    # Constuct matrix over presence of alternative core pairs in genomes
    # The matrix is filled with 3, 2, 1, or 0, depending ong the evidence level for a core-core neighbour pair
    # 3 = Genomic evidence, where two genes are found connected by sequence
    # 2 = Strong suspicion. The consensus core neighbours had no sequence breaks
    # 1 = Weak evidence. at least one consensus core neighbour was next to a sequence break leading to a possible connection between then consensus neoghbours
    # 0 = No evidence for the alternative connection.
    # alt_core_pair_matrix = core_pair_matrix(core_genome_types, alt_core_comp_types,
    #                                         alt_core_pair_count, master_info_total, consensus_core_genome)

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
    # write_consensus_core_gene_synteny(consensus_core_genome)
    # write_core_gene_coverage(core_path_coverage)
    # write_alternative_core_gene_counts(alt_core_pair_count)
    # rearrangement_predictions #TODO - write output for rearrangement prediciton - First possibly make predictions for two alternative core pairs
    # write_core_gene_types(core_genome_types, alt_core_pair_matrix)

    time_calculator(time_start, time.time(), "writing output files")
    #####################

    ## WRITE PLOTS ##
    # Write the consensus core genome coverage plot.
    # consesus_genome_coverage(consensus_core_genome, core_path_coverage)

    ###################
    time_calculator(total_time_start, time.time(), "running the entire program")


if __name__ == "__main__":
    main()
