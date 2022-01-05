'''
Module      : Main
Description : The main entry point for the Corekaburra.
Copyright   : (c) Magnus Ganer Jespersen, 11 Oct 2021 
License     : MIT 
Maintainer  : magnus.ganer.j@gmail.com 
Portability : POSIX

The program reads one or more input FASTA files. For each file it computes a
variety of statistics, and then prints a summary of the statistics as output. # TODO - Change description
'''

import os
import logging
import time
import concurrent.futures

try:
    from Corekaburra.commandline_interface import get_commandline_arguments
except ModuleNotFoundError:
    from commandline_interface import get_commandline_arguments

try:
    from Corekaburra.read_complete_genome_file import parse_complete_genome_file
except ModuleNotFoundError:
    from read_complete_genome_file import parse_complete_genome_file

try:
    from Corekaburra.check_inputs import define_pangenome_program, check_gene_data, check_gff_in_pan
except ModuleNotFoundError:
    from check_inputs import define_pangenome_program, check_gene_data, check_gff_in_pan

try:
    from Corekaburra.parse_gene_presence_absence import read_gene_presence_absence
except ModuleNotFoundError:
    from parse_gene_presence_absence import read_gene_presence_absence

try:
    from Corekaburra.gff_parser import segment_genome_content
except ModuleNotFoundError:
    from gff_parser import segment_genome_content

try:
    from Corekaburra.merge_dicts import merge_dicts_lists, merge_dicts_counts
except ModuleNotFoundError:
    from merge_dicts import merge_dicts_lists, merge_dicts_counts

try:
    from Corekaburra.consesus_core_genome import determine_genome_segments
except ModuleNotFoundError:
    from consesus_core_genome import determine_genome_segments

try:
    from Corekaburra.summary_table import calculate_n_create_summaries
except ModuleNotFoundError:
    from summary_table import calculate_n_create_summaries


from argparse import ArgumentParser
from math import floor
import sys
import pkg_resources

EXIT_INPUT_FILE_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTA_FILE_ERROR = 3
DEFAULT_MIN_LEN = 0
DEFAULT_VERBOSE = False
PROGRAM_NAME = "Corekaburra"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def init_logging(debug_log, quiet, out_path):
    """
    initialise the logging file, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv
    :param debug_log: Bool indicating if log is a debug log
    :param quiet: Bool indicating if logging should be kept minimal
    :param out_path: Output path for the program, and where log will be placed
    :return: Logger object
    """
    if debug_log:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    else:
        level = logging.INFO

    # Construct logger logging to file
    file_logger = logging.getLogger(__name__)
    file_logger.setLevel(level)

    formatter = logging.Formatter('[%(asctime)s] %(levelname)s - %(module)s - %(message)s',
                                  datefmt="%Y-%m-%dT%H:%M:%S%z")

    file_handler = logging.FileHandler(os.path.join(out_path, 'Corekaburra.log'))
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    file_logger.addHandler(file_handler)

    # Log command-line argument and debug line for Corekaburra start
    file_logger.info(f"command line: {' '.join(sys.argv)}")

    return file_logger


def stream_logging(file_logger):
    """
    Function adding in stream logging following initial logging
    :param file_logger: Logger object
    :return: Logger object with added stream logging
    """
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    file_logger.addHandler(stream_handler)

    file_logger.info('Processing started')

    return file_logger


def main():
    """
    This is the main function for running Corekaburra.
    Requires a pan-genome folder from an allowed program, along with GFF3 files used for producing the pan-genome
    """
    total_time_start = time.time()

    # get arguments from the commandline
    args = get_commandline_arguments(sys.argv[1:])

    # TODO - Add in function(s) that will check all files to not be empty. - Andrew?
    # TODO - Make Corekaburra take gzipped inputs

    # Check the presence of provided complete genomes among input GFFs
    if args.comp_genomes is not None:
        comp_genomes = parse_complete_genome_file(args.comp_genomes, args.input_gffs)
    else:
        comp_genomes = None

    # Check source program from pan-genome and presence of nessecary files
    if not args.quiet:
        print("\n----Checking presence of input files in pan genome folder----\n")

    # Check if Panaroo or Roary input folder is given
    source_program, input_pres_abs_file_path = define_pangenome_program(args.input_pan)

    # Check if gene_data file is present if Panaroo input is given an gffs should be annotated
    if args.annotate and source_program == 'Panaroo':
        gene_data_path = check_gene_data(args.input_pan)
    if not args.quiet:
        print(f"Pan genome determined to come from {source_program}")
        print("All files found, let's move on!\n")
        print("--------------------------------------------------------------\n")

    # TODO - Make the program work with less than all files in the pangenome. Just make sure that all gff files supplied can be found in the pan genome. This will make is possible to look at hotspots and segments in different lineages
    check_gff_in_pan(args.input_gffs, input_pres_abs_file_path)


    # Construct output folder
    try:
        os.mkdir(args.output_path)
        if not args.quiet:
            print("Output folder constructed")
    except FileExistsError:
        if not args.quiet:
            print("Output folder exists")

    # Construct temporary folder:
    # TODO - check that the temporary folder does not exist and that the user does not have a folder with same name already. (Maybe use a time stamp for the start to make it unique.)
    tmp_folder_path = os.path.join(args.output_path, 'Corekaburra_tmp')
    os.mkdir(tmp_folder_path)

    ## Read in gene presence absence file
    time_start = time.time()
    # TODO - Add the user specified thresholds for core and low frequency genes.
    core_dict, low_freq_dict, acc_gene_dict, attribute_dict = read_gene_presence_absence(input_pres_abs_file_path,
                                                                                         1, 0.05, source_program,
                                                                                         args.input_gffs,
                                                                                         tmp_folder_path)

    # TODO - Add this into the multiprocessing loop to not doubble files
    # TODO - Add a user command to keep and discard the corrected files (But still using them - Make mutually exclusive with -a option)
    # Add in the refound genes into the gff files and print the corrected GFF files.
    # if source_program == "Panaroo" and args.annotate:
    #     time_start = time.time()
    #     print(f"\n----------Adding in refound annotations for gff files---------")
    #
    #     corrected_folder = correct_gffs(args.input_gffs, gene_data_path, args.output_path, attribute_dict,
    #                                     temp_folder_path)
    #
    #     args.input_gffs = [join(corrected_folder, file) for file in listdir(corrected_folder) if '.gff' in file]
    #     if not args.quiet:
    #         time_calculator(time_start, time.time(), "add refound annotations to gff files")

    # Loop over all gffs and extract info from each of them.
    time_start = time.time()
    # Initialise dictionaries to contain results from all gff files
    core_neighbour_pairs = {}
    core_neighbour_distance = {}
    core_neighbour_accessory_count = {}
    core_neighbour_low_freq = {}
    master_info_total = {}
    non_core_contig_info = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=15) as executor: # TODO - change the max workers to the user specified number
        print(f"\n------Start core region identification of given gff files-----")
        print(f'{len(args.input_gffs)} GFF files to process')

        results = [executor.submit(segment_genome_content, gff, core_dict, low_freq_dict, acc_gene_dict, i, comp_genomes)
                   for i, gff in enumerate(args.input_gffs)]

        for output in concurrent.futures.as_completed(results):
            # Split the outputs
            core_pairs, distance, acc_count, \
            low_freq, master_info_return, \
            core_less_contigs_return = output.result()

            # Merge results into single/master dictionaries
            core_neighbour_pairs = merge_dicts_counts(core_neighbour_pairs, core_pairs)
            core_neighbour_distance = merge_dicts_lists(core_neighbour_distance, distance)
            core_neighbour_accessory_count = merge_dicts_lists(core_neighbour_accessory_count, acc_count)
            core_neighbour_low_freq = merge_dicts_lists(core_neighbour_low_freq, low_freq)
            master_info_total.update(master_info_return)
            non_core_contig_info.update(core_less_contigs_return)
    #
    # time_calculator(time_start, time.time(), "searching gff files for core genes")

    print(f"\n--------------Identifying segments in pan genome--------------")
    time_start = time.time()
    # Count number of unique accessory genes inserted into a core-core region across the genomes
    acc_region_count = {key: len(set(core_neighbour_low_freq[key])) for key in core_neighbour_low_freq}
    # Count number of unique low frequency genes inserted into a core-core region across the genomes
    low_frew_region_count = {key: len(set(core_neighbour_accessory_count[key])) for key in
                             core_neighbour_accessory_count}

    # Combine the accessory and low frequency counts:
    combined_acc_gene_count = {key: low_frew_region_count[key] + acc_region_count[key] for key in low_frew_region_count}

    double_edge_segements, no_acc_segments = determine_genome_segments(core_neighbour_pairs, combined_acc_gene_count,
                                                                       len(args.input_gffs), core_dict)

    # time_calculator(time_start, time.time(), "identifying segments in pan genome")

    # Produce dict containing summarised information from master info.
    master_summary_info = calculate_n_create_summaries(master_info_total)

    ### WRITE OUTPUTS ###
    print(f"\n-----------------------Printing outputs-----------------------")
    # Write master information to output file
    time_start = time.time()
    master_info_writer(master_info_total, args.output_path, args.output_prefix, args.quiet)
    summary_info_writer(master_summary_info, args.output_path, args.output_prefix, args.quiet)
    # TODO - Contruct output for segments - parent column.
    segment_writer(double_edge_segements, args.output_path, args.output_prefix, args.quiet)
    no_acc_segment_writer(no_acc_segments, args.output_path, args.output_prefix, args.quiet)
    # print(non_core_contig_info) TODO - Print core less contigs.
    # TODO - Possibly output core gene graph. with segment annotations?

    # time_calculator(time_start, time.time(), "writing output files")

    # Finish up running
    # time_calculator(total_time_start, time.time(), "running the entire program")

    # Remove temporary database holding gff databases
    # TODO - Implement a nice crash function where the temporary folder is removed not to cause unessecary frustration for the user when trying to rerun the program. - do so in nice exit function
    # print(isdir(temp_folder_path))
    # if isdir(temp_folder_path):
    #     rmdir(temp_folder_path)


if __name__ == '__main__':
    main()
