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
    from Corekaburra.check_inputs import define_pangenome_program, check_gene_data, check_gff_in_pan, check_cutoffs
except ModuleNotFoundError:
    from check_inputs import define_pangenome_program, check_gene_data, check_gff_in_pan, check_cutoffs

try:
    from Corekaburra.parse_gene_presence_absence import read_gene_presence_absence
except ModuleNotFoundError:
    from parse_gene_presence_absence import read_gene_presence_absence

try:
    from Corekaburra.correct_gffs import prepair_for_reannotation
except ModuleNotFoundError:
    from correct_gffs import prepair_for_reannotation

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

try:
    from Corekaburra.output_writer_functions import master_info_writer, summary_info_writer, segment_writer, no_acc_segment_writer, non_core_contig_writer
except ModuleNotFoundError:
    from output_writer_functions import master_info_writer, summary_info_writer, segment_writer, no_acc_segment_writer, non_core_contig_writer

import sys
import pkg_resources

EXIT_INPUT_FILE_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_GFF_REANNOTATION_ERROR = 3
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

    # Check that low-frequency cutoff and core cutoff are as expected
    check_cutoffs(args.low_cutoff, args.core_cutoff)

    # TODO - Make Corekaburra take gzipped inputs
    # TODO - Add so that a single gff file can only be given as input once and not multiple times?

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
    else:
        gene_data_path = None

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
    # TODO - ATM the column with presence of gene in genomes is used to define what is core and not. Is it better to use the number of input gffs instead?
    #   - There are upsides to the current. You can use the same genome to find segments for two different populations with in the dataset using the same reference of core-genes
    #   - Making it depend on the input is not viable for comparing runs, even within the same pan-genome, when using different sets of gff files.
    # TODO - Some day it would be awesome to be able to provide a clustering/population structure which could divide genes into the 13 definitions outlined by Horesh et al. [DOI: 10.1099/mgen.0.000670]
    core_dict, low_freq_dict, acc_gene_dict = read_gene_presence_absence(input_pres_abs_file_path,
                                                                                         args.core_cutoff, args.low_cutoff, source_program,
                                                                                         args.input_gffs,
                                                                                         tmp_folder_path)

    if source_program == "Panaroo" and args.annotate:
        gene_data_dict, corrected_dict, args.input_gffs = prepair_for_reannotation(gene_data_path, args.output_path, args.input_gffs)
    else:
        gene_data_dict = None
        corrected_dict = None

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

    progress_counter = 0
    if len(args.input_gffs) > 10:
        progress_update = len(args.input_gffs) / 10
    else:
        progress_update = 1

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.cpu) as executor:
        print(f"\n------Start core region identification of given gff files-----")
        print(f'{len(args.input_gffs)} GFF files to process')

        results = [executor.submit(segment_genome_content, gff, core_dict, low_freq_dict, acc_gene_dict, comp_genomes,
                                   source_program, args.annotate, gene_data_dict, corrected_dict, tmp_folder_path, args.discard_gffs)
                   for gff in args.input_gffs]

        for output in concurrent.futures.as_completed(results):
            progress_counter += 1
            if progress_counter % progress_update == 0 or progress_counter == 1:
                print(
                    f"GFF file #{progress_counter} has been processed")

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
    master_summary_info = calculate_n_create_summaries(master_info_total, core_dict)

    ### WRITE OUTPUTS ###
    print(f"\n-----------------------Printing outputs-----------------------")
    # Write master information to output file
    time_start = time.time()
    master_info_writer(master_info_total, args.output_path, args.output_prefix, args.quiet)
    summary_info_writer(master_summary_info, args.output_path, args.output_prefix, args.quiet)
    if double_edge_segements is not None:
        segment_writer(double_edge_segements, args.output_path, args.output_prefix, args.quiet)
        no_acc_segment_writer(no_acc_segments, args.output_path, args.output_prefix, args.quiet)
    # TODO - Possibly output core gene graph. with segment annotations?
    #  - Print summary number of genes and names
    #  - Should we print a low-freq, placement?
    if len(non_core_contig_info)> 0:
        non_core_contig_writer(non_core_contig_info, args.output_path, args.output_prefix)

    # time_calculator(time_start, time.time(), "writing output files")

    # Finish up running
    # time_calculator(total_time_start, time.time(), "running the entire program")

    # Remove temporary database holding gff databases
    if os.path.isdir(tmp_folder_path):
        os.rmdir(tmp_folder_path)
    if args.discard_gffs:
        os.rmdir(os.path.join(args.output_path, 'Corrected_gff_files'))


if __name__ == '__main__':
    main()
