'''
Module      : Main
Description : The main entry point for the Corekaburra.
Copyright   : (c) Magnus Ganer Jespersen, 11 Oct 2021 
License     : MIT 
Maintainer  : magnus.ganer.j@gmail.com 
Portability : POSIX

Corekaburra looks at the gene synteny across genomes used to build a pan-genome. Using syntenic information Corekaburra
identifies regions between core gene clusters. Regions are described in terms of their content of accessory gene clusters
and distance between core genes. Information from neighboring core genes is further used to identify stretches of core
gene clusters throughout the pan-genome that appear in all genomes given as input. Corekaburra is compatible with outputs
from standard pan-genome pipelines: Roary and Panaroo.
'''

import os
import logging
import tempfile
import time
import concurrent.futures
from networkx import write_gml

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
EXIT_SEGMENT_IDENTIFICATION_ERROR = 4
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

    formatter = logging.Formatter('[%(asctime)s] - %(levelname)s - %(module)s: %(message)s',
                                  datefmt="%Y-%m-%d %H:%M:%S%z")

    file_handler = logging.FileHandler(os.path.join(out_path, 'Corekaburra.log'))
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    file_logger.addHandler(file_handler)

    # Log command-line argument and debug line for Corekaburra start
    file_logger.info(f"command line: {' '.join(sys.argv)}")


    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    file_logger.addHandler(stream_handler)

    file_logger.info('\n----------------------Processing started----------------------\n')

    return file_logger


def main():
    """
    This is the main function for running Corekaburra.
    Requires a pan-genome folder from an allowed program, along with GFF3 files used for producing the pan-genome
    """
    total_time_start = time.time()

    inital_check_time_start = time.time()

    # get arguments from the commandline
    args = get_commandline_arguments(sys.argv[1:], PROGRAM_VERSION)

    # Construct output folder
    try:
        os.mkdir(args.output_path)
    except FileExistsError:
        pass

    # Run initialisation of logger:
    logger = init_logging(args.log, args.quiet, args.output_path)

    # Check that low-frequency cutoff and core cutoff are as expected
    check_cutoffs(args.low_cutoff, args.core_cutoff, logger)

    # Check the presence of provided complete genomes among input GFFs
    if args.comp_genomes is not None:
        comp_genomes = parse_complete_genome_file(args.comp_genomes, args.input_gffs, logger)
    else:
        comp_genomes = None

    # Check if Panaroo or Roary input folder is given
    source_program, input_pres_abs_file_path = define_pangenome_program(args.input_pan, logger)

    # Check that all GFF files given can be found in the pan-genome
    check_gff_in_pan(args.input_gffs, input_pres_abs_file_path, logger)

    # Construct temporary folder:
    tmp_folder_path = tempfile.TemporaryDirectory()

    logger.info('Initial checks successful\n')
    inital_check_time_end = time.time()

    ## Read in gene presence absence file
    time_start_read_files = time.time()

    # TODO - Add in so that the user can give a list of genes that they wish to use as 'core genes'
    core_dict, low_freq_dict, acc_gene_dict = read_gene_presence_absence(input_pres_abs_file_path, args.core_cutoff,
                                                                         args.low_cutoff, source_program,
                                                                         args.input_gffs, tmp_folder_path.name,
                                                                         logger)


    time_end_read_files = time.time()
    time_start_passing_gffs = time.time()

    # Loop over all Gffs and extract info from each of them.
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
        progress_update = int(len(args.input_gffs) / 10)
    else:
        progress_update = 1

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.cpu) as executor:
        logger.info(f"------Start core region identification of given gff files-----\n")
        logger.info(f'{len(args.input_gffs)} GFF files to process')
        results = [executor.submit(segment_genome_content, gff, core_dict, low_freq_dict, acc_gene_dict, comp_genomes)
                   for gff in args.input_gffs]

        for output in concurrent.futures.as_completed(results):
            progress_counter += 1
            if progress_counter % progress_update == 0 or progress_counter == 1:
                logger.info(f"GFF file #{progress_counter} has been processed")

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

    time_end_passing_gffs = time.time()
    time_start_segments_search = time.time()

    # Count number of unique accessory genes inserted into a core-core region across the genomes
    acc_region_count = {key: len(set(core_neighbour_low_freq[key])) for key in core_neighbour_low_freq}
    # Count number of unique low frequency genes inserted into a core-core region across the genomes
    low_frew_region_count = {key: len(set(core_neighbour_accessory_count[key])) for key in
                             core_neighbour_accessory_count}

    # Combine the accessory and low frequency counts:
    combined_acc_gene_count = {key: low_frew_region_count[key] + acc_region_count[key] for key in low_frew_region_count}

    double_edge_segements, no_acc_segments, core_graph = determine_genome_segments(core_neighbour_pairs, combined_acc_gene_count,
                                                                       len(args.input_gffs), core_dict, args.cpu, logger)

    time_end_segments_search = time.time()

    # Produce dict containing summarised information from master info.
    logger.debug("Commence on calculating summary output")
    master_summary_info = calculate_n_create_summaries(master_info_total, core_dict)

    ### WRITE OUTPUTS ###
    logger.debug("-----------------------Printing outputs-----------------------")
    # Write master information to output file
    time_start = time.time()
    logger.debug("Master outputs")
    master_info_writer(master_info_total, args.output_path, args.output_prefix)

    logger.debug("Summary output")
    summary_info_writer(master_summary_info, args.output_path, args.output_prefix)

    if double_edge_segements:
        logger.debug("Segment output")
        segment_writer(double_edge_segements, args.output_path, args.output_prefix)

        logger.debug("No Accessory segment output")
        no_acc_segment_writer(no_acc_segments, args.output_path, args.output_prefix)

        logger.debug("Writing core gene graph")
        graph_name = f'{args.output_prefix}_core_gene_graph.gml' if args.output_prefix is not None else 'core_gene_graph.gml'
        write_gml(core_graph, path=os.path.join(args.output_path, graph_name))

    if len(non_core_contig_info) > 0:
        logger.debug("Non-core contig output")
        non_core_contig_writer(non_core_contig_info, args.output_path, args.output_prefix)

    # Finish up running
    total_time = round(time.time() - total_time_start, 1)
    initial_time = round(inital_check_time_end - inital_check_time_start, 1)
    read_fies_time = round(time_end_read_files - time_start_read_files, 1)
    passing_gffs_time = round(time_end_passing_gffs - time_start_passing_gffs, 1)
    segment_search_time = round(time_end_segments_search - time_start_segments_search)

    logger.debug("-----------------------Time used in run-----------------------")
    logger.debug(f"Total time used: {total_time}s")
    logger.debug(f"Initial check time: {initial_time}s")
    logger.debug(f"Reading pan-genome files time: {read_fies_time}s")
    logger.debug(f"Passing over Gff files time: {passing_gffs_time}s")
    logger.debug(f"Searching for segments time: {segment_search_time}s")


if __name__ == '__main__':
    main()
