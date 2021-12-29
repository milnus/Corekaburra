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

try:
    from Corekaburra.commandline_interface import get_commandline_arguments
except ModuleNotFoundError:
    from commandline_interface import get_commandline_arguments

try:
    from Corekaburra.read_complete_genome_file import parse_complete_genome_file
except ModuleNotFoundError:
    from read_complete_genome_file import parse_complete_genome_file

try:
    from Corekaburra.check_inputs import define_pangenome_program, check_gene_data
except ModuleNotFoundError:
    from check_inputs import define_pangenome_program, check_gene_data

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
    if args.annotate and source_program is 'Panaroo':
        gene_data_path = check_gene_data(args.input_pan)
    if not args.quiet:
        print(f"Pan genome determined to come from {source_program}")
        print("All files found, let's move on!\n")
        print("--------------------------------------------------------------\n")

    # TODO - Make the program work with less than all files in the pangenome. Just make sure that all gff files supplied can be found in the pan genome. This will make is possible to look at hotspots and segments in different lineages
    check_gff_in_pan(args.input_gffs, input_pres_abs_file_path)


    # Construct output folder
    try:
        mkdir(args.output_path)
        if not args.quiet:
            print("Output folder constructed")
    except FileExistsError:
        if not args.quiet:
            print("Output folder exists")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
