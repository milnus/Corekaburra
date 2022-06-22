import os
import warnings

try:
    from Corekaburra.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error
EXIT_INPUT_FILE_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2


def check_cutoffs(low_cutoff, core_cutoff, logger):
    """
    Function to check the given cutoffs are legal, otherwise provide more info.
    :param low_cutoff: Cutoff for low-frequency genes
    :param core_cutoff: Cutoff for core genes
    :param logger: Program logger
    :return: Nothing
    """
    if 0 <= low_cutoff < core_cutoff <= 1:
        logger.debug(f'User provided cutoffs: {low_cutoff = } and {core_cutoff = } were accepted')
        return
    else:
        exit_with_error('Something is wrong with cutoffs for core and low-frequency genes!\n'
                        'Make sure the cutoff for core genes is larger than for low-frequency, and is >0 or =1.\n'
                        'Also make sure that the low-frequency gene cutoff is either equal to 0 or <1',
                        EXIT_COMMAND_LINE_ERROR, logger)


def define_pangenome_program(folder, logger):
    """
    Function to examine if input pan genome folder stems from Roary or Panaroo.
    :param folder: Input folder provided as pan-genome folder
    :param logger: Program logger
    :return: The name of the program from which the pangenome is suspected to come from
    """

    logger.debug("----Checking presence of input files in pan genome folder----")

    try:
        if os.path.isfile(os.path.join(folder, 'gene_presence_absence.csv')):
            # See if input is from Roary
            with open(os.path.join(folder, 'gene_presence_absence.csv'), 'r') as gene_pres_abs:
                if '"' in gene_pres_abs.readline():
                    gene_pres_abs_file_path = os.path.join(folder, 'gene_presence_absence.csv')
                    logger.debug('Pan-genome_program detected to be: Roary')
                    return "Roary", gene_pres_abs_file_path

        # See if input is from Panaroo
        gene_pres_abs_file_path = os.path.join(folder, 'gene_presence_absence_roary.csv')
        if os.path.isfile(gene_pres_abs_file_path):
            logger.debug('Pan-genome_program detected to be: Panaroo')
            return "Panaroo", gene_pres_abs_file_path
        else:
            exit_with_error('No gene presence/absence file was found in given pan-genome folder', EXIT_INPUT_FILE_ERROR, logger)

    except FileNotFoundError:
        exit_with_error('No gene presence/absence file was found in given pan-genome folder', EXIT_INPUT_FILE_ERROR)


def check_gene_data(folder, logger):
    """
    Check if the gene_data.csv file is present in the folder from a Panaroo pan-genome run
    :param folder: Input folder provided as pan-genome folder
    :param logger: Program logger
    :return: Path to the identified gene_data.csv file
    """

    logger.debug('Identify gene_data.csv from Panaroo')

    if os.path.isfile(os.path.join(folder, 'gene_data.csv')):
        logger.debug('Gene_data.csv from Panaroo found!')
        return os.path.join(folder, 'gene_data.csv')
    else:
        exit_with_error('gene_data.csv file could not be located in the given pan genome input folder.\n'
                        'Please give the -a flag to omit this step or locate the gene_data.csv file.',
                        EXIT_INPUT_FILE_ERROR, logger)


def check_gff_in_pan(file_list, gene_presence_absence_path, logger):
    """
    Function to check if all given Gff files are in the given pan-genome
    :param file_list: List of Gff file paths
    :param gene_presence_absence_path: File path to the identified gene_presence_absence_file
    :param logger: Program logger
    :return: Bool used for unit testing
    """
    with open(gene_presence_absence_path, 'r') as pan_file:
        # Read the first line of the gene_presence_absence and extract the genome names
        pan_header_line = pan_file.readline()
        pan_header_line = pan_header_line.strip().split(',')
        genome_names = pan_header_line[14:]
        genome_names = [name.replace('"', '') for name in genome_names]

    file_list = [os.path.basename(file) for file in file_list]
    file_list_no_suffix = [file.rstrip('.gz') for file in file_list]
    file_list_no_suffix = [file.rstrip('.gff') for file in file_list_no_suffix]

    # Check if all or subset of GFFs from pan genome have been supplied,
    # if only a subset then raise warning
    if set(file_list).issubset(genome_names) or set(file_list_no_suffix).issubset(genome_names):
        if len(file_list) < len(genome_names):
            logger.info("\nNot all gff in pan genome given as input. I will run with it but are you sure this is deliberate?\n")

        return True  # True used for unit testing

    # Exit with error is not all inputs can be found in the pan-genome presence absence file
    exit_with_error('Unexpected occurrence in the matching of input GFF files and the pan genome presence/absence file',
                    EXIT_INPUT_FILE_ERROR, logger)
