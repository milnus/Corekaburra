import os
import warnings

try:
    from Corekaburra.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error
EXIT_INPUT_FILE_ERROR = 1


def check_gff_files(file_list):
    for file in file_list:
        if not os.path.isfile(file):
            raise FileNotFoundError(f'{file} can not be found!')
    return True


def check_gff_in_pan(file_list, gene_presence_absence_path):
    with open(gene_presence_absence_path, 'r') as pan_file:
        # Read the first line of the gene_presence_absence and extract the genome names
        pan_header_line = pan_file.readline()
        pan_header_line = pan_header_line.strip().split(',')
        genome_names = pan_header_line[14:]


    file_list = [os.path.basename(file) for file in file_list]
    file_list_no_suffix = [file.rstrip('.gff') for file in file_list]

    # Check if all or subset of GFFs from pan genome have been supplied,
    # if only a subset then raise warning
    if set(file_list).issubset(genome_names) or set(file_list_no_suffix).issubset(genome_names):
        if len(file_list) < len(genome_names):
            warnings.warn(
                "Not all gff in pan genome given as input. I will run with it but are you sure this is deliberate?")

        return True  # True used for unit testing

    raise FileNotFoundError('Unexpected occurrence in the matching of input GFF files and the pan genome presence/absence file')


def define_pangenome_program(folder):
    """
    Function to examine if input pan genome folder stems from Roary or Panaroo.
    :param folder: Input folder provided as pan-genome folder.
    :return: The name of the program from which the pangenome is suspected to come from
    """

    try:
        if os.path.isfile(os.path.join(folder, 'gene_presence_absence.csv')):
            # See if input is from Roary
            with open(os.path.join(folder, 'gene_presence_absence.csv'), 'r') as gene_pres_abs:
                if '"' in gene_pres_abs.readline():
                    gene_pres_abs_file_path = os.path.join(folder, 'gene_presence_absence.csv')
                    return "Roary", gene_pres_abs_file_path

        # See if input is from Panaroo
        gene_pres_abs_file_path = os.path.join(folder, 'gene_presence_absence_roary.csv')
        if os.path.isfile(gene_pres_abs_file_path):
            return "Panaroo", gene_pres_abs_file_path
        else:
            exit_with_error('No gene presence/absence file was found in given pan-genome folder', EXIT_INPUT_FILE_ERROR)

    except FileNotFoundError:
        exit_with_error('No gene presence/absence file was found in given pan-genome folder', EXIT_INPUT_FILE_ERROR)


def check_gene_data(folder):
    """ Check if the gene_data.csv file is present in the folder from a Panaroo pan genome run. """
    if os.path.isfile(os.path.join(folder, 'gene_data.csv')):
        return os.path.join(folder, 'gene_data.csv')
    else:
        raise FileNotFoundError('gene_data.csv file could not be located in the given pan genome input folder.\n'
                                'Please give the -a flag to omit this step or locate the gene_data.csv file.')
