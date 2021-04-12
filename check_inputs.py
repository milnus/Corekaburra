import os
import warnings


def check_gff_files(file_list):
    for file in file_list:
        if not os.path.isfile(file):
            raise FileNotFoundError(f'{file} can not be found!')
    return True


def define_input_source(folder):
    """ Function to examine if input pan genome folder stems from Roary or Panaroo.
    Returns the program that is the source and the path to the right gene presence/absence file """
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

    except FileNotFoundError:
        raise FileNotFoundError('No gene presence absence file was found in given pan genome folder')


def check_gene_data(folder):
    """ Check if the gene_data.csv file is present in the folder from a Panaroo pan genome run. """
    if os.path.isfile(os.path.join(folder, 'gene_data.csv')):
        return os.path.join(folder, 'gene_data.csv')
    else:
        raise FileNotFoundError('gene_data.csv file could not be located in the given pan genome input folder.\n'
                                'Please give the -a flag to omit this step or locate the gene_data.csv file.')


def check_gene_alignments(folder, core_gene_dict):
    """ Check if the folder containing alignments for genes is available in the Panaroo folder """
    alignment_folder = os.path.join(folder, 'aligned_gene_sequences')
    if os.path.isdir(alignment_folder):

        # Get a list of unique core genes
        genome_dicts = [genome for genome in core_gene_dict.values()]
        core_genes = [genome.values() for genome in genome_dicts]
        core_genes = set([core_gene for genome in core_genes for core_gene in genome])

        # Check that all core genes has an alignment
        alignments_missing = [core_gene for core_gene in core_genes
                              if not os.path.isfile(os.path.join(alignment_folder, f'{core_gene}.aln.fas'))]

        if len(alignments_missing) == 0:
            return alignment_folder
        else:
            warnings.warn(f'Not all core genes have alignments. '
                          f'Genes missing alignments: {alignments_missing}')
            return False

    else:
        return False
