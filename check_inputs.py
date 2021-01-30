import os


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

            gene_pres_abs.close()

        # See if input is from Panaroo
        print('PAN')
        gene_pres_abs_file_path = os.path.join(folder, 'gene_presence_absence_roary.csv')
        if os.path.isfile(gene_pres_abs_file_path):
            return "Panaroo", gene_pres_abs_file_path

    except FileNotFoundError:
        raise FileNotFoundError('No gene presence absence file was found in given pan genome folder')


def check_gene_data(folder):
    """ Check if the gene_data.csv file is pressent in the folder from a Panaroo pan genome run. """
    if os.path.isfile(os.path.join(folder, 'gene_data.csv')):
        return True
    else:
        raise FileNotFoundError('gene_data.csv file could not be located in the given pan genome input folder.\n'
                                'Please give the -a flag to omit this step or locate the gene_data.csv file.')
