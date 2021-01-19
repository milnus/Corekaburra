import os


def check_gff_files(file_list):
    for file in file_list:
        if not os.path.isfile(file):
            raise FileNotFoundError(f'{file} can not be found!')
    return True


def define_input_source(folder):
    try:
        with open(os.path.join(folder, 'gene_presence_absence.csv'), 'r') as gene_pres_abs:
            if '"' in gene_pres_abs.readline():
                gene_pres_abs_file_path = os.path.join(folder, 'gene_presence_absence.csv')
                return "Roary", gene_pres_abs_file_path
            else:
                gene_pres_abs.close()
                gene_pres_abs_file_path = os.path.join(folder, 'gene_presence_absence_roary.csv')
                if os.path.isfile(gene_pres_abs_file_path):
                    return "Panaroo", gene_pres_abs_file_path

    except FileNotFoundError:
        raise FileNotFoundError('No gene presence absence file was found in given pan genome folder')


def check_gene_data(folder):
    if os.path.isfile(os.path.join(folder, 'gene_data.csv')):
        return True
    else:
        raise FileNotFoundError('gene_data.csv file could not be located in the given pan genome input folder.\n'
                                'Please give the -a flag to omit this step or locate the gene_data.csv file.')
