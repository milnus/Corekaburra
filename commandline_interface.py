import argparse


def get_commandline_arguments(args):
    # Set up parser
    parser = argparse.ArgumentParser(description='Welcome to Coredial!\n '
                                                 'Program to determine consensus core sequence from multiple genomes.\n'
                                                 'Outputs consensus core gene alignment, distance between core genes, '
                                                 'number of accessory genes between core genes and low frequency genes '
                                                 'between core genes')

    parser.add_argument('-i_gffs',
                        '--input_gffs',
                        help='Path to gff files used for pan-genome',
                        required=True,
                        dest='input_gffs',
                        nargs='+')

    parser.add_argument('-i_pres_abs',
                        '--input_presence_absence',
                        help='Path to gene presence/absence file from a pan-genome',
                        required=True,
                        dest='input_pres_abs')

    parser.add_argument('-o',
                        help='Path to where output files will be placed',
                        required=True,
                        dest='output_path')

    parser.add_argument('-p',
                        '--prefix',
                        help='Prefix for output files, if any is desired',
                        required=False,
                        dest='output_prefix')

    args = parser.parse_args(args)

    return args
