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

    parser.add_argument('-i_pan',
                        '--input_pangenome_folder',
                        help='Path to the folder produced by Panaroo or Roary',
                        required=True,
                        dest='input_pan')

    parser.add_argument('-o',
                        help='Path to where output files will be placed',
                        required=True,
                        dest='output_path')

    parser.add_argument('-p',
                        '--prefix',
                        help='Prefix for output files, if any is desired',
                        required=False,
                        dest='output_prefix')

    parser.add_argument('-a',
                        '--annotate_refound',
                        help='Flag to toggle off creation of new gff files, with annotation of refound genes.\n'
                             'Only done if input pangenome is detected as comming from Panaroo',
                        required=False,
                        default=True,
                        action='store_false',
                        dest='annotate')

    parser.add_argument('-q',
                        '--quiet',
                        help='Flag to toggle off printed info about the run',
                        required=False,
                        default=False,
                        action='store_true',
                        dest='quiet')

    args = parser.parse_args(args)

    return args
