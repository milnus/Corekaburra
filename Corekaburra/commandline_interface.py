import argparse
import sys

EXIT_COMMAND_LINE_ERROR = 2


def get_commandline_arguments(args):
    """
    Function that takes the input given to the commandline and passes it.
    will check for no input and '-help'
    :param args: List of input arguments given to the commandline
    :return: matched argument object for passing in main function.
    """
    # Set up parser
    parser = argparse.ArgumentParser(description='Welcome to Corekaburra!\n '
                                                 'Program to determine consensus core sequence from multiple genomes.\n'
                                                 'Outputs consensus core gene alignment, distance between core genes, '
                                                 'number of accessory genes between core genes and low frequency genes '
                                                 'between core genes') #TODO - Change

    parser.add_argument('-ig',
                        '--input_gffs',
                        help='Path to gff files used for pan-genome',
                        required=True,
                        metavar='file_1.gff ... file_n.gff',
                        dest='input_gffs',
                        nargs='+')

    parser.add_argument('-ip',
                        '--input_pangenome',
                        help='Path to the folder produced by Panaroo or Roary',
                        metavar='path/to/pan_genome',
                        required=True,
                        dest='input_pan')

    parser.add_argument('-cg',
                        '--complete_genomes',
                        help='text file containing names of genomes that are to be handled as complete genomes',
                        required=False,
                        metavar='complete_genomes.txt',
                        default=None,
                        dest='comp_genomes')

    parser.add_argument('-o',
                        '--output',
                        help='Path to where output files will be placed [default: current folder]',
                        required=False,
                        type=str,
                        metavar='path/to/output',
                        default='.',
                        dest='output_path')

    parser.add_argument('-p',
                        '--prefix',
                        help='Prefix for output files, if any is desired',
                        required=False,
                        default=None,
                        dest='output_prefix')

    parser.add_argument('-a',
                        '--no_annotate_refound',
                        help='Flag to toggle off the creation of new gff files, with annotation of refound genes.\n'
                             'Only done if input pangenome is detected as comming from Panaroo',
                        required=False,
                        default=True,
                        action='store_false',
                        dest='annotate')

    parser.add_argument('-c',
                        '--cpu',
                        help='Give max number of CPUs [default: 1]',
                        required=False,
                        metavar='int',
                        default=1,
                        type=int,
                        dest='cpu')

    logger_level = parser.add_mutually_exclusive_group()
    logger_level.add_argument('-l',
                              '--log',
                              help='Record program progress in for debugging purpose',
                              action='store_true',
                              default=False,
                              required=False)

    logger_level.add_argument('-q',
                              '--quiet',
                              help='Only print warnings',
                              action='store_true',
                              default=False,
                              required=False)

    # Check if any thing is given as input otherwise warn and print help
    if len(args) < 1:
        parser.print_help()
        sys.exit(EXIT_COMMAND_LINE_ERROR)
    elif '-help' in args:
        parser.print_help()
        sys.exit(0)
    if '--check' in args:
        sys.exit(1) #TODO write script that checks for dependencies!

    args = parser.parse_args(args)

    return args
