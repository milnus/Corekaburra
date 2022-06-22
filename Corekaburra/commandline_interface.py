import argparse
import sys

EXIT_COMMAND_LINE_ERROR = 2


def get_commandline_arguments(args, version):
    """
    Function that takes the input given to the commandline and passes it.
    will check for no input and '-help'
    :param args: List of input arguments given to the commandline
    :param version: Version of Corekaburra
    :return: matched argument object for passing in main function.
    """
    # Set up parser
    parser = argparse.ArgumentParser(description='Welcome to Corekaburra! '
                                                 'An extension to pan-genome analyses that summarise genomic regions '
                                                 'between core genes and segments of neighbouring core genes using '
                                                 'gene synteny from a set of input genomes and a pan-genome folder.',
                                     add_help=False)

    required = parser.add_argument_group('Required arguments')
    run_mods = parser.add_argument_group('Analysis modifiers')
    output_control = parser.add_argument_group('Output control')
    rem_args = parser.add_argument_group('Other arguments')

    required.add_argument('-ig',
                          '--input_gffs',
                          help='Path to gff files used for pan-genome',
                          metavar='file.gff',
                          required=True,
                          dest='input_gffs',
                          nargs='+')

    required.add_argument('-ip',
                          '--input_pangenome',
                          help='Path to the folder produced by Panaroo or Roary',
                          metavar='path/to/pan_genome',
                          required=True,
                          dest='input_pan')

    run_mods.add_argument('-cg',
                          '--complete_genomes',
                          help='text file containing names of genomes that are to be handled as complete genomes',
                          required=False,
                          metavar='complete_genomes.txt',
                          default=None,
                          dest='comp_genomes')

    # run_mods.add_argument('-a',
    #                       '--no_annotate_refound',
    #                       help='Flag to toggle off the creation of new gff files, with annotation of refound genes.\n'
    #                            'Only done if input pangenome is detected as coming from Panaroo',
    #                       required=False,
    #                       default=True,
    #                       action='store_false',
    #                       dest='annotate')

    run_mods.add_argument('-cc',
                          '--core_cutoff',
                          help='Percentage of isolates in which a core gene must be present [default: 1.0]',
                          metavar='1.0',
                          type=float,
                          default=1.0,
                          dest='core_cutoff')

    run_mods.add_argument('-lc',
                          '--low_cutoff',
                          help='Percentage of isolates where genes found in less than these are seen as low-frequency genes [default: 0.05]',
                          metavar='0.05',
                          type=float,
                          default=0.05,
                          dest='low_cutoff')

    output_control.add_argument('-o',
                                '--output',
                                help='Path to where output files will be placed [default: current folder]',
                                required=False,
                                type=str,
                                metavar='path/to/output',
                                default='.',
                                dest='output_path')

    output_control.add_argument('-p',
                                '--prefix',
                                help='Prefix for output files, if any is desired',
                                required=False,
                                default=None,
                                dest='output_prefix')

    # output_control.add_argument('-d',
    #                             '--discard_corrected',
    #                             help='Discard gff files corrected with refound genes identified by Panaroo - Only compativle if pan-genome comes from Panaroo [Default: Corrected files are kept]',
    #                             required=False,
    #                             default=False,
    #                             action='store_true',
    #                             dest='discard_gffs')

    rem_args.add_argument('-c',
                          '--cpu',
                          help='Give max number of CPUs [default: 1]',
                          required=False,
                          metavar='int',
                          default=1,
                          type=int,
                          dest='cpu')

    logger_level = rem_args.add_mutually_exclusive_group()
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

    rem_args.add_argument('-h',
                          '--help',
                          action='help',
                          help='Show help function')

    rem_args.add_argument('-v',
                          '--version',
                          action='version',
                          version=f'Corekaburra {version}')

    # Check if any thing is given as input otherwise warn and print help
    if len(args) < 1:
        parser.print_help()
        sys.exit(EXIT_COMMAND_LINE_ERROR)
    elif '-help' in args:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args(args)

    return args


if __name__ == '__main__':
    get_commandline_arguments([], 666)