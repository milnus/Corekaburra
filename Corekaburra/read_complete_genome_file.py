import os

try:
    from Corekaburra.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error
EXIT_INPUT_FILE_ERROR = 1


def parse_complete_genome_file(complete_genome_file, gff_files, logger):
    """
    Function to check if all genomes given as complete genomes can be found in the pan genome.
    :param complete_genome_file:
    :param gff_files:
    :param logger: Logger for program
    :return: a list of the base name of the complete genomes.
    """

    # Read the file and all lines (complete genomes given)
    with open(complete_genome_file, 'r') as genome_file:
        complete_genomes = genome_file.readlines()
        complete_genomes = [name.strip().replace('.gz', '') for name in complete_genomes]
        complete_genomes = [name.replace('.gff', '') for name in complete_genomes]
        complete_genomes = [os.path.basename(name) for name in complete_genomes]

    # Take input gffs and remove path to the file
    gffs = [os.path.basename(gff).replace('.gff', '').replace('.gz', '') for gff in gff_files]

    # check that all complete genomes are in the input gffs
    complete_genome_status = all(complete_genome in gffs for complete_genome in complete_genomes)

    # If the complete genomes are found, return a list of complete genomes
    if complete_genome_status:
        logger.debug(f'complete genomes: {complete_genomes = } were identified and accepted')
        return complete_genomes
    else:
        logger.debug(f'complete genomes: {complete_genomes = } were identified but not accepted!')
        exit_with_error('Genome given in Complete genomes was not identified in pan-genome!', EXIT_INPUT_FILE_ERROR, logger)


if __name__ == '__main__':
    pass
