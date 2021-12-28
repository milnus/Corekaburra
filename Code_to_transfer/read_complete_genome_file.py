import os


def parse_complete_genome_file(complete_genome_file, gff_files):
    # Read the file and all lines (complete genomes given)
    with open(complete_genome_file, 'r') as genome_file:
        complete_genomes = genome_file.readlines()
        complete_genomes = [name.strip().replace('.gff', '') for name in complete_genomes]
        complete_genomes = [os.path.basename(name) for name in complete_genomes]

    # Take input gffs and remove path to the file
    gffs = [os.path.basename(gff).replace('.gff', '') for gff in gff_files]

    # check that all complete genomes are in the input gffs
    complete_genome_status = all(complete_genome in gffs for complete_genome in complete_genomes)

    # If the complete genomes are found, return. Else try to remove the extension of input files and compare
    if complete_genome_status:
        return complete_genomes

    NotImplementedError('Reading the file of complete genomes encountered an error that has not been expected!\n'
                        'Please report this and give the file you passed. Cheers!')

if __name__ == '__main__':
    pass
