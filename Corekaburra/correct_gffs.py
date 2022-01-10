from Bio import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
from gffutils.gffwriter import GFFWriter
import os
import concurrent.futures
# from Bio.Blast import
from time import time

try:
    from Corekaburra.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error
EXIT_GFF_REANNOTATION_ERROR = 3


def read_gene_data(gene_data_file):
    """
    Function to read the gene_data.csv file outputted by Panaroo and
    :param gene_data_file: File path to the gene_data.csv file
    :return: A dict of genomes with their refound genes
    """

    # Construct dictionary to hold refound genes and sequences for these
    gene_data_dict = {}

    # Read the gene_data.csv file and record all refound genes
    with open(gene_data_file, 'r') as gene_data:

        for line in gene_data.readlines():
            # Split read line at commas
            line = line.split(',')

            # Check if refound gene
            if 'refound' in line[2]:
                # Try to add the refound gene to the gene_data dict as a second key, value being the DNA sequence,
                # if the first key (genome) is not found in gene_data dict,
                # then construct dict for the genome and add the gene
                try:
                    gene_data_dict[line[0]][line[2]] = [line[5], line[6], line[7].strip()]
                except KeyError:
                    gene_data_dict[line[0]] = {line[2]: [line[5], line[6], line[7].strip()]}

    return gene_data_dict


def prepair_for_reannotation(gene_data_path, output_folder, gffs, logger):
    """
    Function for creating an output folder for corrected genomes, check if any are present, and if then which.
    :param gene_data_path: Path to the gene_data.csv file from Panaroo
    :param output_folder: Folder designated as the output folder for Corekaburra
    :param gffs: List of file-paths to gff files.
    :param logger: Program logger

    :return gene_data_dict: Dict containing the information expected from the gene_data.csv file
    :return corrected_gff_out_dir: File path to the created or identified directory of corrected gff files
    :return gffs: List of gff files, some may be altered to be the corrected verison from prior runs/
    """

    logger.debug('Initialise structures for reannotating genes found by Panaroo')

    # Read Gene_data.csv file into dict with a dict of refound genes for each genome
    gene_data_dict = read_gene_data(gene_data_path)

    # Construct directory to hold corrected gff files:
    corrected_gff_out_dir = os.path.join(output_folder, 'Corrected_gff_files')
    # Try and construct folder,
    # if present check if content matches input to avoid process
    try:
        os.mkdir(corrected_gff_out_dir)
    except FileExistsError:
        corrected_folder_content = os.listdir(corrected_gff_out_dir)

        gff_names = [os.path.basename(gff) for gff in gffs]

        corrected_files = [file for file in corrected_folder_content if
                           f'{file.split("_corrected")[0]}.gff' in gff_names]

        if len(corrected_files) > 0:
            gffs = [file for file in gff_names if f'{file.replace(".gff", "")}_corrected.gff' not in corrected_files]
            gffs = gffs + corrected_files

    return gene_data_dict, corrected_gff_out_dir, gffs


def extract_genome_fasta(gff_name):
    """
    Function to read and extract information from a gff3 file
    :param gff_name: File path to a gff file

    :return genome_fasta_dict: Dict over contig names as keys, and values being the contig sequence
    :return largest_locus_tag: the largest locus_tag identified in gff file
    :return header_lines: the header lines proceeding annotations.
    """

    # Initialise the two return variables a dict for refound gene's annotaton and the largest locus_tag
    genome_fasta_dict = {}
    largest_locus_tag = ''
    header_lines = []

    # Open the gff file and indicate that the FASTA sequence has not beed reached
    with open(gff_name, 'r') as gff_file:
        found_fasta = False

        # Go through gff file and find and read the fasta seuqence at the end after the ##FASTA mark
        for line in gff_file.readlines():
            if found_fasta:
                # Check if line is fasta header,
                # If then construct new dict key if not append sequence to current key
                if '>' in line:
                    line = line.split(' ')[0]
                    try:
                        current_contig = line.strip()
                        current_contig = current_contig.split(">")[1]
                        genome_fasta_dict[current_contig] = ''
                    except KeyError:
                        raise KeyError("Some contig names contain redundant names when seperated by white space")
                else:
                    line = line.split('\n')[0]
                    genome_fasta_dict[current_contig] = genome_fasta_dict[current_contig] + line

            # Check if FASTA part of gff file has been found,
            # if then indicate to start recording sequences,
            # else compare locus_tag of the line
            elif '##FASTA' in line:
                found_fasta = True

            else:
                # Compare the locus_tag with the previously largest locus_tag,
                # save the largest of the two
                if '#' not in line:
                    line = line.split('\t')
                    line = line[8].split(';')
                    line = [element for element in line if 'locus_tag' in element]
                    # Examine non empty locus_tags
                    if len(line) > 0:
                        line = line[0]
                        line = line.split('locus_tag=')[1]
                        line = line.strip()
                        if line > largest_locus_tag:
                            largest_locus_tag = line
                # Save header lines
                else:
                    header_lines.append(line)

    return genome_fasta_dict, largest_locus_tag, header_lines


def add_gene_to_gff(tmp_gff, gene_oi, genome_oi, contig, strand, refound_gene_tag, annotation, largest_locus_tag):
    """
    Function to construct and append a line to a file.
    :param tmp_gff: An open file to append a line to
    :param gene_oi: Gene in question
    :param genome_oi: Genome in question
    :param contig: Contig of the genome
    :param strand: Strand of the gene in question
    :param annotation: Any annotation found in Panaroo
    :param refound_gene_tag: Tag given by Panaroo
    :param largest_locus_tag: The current largest locus_tag
    :return: The new largest locus_tag
    """
    gene_start = genome_oi.find(gene_oi) + 1
    gene_end = gene_start + len(gene_oi) - 1
    locus_tag_parts = largest_locus_tag.rsplit('_', maxsplit=1)
    tag_length = len(locus_tag_parts[1])

    locus_tag_parts[1] = int(locus_tag_parts[1]) + 1
    preceding_zeros = str(0) * (tag_length - len(str(locus_tag_parts[1])))
    # Add preceding zeros
    locus_tag_parts[1] = f'{preceding_zeros}{locus_tag_parts[1]}'
    new_locus_tag = f'{locus_tag_parts[0]}_{locus_tag_parts[1]}'
    description = annotation[1]
    name = annotation[0]

    # Construct gff field 9
    info_field = f'ID={new_locus_tag};locus_tag={new_locus_tag};old_locus_tag={refound_gene_tag}'

    if name != '':
        info_field += f';name={name}'
    if description != '':
        info_field += f';annotation={description}'

    # construct tab delimited string containing the features of the gene
    gff_line = f'{contig}\t' \
               f'Panaroo\t' \
               f'CDS\t' \
               f'{gene_start}\t' \
               f'{gene_end}\t' \
               f'.\t' \
               f'{strand}\t' \
               f'0\t' \
               f'{info_field}'

    tmp_gff.write(gff_line + '\n')

    return new_locus_tag


def write_contig(file, contig_name, sequence):
    """
    Write contig into file
    :param file: Open file for appending contig to
    :param contig_name: Name of the contig to be added
    :param sequence: Sequence of contig to be added
    :return: Nothing
    """
    # Write contig name
    file.write(f'>{contig_name}\n')

    # Write bulk of sequence
    for i in range(len(sequence) // 60):
        file.write(sequence[0+60*i:60+60*i] + '\n')

    # Write remainder of sequence
    remainder = len(sequence) % 60
    genome_length = len(sequence)
    file.write(sequence[len(sequence) - remainder:genome_length+1] + '\n')


def annotate_refound_genes(gff_name, gene_data_dict, tmp_folder_path, corrected_gff_out_dir, logger):
    """
    Function to add back in genes that are refound by Panaroo into gff files.
    :param gff_name: File path of gff to be corrected
    :param gene_data_dict: Dict of refound genes identified from gene_presence_absence.csv file
    :param tmp_folder_path: File path to the temporary folder
    :param corrected_gff_out_dir: File path to the folder where corrected genomes should be place
    :param logger: Program logger
    :return: Nothing.
    """
    """ Function to annotate the genes refound by Panaroo in a gff3 file"""
    # Read in a gff file
    # Get base name of gff file and construct path to database in temporary folder
    gff_file_name = os.path.basename(gff_name)
    data_base = os.path.join(tmp_folder_path, f'{gff_file_name}_db')

    # Create a database for the gff file
    gffutils.create_db(gff_name, data_base)
    # Attach database
    gff_db = gffutils.FeatureDB(data_base)

    # Write quick tmp database for appending new genes.
    tmp_gff = os.path.join(tmp_folder_path, f'{gff_file_name.split(".gff")[0]}_tmp.gff')
    with open(tmp_gff, 'w') as gff_file:
        for feature in gff_db.all_features():
            gff_file.writelines(str(feature) + '\n')

    # Pass the gff file manually to extract the genome fasta sequence(s) and the largest locus_tag
    fasta_genome, largest_locus_tag, header_lines = extract_genome_fasta(gff_name)

    # Find all refound genes for given genome in gene data file
    genome_name = gff_file_name.split('.')[0]

    # Search for the refound genes and record their coordinate, strand and add them to the gff file
    with open(tmp_gff, 'a') as tmp_gff_file:
        for refound_gene in gene_data_dict[genome_name]: # .keys()
            gene_oi = gene_data_dict[genome_name][refound_gene][0]

            strand = None
            contig_counter = 0
            contigs = list(fasta_genome) # .keys()
            while strand is None and contig_counter < len(contigs) :
                contig = contigs[contig_counter]
                genome_oi = fasta_genome[contig]

                if gene_oi in genome_oi:
                    strand = '+'

                else:
                    # get reverse complement of the gene
                    gene_oi = Seq.reverse_complement(gene_oi)
                    if gene_oi in genome_oi:
                        strand = '-'

                contig_counter += 1

            if strand is not None:
                # Add the gene to the gff file.
                largest_locus_tag = add_gene_to_gff(tmp_gff_file, gene_oi, genome_oi, contig, strand,
                                                    refound_gene, gene_data_dict[genome_name][refound_gene][1:], largest_locus_tag)
            else:
                exit_with_error(f"When correcting gff {gff_name}, the gene: {refound_gene} "
                                 f"did not have any hit in the genome!", EXIT_GFF_REANNOTATION_ERROR, logger)

    # Construct a database from the temporary gff that contain the added annotations
    path_tmp_gff_db = os.path.join(tmp_folder_path, f'{gff_file_name}_tmp_db')
    # make database
    gffutils.create_db(tmp_gff, path_tmp_gff_db)
    # Attach database
    tmp_gff_db = gffutils.FeatureDB(path_tmp_gff_db)

    # Print the final GFF3 file
    corrected_gff_file = os.path.join(corrected_gff_out_dir, f'{gff_file_name.split(".gff")[0]}_corrected.gff')
    with open(corrected_gff_file, 'w') as gff_file:
        # Write initial lines
        for line in header_lines:
            gff_file.write(line)

        # ADD the gff lines
        for feature in tmp_gff_db.all_features(order_by=('seqid', 'start')):
            gff_file.writelines(str(feature) + '\n')

        # Write line to separate genome fasta
        gff_file.write("##FASTA\n")

        # Write the genome fasta
        for contig_name in fasta_genome.keys():
            write_contig(gff_file, contig_name, fasta_genome[contig_name])

    # remove database for gff file and the temporary gff file in temporary folder
    os.remove(data_base)
    os.remove(tmp_gff)
    os.remove(path_tmp_gff_db)

    return corrected_gff_file


if __name__ == '__main__':
    pass
    # _, _, attribute_dict = read_gene_presence_absence('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_pan_split_paralogs/gene_presence_absence_roary.csv',
    #                                                   1, 0.05)
    #
    # correct_gffs(['/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_genomes/GCA_008694005.gff'], '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_pan_split_paralogs/gene_data.csv',
    #                  "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests", attribute_dict)
    # # genome_dict = extract_genome_fasta('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_genomes/GCA_000006785.gff')
