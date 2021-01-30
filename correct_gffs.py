from Bio import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
import os
import concurrent.futures
# from Bio.Blast import

# REMOVE!
from parse_gene_presence_absence import read_gene_presence_absence


def read_gene_data(gene_data_file):
    """ Function to read the gene_data.csv file outputted by Panaroo """
    # Construct dictionary to hold refound genes and sequences for these
    gene_data_dict = {}

    # Read the gene_data.csv file and record all refound genes
    with open(gene_data_file, 'r') as gene_data:

        for line in gene_data.readlines():
            # Split read line at commas
            line = line.split(',')

            # Check if refound gene
            if 'refound' in line[2]:
                # Try to add the refound gene to the gene_data dict,
                # if the genome is not found gene gene_data dict, then construct dict for the genome and add the gene
                try:
                    gene_data_dict[line[0]][line[2]] = line[5]
                except KeyError:
                    gene_data_dict[line[0]] = {}
                    gene_data_dict[line[0]][line[2]] = line[5]

                # TODO - Find example of gene with differnt length and a premature stop codon.
                # if 'len' in line[2]:
                #     print(line)
                # if 'stop' in line[2]:
                #     print(line)
    gene_data.close()

    return gene_data_dict


def extract_genome_fasta(gff_name):
    """ Function to read a gff3 file and extract the fasta seuqnce at the end
    each contig is returned as an entry in a dictionary """
    # Initialise the two return variables a dict for refound gene's annotaton and the largest locus_tag
    genome_fasta_dict = {}
    largest_locus_tag = ''

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
                    # Examien non empty locus_tags
                    if len(line) > 0:
                        line = line[0]
                        line = line.split('locus_tag=')[1]
                        if line > largest_locus_tag:
                            largest_locus_tag = line

    return genome_fasta_dict, largest_locus_tag


def construct_gff_line(gene_oi, genome_oi, contig, strand, annotation, refound_gene_tag, largest_locus_tag):
    # Prepare the start, end, and locus_tag of the gene feature line.
    gene_start = genome_oi.find(gene_oi) + 1
    gene_end = gene_start + len(gene_oi)
    locus_tag_parts = largest_locus_tag.split('_')
    locus_tag_parts[1] = int(locus_tag_parts[1]) + 1
    new_locus_tag = f'{locus_tag_parts[0]}_{locus_tag_parts[1]}'

    # construct tab delimited string containing the features of the gene
    gff_line = f'{contig}\t' \
               f'Panaroo\t' \
               f'CDS\t' \
               f'{gene_start}\t' \
               f'{gene_end}\t' \
               f'.\t' \
               f'{strand}\t' \
               f'0\t' \
               f'ID={new_locus_tag};annotation={annotation};locus_tag={new_locus_tag};old_locus_tag={refound_gene_tag}'

    # Make gff line into gff feature line
    gene_feature_line = gffutils.feature.feature_from_line(gff_line) # TODO: add the feature line to the gff database.

    return gene_feature_line, new_locus_tag


def annotate_refound_genes(gff_name, gene_data_dict, gene_pres_abs, temp_folder_path, annotation_dict):
    """ Function to annotate the genes refound by Panaroo in a gff3 file"""
    # Read in a gff file
    # Get base name of gff file and construct path to database in temporary folder
    gff_file_name = os.path.basename(gff_name)
    data_base = os.path.join(f'{temp_folder_path}{gff_file_name}_db')

    # Create a database for the gff file
    gffutils.create_db(gff_name, data_base)
    # Attach database
    gff_db = gffutils.FeatureDB(data_base)

    # Pass the gff file manually to extract the genome fasta sequence(s) and the largest locus_tag
    fasta_genome, largest_locus_tag = extract_genome_fasta(gff_name)

    # Find all refound genes for given genome in gene data file
    genome_name = gff_file_name.split('.')[0]

    # Search for the refound genes and record their coordinate, strand,
    for refound_gene in gene_data_dict[genome_name].keys():
        gene_oi = gene_data_dict[genome_name][refound_gene]

        for contig in fasta_genome.keys():
            genome_oi = fasta_genome[contig]
            strand = None

            if gene_oi in genome_oi:
                print("Found in contig")
                strand = '+'

            else:
                # get reverse complement of the gene
                gene_oi = Seq.reverse_complement(gene_oi)
                if gene_oi in genome_oi:
                    print("Reversed - Found in contig")
                    strand = '-'

            if strand is not None:
                gff_line, largest_locus_tag = construct_gff_line(gene_oi, genome_oi, contig, strand,
                                   annotation_dict[refound_gene], refound_gene, largest_locus_tag)
                print(gff_line)
                # gff_db = gff_db.update(data=gff_line, make_backup=False)
                gff_db = gffutils.FeatureDB.update(self=gff_db, data=gff_line)

            else:
                raise ValueError(f"When correcting gff {gff_name}, a gene did not have any hit in the genome!")

    # TODO - Add the annotation of the refound gene to the GFF file for the given genome. with contig, coordinates, annotation, locus_tag as .5 and old_locus_tag being the refound identifyer.
        # TODO - if no contig is found then it's a complete genome. Find what the genome is called.

    # TODO - Write the new GFF file.

    # remove database for gff file in temporary folder
    os.remove(data_base)


def correct_gffs(gffs, gene_data_file, gene_pres_abs, output_folder, annotation_dict):
    # Read Gene_data.csv file into dict with a dict of refound genes for each genome
    gene_data_dict = read_gene_data(gene_data_file)

    # Construct temporary folder with write permission to store gff databases
    temp_folder_path = os.path.join(output_folder, 'genome_corer_tmp')
    os.mkdir(temp_folder_path)

    # TODO - multi process or thread commands below.
    annotate_refound_genes('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_genomes/GCA_000006785.gff',
        gene_data_dict, None,
        '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/', annotation_dict)

    # Remove temporary database holding gff databases
    os.rmdir(temp_folder_path)


if __name__ == '__main__':
    _, _, attribute_dict = read_gene_presence_absence('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_pan_split_paralogs/gene_presence_absence_roary.csv',
                                                      1, 0.05)

    correct_gffs(None, '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_pan_split_paralogs/gene_data.csv',
                     None, "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests", attribute_dict)
    # genome_dict = extract_genome_fasta('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_genomes/GCA_000006785.gff')
    # print(genome_dict)
