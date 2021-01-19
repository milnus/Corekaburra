from Bio import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
import os
import concurrent.futures
# from Bio.Blast import


def read_gene_data(gene_data_file):
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


def annotate_refounds(gff_name, gene_data_dict, gene_pres_abs, temp_folder_path):
    # TODO - Read in a gff file
    # Get base name of gff file and construct path to database in temporary folder
    gff_file_name = os.path.basename(gff_name)
    data_base = os.path.join(f'{temp_folder_path}{gff_file_name}_db')

    # Create a database for the gff file
    gffutils.create_db(gff_name, data_base)

    # TODO - Find all refound genes for given genome in gene data file

    # TODO - Identify the pan genome cluster for each refound gene

    # TODO - Get the annotation for each refound gene

    # TODO - BLAST refound gene against the given genome to determine its position in the genome

    # TODO - Add the annotation of the refound gene to the GFF file for the given genome. with contig, coordinates, annotation, locus_tag as .5 and old_locus_tag being the refound identifyer.
        # TODO - if no contig is found then it's a complete genome. Find what the genome is called.

    # TODO - Write the new GFF file.

    # remove database for gff file in temporary folder
    os.remove(data_base)


def annotate_refound(gffs, gene_data_file, gene_pres_abs, output_folder):
    # Read Gene_data.csv file into dict with a dict of refound genes for each genome
    gene_data_dict = read_gene_data(gene_data_file)

    # Construct temporary folder with write permission to store gff databases
    temp_folder_path = os.path.join(output_folder, 'genome_corer_tmp')
    os.mkdir(temp_folder_path)

    # TODO - multi process or thread commands below.

    # Remove temporary database holding gff databases
    os.rmdir(temp_folder_path)

annotate_refounds('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_genomes/GCA_000006785.gff',None, None,
                  '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/')
