from Bio import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
from gffutils.gffwriter import GFFWriter
import os
import concurrent.futures
# from Bio.Blast import
from time import time

# REMOVE!
from parse_gene_presence_absence import read_gene_presence_absence


def read_gene_data(gene_data_file):
    """ Function to read the gene_data.csv file outputted by Panaroo and
    return a dict of genomes with their refound genes"""
    # Construct dictionary to hold refound genes and sequences for these
    gene_data_dict = {}

    # Read the gene_data.csv file and record all refound genes
    with open(gene_data_file, 'r') as gene_data:

        for line in gene_data.readlines():
            # Split read line at commas
            line = line.split(',')

            # Check if refound gene
            if 'refound' in line[2]:
                # Try to add the refound gene to the gene_data dict as a key to value being the DNA sequence,
                # if the genome is not found in gene_data dict, then construct dict for the genome and add the gene
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
                    # Examien non empty locus_tags
                    if len(line) > 0:
                        line = line[0]
                        line = line.split('locus_tag=')[1]
                        if line > largest_locus_tag:
                            largest_locus_tag = line
                # Save header lines
                else:
                    header_lines.append(line)

    return genome_fasta_dict, largest_locus_tag, header_lines


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

    return gff_line, new_locus_tag


def add_gene_to_gff(tmp_gff, gene_oi, genome_oi, contig, strand, annotation, refound_gene_tag, largest_locus_tag):
    gene_start = genome_oi.find(gene_oi) + 1
    gene_end = gene_start + len(gene_oi)
    locus_tag_parts = largest_locus_tag.rsplit('_', maxsplit=1)
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

    with open(tmp_gff, 'a') as tmp_gff_file:
        tmp_gff_file.write(gff_line + '\n')
    tmp_gff_file.close()

    return new_locus_tag


def write_contig(file, contig_name, sequnce):
    # Write contig name
    file.write(f'>{contig_name}\n')

    # Write bulk of sequence
    for i in range(len(sequnce) // 60):
        file.write(sequnce[0+60*i:60+60*i] + '\n')

    # Write remainder of sequence
    remainder = len(sequnce) % 60
    genome_length = len(sequnce)
    file.write(sequnce[len(sequnce) - remainder:genome_length+1] + '\n')


def annotate_refound_genes(gff_name, gene_data_dict, temp_folder_path, annotation_dict, corrected_gff_out_dir, i):
    """ Function to annotate the genes refound by Panaroo in a gff3 file"""
    # Print info on progress
    if (i+1) % 25 == 0 or i == 0:
        print(f"Correcting GFF file #{i+1}")

    # Read in a gff file
    # Get base name of gff file and construct path to database in temporary folder
    gff_file_name = os.path.basename(gff_name)
    data_base = os.path.join(temp_folder_path, f'{gff_file_name}_db')

    # Create a database for the gff file
    gffutils.create_db(gff_name, data_base)
    # Attach database
    gff_db = gffutils.FeatureDB(data_base)

    # Write quick tmp database for appending new genes.
    tmp_gff = os.path.join(temp_folder_path, f'{gff_file_name.split(".gff")[0]}_tmp.gff')
    with open(tmp_gff, 'w') as gff_file:
        for feature in gff_db.all_features():
            gff_file.writelines(str(feature) + '\n')
    gff_file.close()

    # Pass the gff file manually to extract the genome fasta sequence(s) and the largest locus_tag
    fasta_genome, largest_locus_tag, header_lines = extract_genome_fasta(gff_name)

    # Find all refound genes for given genome in gene data file
    genome_name = gff_file_name.split('.')[0]

    # Search for the refound genes and record their coordinate, strand,
    for refound_gene in gene_data_dict[genome_name].keys():
        gene_oi = gene_data_dict[genome_name][refound_gene]

        strand = None
        contig_counter = 0
        contigs = list(fasta_genome.keys())
        while strand is None:
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
            largest_locus_tag = add_gene_to_gff(tmp_gff, gene_oi, genome_oi, contig, strand,
                                   annotation_dict[refound_gene], refound_gene, largest_locus_tag)

        else:
            raise ValueError(f"When correcting gff {gff_name}, the gene: {refound_gene} "
                             f"did not have any hit in the genome!")

    # Construct a database from the temporary gff that contain the added annotations
    path_tmp_gff_db = os.path.join(temp_folder_path, f'{gff_file_name}_tmp_db')
    # make database
    gffutils.create_db(tmp_gff, path_tmp_gff_db)
    # Attach database
    tmp_gff_db = gffutils.FeatureDB(path_tmp_gff_db)

    # Print the final GFF3 file
    with open(os.path.join(corrected_gff_out_dir, f'{gff_file_name.split(".gff")[0]}_corrected.gff'), 'w') as gff_file:
        # Write initial lines
        for line in header_lines:
            gff_file.write(line)

        # ADD the gff lines
        for feature in tmp_gff_db.all_features(order_by=('seqid', 'start')):
            gff_file.writelines(str(feature) + '\n')

        # Write line to seperate genome fasta
        gff_file.write("##FASTA\n")

        # Write the genome fasta
        for contig_name in fasta_genome.keys():
            write_contig(gff_file, contig_name, fasta_genome[contig_name])

    # remove database for gff file and the temporary gff file in temporary folder
    os.remove(data_base)
    os.remove(tmp_gff)
    os.remove(path_tmp_gff_db)


def correct_gffs(gffs, gene_data_file, output_folder, annotation_dict, temp_folder_path):
    # TODO - Make verbose
    print("Reading gene_data.csv")
    # Read Gene_data.csv file into dict with a dict of refound genes for each genome
    gene_data_dict = read_gene_data(gene_data_file)

    # Construct directory to hold corrected gff files:
    corrected_gff_out_dir = os.path.join(output_folder, 'Corrected_gff_files')
    # Try and construct folder, if present check if content matches input to avoid process
    try:
        os.mkdir(corrected_gff_out_dir)
    except FileExistsError:
        corrected_folder_content = os.listdir(corrected_gff_out_dir)

        gff_names = [os.path.basename(gff) for gff in gffs]

        corrected_files = [file for file in corrected_folder_content if f'{file.split("_corrected")[0]}.gff' in gff_names]
        if len(corrected_files) == len(gffs):
            os.rmdir(temp_folder_path)
            return corrected_gff_out_dir

        else:
            os.rmdir('genome_corer_tmp')
            raise FileNotFoundError
            #TODO - RAISE ERROR THAT SOME FILES ARE NOT IN THE corrected
            # Find the files and process them
            # If none of the files match then abort.


    # TODO - Make verbose
    print("Start correcting GFF files with refound genes from Panaroo")
    # Multi process the annotation of the genomes. (Process has been found to be better than threads)
    total_time = time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=15) as executor:
        result = [executor.submit(annotate_refound_genes, gff, gene_data_dict, temp_folder_path,
                                  annotation_dict, corrected_gff_out_dir, i) for i, gff in enumerate(gffs)]

        for f in concurrent.futures.as_completed(result):
            f.result()


    print(f"total time for correction {time() - total_time}")

    return corrected_gff_out_dir


if __name__ == '__main__':
    _, _, attribute_dict = read_gene_presence_absence('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_pan_split_paralogs/gene_presence_absence_roary.csv',
                                                      1, 0.05)

    correct_gffs(['/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_genomes/GCA_008694005.gff'], '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_pan_split_paralogs/gene_data.csv',
                     "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests", attribute_dict)
    # genome_dict = extract_genome_fasta('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/50_refseq_genomes/GCA_000006785.gff')
