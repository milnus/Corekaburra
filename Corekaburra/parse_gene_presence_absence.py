import os
import csv
from math import ceil, floor
import gffutils

try:
    from Corekaburra.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error


def add_gene_to_dict(main_dict, gene, pan_gene_name, genome):
    """
    Function to add a gene to a given dictionary
    :param main_dict: Dict of genes from genomes. A dict of dicts, with first set of keys being genomes, second is locus_tags with pan-genome gene being the key.
    :param gene: The gene in question from a specific genome (locus_tag)
    :param pan_gene_name: The name of the pan-genome gene (cluster) to which the above gene belongs.
    :param genome: The name of the genome in question
    :return: returns the dict to be used further
    """
    if ';' in gene:
        for gene_part in gene.split(';'):
            main_dict[genome][gene_part] = pan_gene_name
    else:
        main_dict[genome][gene] = pan_gene_name

    return main_dict


def check_fragmented_gene(fragment_info, input_gffs, tmp_folder_path, logger):
    """
    Function that check for that placement of fragmented gene parts, to determine if they are neighbouring or have some genomic feature between them
    :param fragment_info: List of genes that are found to be fragmented, one composite of fragments for each index
    :param input_gffs: A list of file-paths to the gff files given as input
    :param tmp_folder_path: A file-path to the temporary folder of the Corekaburra run
    :return: A List of booleans indicating if a fragments has nothing in between fragments (True) or not (False)
    """

    fragments_close = []
    for fragment in fragment_info:
        # split the two fragments
        fragment_pieces = fragment[0].split(';')

        # Get the name of the genome
        genome = fragment[1]

        # Get the gff and its path
        if '.gff' not in genome:
            try:
                gff_file = [file for file in input_gffs if f'{genome}.gff' in file][0]
                db_name = os.path.join(tmp_folder_path, f'{genome}_db')
            except IndexError:
                raise NotImplementedError(f'No gff match was found when searching fragments for genome: {genome}')
        else:
            gff_file = genome
            db_name = f"{os.path.basename(genome)}_db"
            db_name = os.path.join(tmp_folder_path, db_name)

        # Construct gff database to be searched
        if not os.path.isfile(db_name):
            gffutils.create_db(gff_file, db_name, force_gff=True, id_spec=['old_locus_tag', 'ID'])

        # Attach database
        gff_database = gffutils.FeatureDB(db_name)

        # Check that all fragments are on the same contig.
        first_fragment_contig = gff_database[fragment_pieces[0]][0]

        frag_same_contig = all([first_fragment_contig == gff_database[fragment][0] for fragment in fragment_pieces])
        if frag_same_contig:
            # Get all coordinates
            frag_coors = []
            for frag in fragment_pieces:
                frag_coors.append(gff_database[frag][3])
                frag_coors.append(gff_database[frag][4])

            # Construct region to be searched for annotations between fragments:
            max_frag_coor = max(frag_coors)
            min_frag_coor = min(frag_coors)
            region = (first_fragment_contig, min_frag_coor, max_frag_coor)

            # Find all features that are completely within the region
            region_features = gff_database.region(region=region, completely_within=True)

            # Find if some pieces are refound and change old_locus_tag to ID
            refound_pieces = [[i, fragment_piece] for i, fragment_piece in enumerate(fragment_pieces) if 'refound' in fragment_piece]
            if refound_pieces:
                for i, piece in refound_pieces:
                    fragment_pieces[i] = gff_database[piece]['ID'][0]
            # find all genes that are not part of the fragmented gene
            region_locus_tags = set([feature[8]['locus_tag'][0] for feature in region_features])
            excess_genes = region_locus_tags.difference(fragment_pieces)

            # check the number of excess genes, if any then False to being core
            if len(excess_genes) > 0:
                fragments_close.append(False)
            else:
                fragments_close.append(True)
        else:
            fragments_close.append(False)

    return fragments_close


def read_gene_presence_absence(pres_abs_file, core_gene_presence, low_freq_gene, source_program, input_gffs, tmp_folder_path, logger):
    """
    Function that pass a Roary style gene presence/absence file.
    :param pres_abs_file: File path to the gene presence/absence file identified
    :param core_gene_presence: The ratio of genomes in which a gene must present, to be seen as a core gene
    :param low_freq_gene: The ratio of genomes in which a gene must not surpass, to be seen as a low-frequency gene
    :param source_program: The program from which the pan-genome was produced
    :param input_gffs: A list of file-paths to the gff files given as input
    :param tmp_folder_path: A file-path to the temporary folder of the Corekaburra run
    :param logger: Program logger
    :return: Directories of directories of core and low frequency genes, and a directory of pan genome clusters and their annotation.
    """

    # Open the presence/absense file to index gene into core, accessory, or low-frequency genes
    with open(pres_abs_file, 'r', newline='', ) as gene_presence_absence:
        # Read column header line
        gff_file_names = gene_presence_absence.readline()
        # Strip for whitespace
        gff_file_names = gff_file_names.strip()
        # split column names
        gff_file_names = gff_file_names.split(',')

        # Remove the quotes from Rorary input
        if source_program == 'Roary':
            gff_file_names = [filename.replace('"', '') for filename in gff_file_names]

        # Index gff filenames and column position in dict for better search
        gff_file_dict = {}
        for i, gff_name in enumerate(gff_file_names[14:]):
            gff_file_dict[gff_name] = i

        # Read remaining lines and construct a nested dicts one dict for each genome and its core genes,
        # and a dict for low frequency genes found in less than set percent of isolates

        # Initialise reader object to read remaining lines
        reader = csv.reader(gene_presence_absence, delimiter=',')
        # Counters
        core_gene_number = 0
        low_freq_gene_number = 0
        acc_gene_number = 0

        # Determine number of isolates that represent core and low frequency genes
        core_gene_isolate_presence = floor(len(gff_file_dict.keys()) * core_gene_presence)
        low_freq_gene_isolate_presence = ceil(len(gff_file_dict.keys()) * low_freq_gene)

        logger.info(f"------------Opening the gene presence/absence file------------\n"
                    f"Core genes must be found in {core_gene_isolate_presence} or more genomes\n"
                    f"Low frequency genes must be found in less than {low_freq_gene_isolate_presence} genomes\n")

        # initialise dict of dicts to hold genes from each gffs and to be returned
        core_gene_dict = {item: {} for item in gff_file_names[14:]}
        low_freq_gene_dict = {item: {} for item in gff_file_names[14:]}
        acc_gene_dict = {item: {} for item in gff_file_names[14:]}

        # Read lines from file and determine if core, low frequency or 'regular' accessory and record annotations
        for line in reader:
            # Remove quotes if Roary
            if source_program == 'Roary':
                line = [element.replace('"', '') for element in line]

            # Get number of genes in line and average presence of genes in genomes
            gene_isolate_presence = int(line[3])
            no_seq_presence = int(line[4])

            # Check if core gene, if then add annotations to genomes
            # Check if gene is present in all genomes and no one gene is fragmented
            if core_gene_isolate_presence <= gene_isolate_presence == no_seq_presence:
                # Add gene cluster to genomes
                for genome in core_gene_dict.keys():
                    # Check if there is an annotation for the given genome
                    if len(line[14 + gff_file_dict[genome]]) > 0:
                        core_gene_dict[genome][line[14+gff_file_dict[genome]]] = line[0]
                core_gene_number += 1

            # Check if gene is present in all genomes, but more than one copy is present
            elif core_gene_isolate_presence <= gene_isolate_presence:
                # Identify annotations for genomes that are fragmented genes
                fragment_info = [[genes, gff] for genes, gff in zip(line[14:], gff_file_names[14:]) if ';' in genes]

                # Check that each annotation is neighboring the other annotation.
                fragments_close = check_fragmented_gene(fragment_info, input_gffs, tmp_folder_path, logger)
                # Check if gene was found to be a core gene
                if all(fragments_close):
                    # Add the gene to the annotation dict
                    for genome in core_gene_dict:
                        # Get the annoations for a specific genome
                        genes_in_genome = line[14 + gff_file_dict[genome]]
                        # If there is an annotation add id
                        if len(genes_in_genome) > 0:
                            # Check if genome has fragments of genes,
                            # if then add them all to the annotation dict,
                            # if not then just ad the single annotation
                            add_gene_to_dict(core_gene_dict, genes_in_genome, line[0], genome)
                    core_gene_number += 1

                else:
                    # Check if low frequency, if then add else then add as normal accessory
                    if low_freq_gene_isolate_presence >= gene_isolate_presence == no_seq_presence:
                        for genome in low_freq_gene_dict.keys():
                            if len(line[14 + gff_file_dict[genome]]) > 0:
                                add_gene_to_dict(low_freq_gene_dict, line[14 + gff_file_dict[genome]], line[0], genome)
                        low_freq_gene_number += 1
                    else:
                        for genome in acc_gene_dict.keys():
                            if len(line[14 + gff_file_dict[genome]]) > 0:
                                add_gene_to_dict(acc_gene_dict, line[14 + gff_file_dict[genome]], line[0], genome)
                        acc_gene_number += 1

            # Check if accessory if then add annotation to genomes
            elif low_freq_gene_isolate_presence >= gene_isolate_presence == no_seq_presence:
                for genome in low_freq_gene_dict.keys():
                    if len(line[14+gff_file_dict[genome]]) > 0:
                        add_gene_to_dict(low_freq_gene_dict, line[14 + gff_file_dict[genome]], line[0], genome)
                low_freq_gene_number += 1

            # If not core or low frequency count as regular accessory
            else:
                for genome in acc_gene_dict.keys():
                    if len(line[14+gff_file_dict[genome]]) > 0:
                        add_gene_to_dict(acc_gene_dict, line[14 + gff_file_dict[genome]], line[0], genome)
                acc_gene_number += 1

    logger.info("A total of:\n"
                f"{core_gene_number} core gene clusters were identified\n"
                f"{low_freq_gene_number} low frequency gene clusters were identified\n"
                f"{acc_gene_number} intermediate accessory gene clusters were identified\n")

    return core_gene_dict, low_freq_gene_dict, acc_gene_dict


if __name__ == '__main__':
    pass
