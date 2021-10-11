import os
import csv
from math import ceil, floor
import gffutils


def check_fragmented_gene(fragments_in_line, input_gffs, temp_folder_path):
    return_list = []
    for fragment in fragments_in_line:
        # TODO - write unit test for function.
        # split the two fragments
        fragments = fragment.split(';')

        # Get the name of the genome
        genome = fragments[0].rsplit("_", 1)[0]

        # Get the gff and its path
        try:
            gff_file = [file for file in input_gffs if f'{genome}.gff' in file][0]
        except IndexError:
            raise NotImplementedError(f'No gff match was found when searching fragments for genome: {genome}')

        # Construct gff database to be searched
        db_name = os.path.join(temp_folder_path, f'{genome}_db')
        if not os.path.isfile(db_name):
            gffutils.create_db(gff_file, db_name)

        # Attach database
        gff_database = gffutils.FeatureDB(db_name)

        # Check that all fragments are on the same contig.
        first_fragment_contig = gff_database[fragments[0]][0]
        frag_same_contig = all([first_fragment_contig == gff_database[fragment][0] for fragment in fragments])
        if frag_same_contig:
            # TODO - Get the coordinate of the fragments
            # Get all coordinates
            frag_coors = []
            for frag in fragments:
                frag_coors.append(gff_database[frag][3])
                frag_coors.append(gff_database[frag][4])

            # Construct region to be searched for annotations between fragments:
            max_frag_coor = max(frag_coors)
            min_frag_coor = min(frag_coors)
            region = (first_fragment_contig, min_frag_coor, max_frag_coor)

            # Find all features that are completly within the region
            region_features = gff_database.region(region=region, completely_within=True)
            # print(list(region_features))

            # find all genes that are not part of the fragmented gene
            region_locus_tags = set([feature[8]['locus_tag'][0] for feature in region_features])
            # print(region_locus_tags)
            # print(fragments)
            excess_genes = region_locus_tags.difference(fragments)


            # check the number of excess genes, if any then False to being core
            if len(excess_genes) > 0:
                return_list.append(False)
            else:
                return_list.append(True)

    return return_list

    # TODO - Find out how the gff parser handles this? Does there need to be a check if a gene cluster is being paired to it self and if then drop it and change the end coordinates.


def read_gene_presence_absence(file_name, core_gene_presence, low_freq_gene, source_program, input_gffs, temp_folder_path, verbose=True):
    """Function that pass a Roary style gene presence/absence file.
    Returns directories of core and low frequency genes, and a directory of pan genome clusters and their annotation"""

    file = os.path.join("", file_name)

    # Check if file exists, if the read otherwise raise error.
    if os.path.isfile(file):  # TODO - THIS CHECK IS IRRELEVANT WITH THE CHECK OF THE SOURCE PROGRAM!
        with open(file, 'r', newline='', ) as gene_presence_absence:
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

            if verbose:
                print(f"\n------------Opening the gene presence/absence file------------\n")
                print(f"Core genes must be found in {core_gene_isolate_presence} or more isolates")
                print(f"Low frequency genes must be found in {low_freq_gene_isolate_presence} or fewer isolates\n")

            # initialise dict of dicts to hold genes from each gffs and to be returned
            core_gene_dict = {item: {} for item in gff_file_names[14:]}
            low_freq_gene_dict = {item: {} for item in gff_file_names[14:]}
            acc_gene_dict = {item: {} for item in gff_file_names[14:]}

            # Initialise dict that contain annotations
            annotation_dict = {}

            # Read lines from file and determine if core, low frequency or 'regular' accessory and record annotations
            for line in reader:
                # Remove quotes if Roary
                if source_program == 'Roary':
                    line = [element.replace('"', '') for element in line]

                # Record annotations of refound genes
                if any(['refound' in gene for gene in line[14:]]):
                    refound_genes = [gene for gene in line[14:] if 'refound' in gene]
                    for gene in refound_genes:
                        annotation_dict[gene] = line[2]

                # Get number of genes in line and average presence of genes in genomes
                gene_isolate_presence = int(line[3])
                avg_gene_presence = int(line[4])

                # Check if core gene, if then add annotations to genomes
                # TODO - Handle genes that have a paralog and are concatenated by ';', and check if neighbours
                # Check if gene is present in all genomes and no one gene is fragmented
                if core_gene_isolate_presence <= gene_isolate_presence == avg_gene_presence:
                    # Add gene cluster to genomes
                    for genome in core_gene_dict.keys():
                        # Check if there is an annotation for the given genome
                        if len(line[14 + gff_file_dict[genome]]) > 0:
                            core_gene_dict[genome][line[14+gff_file_dict[genome]]] = line[0]
                    core_gene_number += 1

                # Check if gene is present in all genomes, but more than one copy is pressent
                elif core_gene_isolate_presence <= gene_isolate_presence:
                    # Identify annotations for genomes that are fragmented genes
                    fragments_in_line = [genes for genes in line[14:] if ';' in genes]

                    # Check that each annotation is neighboring the other annotation.
                    return_list = check_fragmented_gene(fragments_in_line, input_gffs, temp_folder_path)

                    # Check if gene was found to be a core gene
                    if all(return_list):
                        # Add the gene to the annotation dict
                        for genome in core_gene_dict.keys():
                            # Get the annoations for a specific genome
                            genes_in_genome = line[14 + gff_file_dict[genome]]
                            # If there is an annotation add id
                            if len(genes_in_genome) > 0:
                                # Check if genome has fragments of genes,
                                # if then add them all to the annotation dict,
                                # if not then just ad the single annotation
                                if ';' in genes_in_genome:
                                    for gene in genes_in_genome.split(';'):
                                        core_gene_dict[genome][gene] = line[0]
                                else:
                                    core_gene_dict[genome][genes_in_genome] = line[0]
                        core_gene_number += 1

                    else:
                        # Check if low frequency, if then add else then add as normal accessory
                        if low_freq_gene_isolate_presence >= gene_isolate_presence == avg_gene_presence:  # TODO - review this == statement, should it be there?
                            for genome in low_freq_gene_dict.keys():
                                if len(line[14 + gff_file_dict[genome]]) > 0:
                                    low_freq_gene_dict[genome][line[14 + gff_file_dict[genome]]] = line[0]
                            low_freq_gene_number += 1
                        else:
                            for genome in acc_gene_dict.keys():
                                if len(line[14 + gff_file_dict[genome]]) > 0:
                                    acc_gene_dict[genome][line[14 + gff_file_dict[genome]]] = line[0]
                            acc_gene_number += 1

                # Check if accessory if then add annotation to genomes
                elif low_freq_gene_isolate_presence >= gene_isolate_presence == avg_gene_presence: # TODO - review this == statement, should it be there?
                    for genome in low_freq_gene_dict.keys():
                        if len(line[14+gff_file_dict[genome]]) > 0:
                            low_freq_gene_dict[genome][line[14+gff_file_dict[genome]]] = line[0]
                    low_freq_gene_number += 1

                # If not core or low frequency count as regular accessory
                else:
                    for genome in acc_gene_dict.keys():
                        if len(line[14+gff_file_dict[genome]]) > 0:
                            acc_gene_dict[genome][line[14+gff_file_dict[genome]]] = line[0]
                    acc_gene_number += 1

        if verbose:
            print("A total of:")
            print(f"{core_gene_number} core gene clusters were identified")
            print(f"{low_freq_gene_number} low frequency gene clusters were identified")
            print(f"{acc_gene_number} intermediate accessory gene clusters were identified\n")
    else:
        raise FileNotFoundError('Given gene presence absence file not found. Please check and try again.')

    # Remove gff databases
    files_in_tmp = os.listdir(temp_folder_path)
    gff_dbs = [file for file in files_in_tmp if '_db' in file]
    [os.remove(os.path.join(temp_folder_path, db)) for db in gff_dbs]

    return core_gene_dict, low_freq_gene_dict, acc_gene_dict, annotation_dict


if __name__ == '__main__':
    core_dict, *_ = read_gene_presence_absence('/Users/mjespersen/Downloads/gene_presence_absence.csv', 1, 0.05, 'Roary')
