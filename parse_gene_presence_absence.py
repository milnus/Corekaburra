import os
import csv
from math import ceil, floor


def read_gene_presence_absence(file_name, core_gene_presence, low_freq_gene, verbose=True):
    """Function that pass a Roary gene presence/absence file and returns directories of core and low frequency genes"""

    file = os.path.join("", file_name)

    # Check if file exists, if the read otherwise raise error.
    if os.path.isfile(file):
        with open(file, 'r', newline='', ) as gene_presence_absence:
            # Read column header line
            gff_file_names = gene_presence_absence.readline()
            # Strip for whitespace
            gff_file_names = gff_file_names.strip()
            # split column names
            gff_file_names = gff_file_names.split(',')

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
                print(f"\n------------Opening the gene presence/absence file------------")
                print(f"Core genes must be found in {core_gene_isolate_presence} or more isolates")
                print(f"Low frequency genes must be found in {low_freq_gene_isolate_presence} or fewer isolates\n")

            # initialise dict of dicts to hold genes from each gffs and to be returned
            core_gene_dict = {item: {} for item in gff_file_names[14:]}
            low_freq_gene_dict = {item: {} for item in gff_file_names[14:]}

            # Read lines from file and determine if core, low frequency or 'regular' accessory.
            for line in reader:
                # Get number of genes in line and average presence of genes in genomes
                gene_isolate_presence = int(line[3])
                avg_gene_presence = int(line[4])

                # Check if core gene, if then add annotations to genomes
                # TODO - Handle genes that are refound
                # TODO - Handle genes that have a paralog and are concatenated by ';', and check if neighbours
                if core_gene_isolate_presence <= gene_isolate_presence == avg_gene_presence:
                    # Add gene cluster to genomes
                    for genome in core_gene_dict.keys():
                        if len(line[14 + gff_file_dict[genome]]) > 0:
                            core_gene_dict[genome][line[14+gff_file_dict[genome]]] = line[0]
                    core_gene_number += 1

                # Check if accessory if then add annotation to genomes
                elif low_freq_gene_isolate_presence >= gene_isolate_presence == avg_gene_presence:
                    for genome in low_freq_gene_dict.keys():
                        if len(line[14+gff_file_dict[genome]]) > 0:
                            low_freq_gene_dict[genome][line[14+gff_file_dict[genome]]] = line[0]
                    low_freq_gene_number += 1

                # If not core or low frequency count as regular accessory
                else:
                    acc_gene_number += 1

        if verbose:
            print("A total of:")
            print(f"{core_gene_number} core gene clusters were identified")
            print(f"{low_freq_gene_number} low frequency gene clusters were identified")
            print(f"{acc_gene_number} intermediate accessory gene clusters were identified\n")
    else:
        raise FileNotFoundError('Given gene presence absence file not found. Please check and try again.')

    return core_gene_dict, low_freq_gene_dict
