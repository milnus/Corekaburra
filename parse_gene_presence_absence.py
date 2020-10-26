import os
import csv
from math import ceil, floor


def read_gene_presence_absence(file_name, core_gene_presence, low_freq_gene, input_path=""):
    """Function that pass a Roary gene presence/absence file and returns directories of core and low frequency genes"""

    file = os.path.join(input_path, file_name)

    if os.path.isfile(file):
        with open(file, 'r', newline='', ) as gene_presence_absence:
            # Read column header line
            gff_files = gene_presence_absence.readline()
            # Strip for whitespace
            gff_files = gff_files.strip()
            # split column names
            gff_files = gff_files.split(',')

            # Index genomes and column position in dict for better search
            gff_file_dict = {}
            for i, file_name in enumerate(gff_files[14:]):
                gff_file_dict[file_name] = i

            # Read remaining lines and construct a nested dicts one dict for each genome and its core genes,
            # and a dict for low frequency genes found in less than set percent of isolates

            # Initialise reader object
            reader = csv.reader(gene_presence_absence, delimiter=',')
            # Counters
            core_gene_number = 1
            low_freq_gene_number = 1
            # Determine number of isolates that represent core and low frequency genes
            core_gene_isolate_presence = floor(len(gff_file_dict.keys()) * core_gene_presence)
            low_freq_gene_isolate_presence = ceil(len(gff_file_dict.keys()) * low_freq_gene)

            # initialise data structure to return
            core_gene_dict = dict.fromkeys(gff_files[14:])
            low_freq_gene_dict = dict.fromkeys(gff_files[14:])
            gene_cluster_lookup = {}
            # Set keys
            for key in core_gene_dict.keys():
                core_gene_dict[key] = {}
                low_freq_gene_dict[key] = {}

            # Read lines from file and determine if core, low frequency or 'regular' accessory.
            for line in reader:
                gene_isolate_presence = int(line[3])
                avg_gene_presence = int(line[4])

                # Check if core gene, if then add annotations to genomes
                if core_gene_isolate_presence <= gene_isolate_presence == avg_gene_presence:
                    for genome in core_gene_dict.keys():
                        core_gene_dict[genome][line[14+gff_file_dict[genome]]] = line[0]
                    core_gene_number += 1

                if low_freq_gene_isolate_presence >= gene_isolate_presence == avg_gene_presence:
                    for genome in low_freq_gene_dict.keys():
                        if len(line[14+gff_file_dict[genome]]) > 0:
                            low_freq_gene_dict[genome][line[14+gff_file_dict[genome]]] = line[0]
                    low_freq_gene_number += 1

    return core_gene_dict, low_freq_gene_dict
