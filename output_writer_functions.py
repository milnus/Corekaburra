import csv
import time


def master_info_writer(master_info, verbose=False):
    if verbose:
        print("Printing master output")


    with open('low_frequency_gene_placement.tsv', 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file, delimiter="\t")

        # Create header
        header = ['Gff', 'Core_gene_1', 'Core_gene_2', 'Core_region_size',
                  'Core_region_accessory_count', 'Core_region_low_frequency_count']
        writer.writerow(header)

        # Write remaining rows:
        for key in master_info.keys():
            if len(master_info[key][5]) > master_info[key][4]:
                print(f"{master_info[key]} - Accessory smaller than low frequency")
            row_info = master_info[key][0:3] + \
                       [master_info[key][3]] + \
                       [master_info[key][4]] + \
                       [len(master_info[key][5])]

            writer.writerow(row_info)
    out_file.close()


def write_consensus_core_gene_synteny(core_gene_synteny, core_path_coverage):
    with open('consesus_core_gene_synteny.txt', 'w', newline='', encoding='utf-8') as out_file:
        for i, gene in enumerate(core_gene_synteny):
            for entery in list(map(lambda e: [gene, e], core_path_coverage[i])):
                out_file.write(f'{entery[0].strip()}\t{entery[1]}\n')

    out_file.close()


def write_alternative_core_gene_counts(alternative_core_gene_counts):
    with open('alternative_core_pairs_count.tsv', 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file, delimiter="\t")

        header = ['Core_gene_1', 'Core_gene_2', 'Num._connections']
        writer.writerow(header)

        for key in alternative_core_gene_counts.keys():
            split_key = key.split('--')

            row_info = [split_key[0].strip(), split_key[1].strip(), alternative_core_gene_counts[key]]

            writer.writerow(row_info)
    out_file.close()


def write_core_gene_types(core_genome_types, alt_core_pair_matrix):
    with open('core_genome_synteny_types.csv', 'w', newline='', encoding='utf-8') as out_file:
        header = ['Genome', 'Type']

        writer = csv.writer(out_file, delimiter=',')
        writer.writerow(header)
        for key in core_genome_types:
            writer.writerow([key, core_genome_types[key]])

        out_file.close()


    with open('core_pair_matrix.csv', 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.DictWriter(out_file, fieldnames=alt_core_pair_matrix[1])
        writer.writeheader()
        writer.writerows(alt_core_pair_matrix[0])


