import csv
import time


def master_info_writer(master_info, verbose=False):
    if verbose:
        print("Printing master output")

    # Write general content
    with open('low_frequency_gene_placement.tsv', 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file, delimiter="\t")

        # Create header
        header = ['Gff', 'Core_gene_1', 'Core_gene_2', 'Core_region_size',
                  'Core_region_accessory_count']
        writer.writerow(header)

        # Write remaining rows:
        for key in master_info.keys():
            info = master_info[key][0:5]

            writer.writerow(info)
    out_file.close()

    # Write gene content in long format
    with open('core_core_accessory_gene_content.tsv', 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file, delimiter="\t")

        # Create header
        header = ['Gff', 'Core_gene_1', 'Core_gene_2', 'gene', 'type']
        writer.writerow(header)

        # Write remaining rows:
        for key in master_info.keys():
            core_core_region = master_info[key]
            if len(core_core_region[5]):
                for gene in core_core_region[5]:
                    row = [core_core_region[0],
                           core_core_region[1],
                           core_core_region[2],
                           gene,
                           'low_frequency']
                    writer.writerow(row)

            if len(core_core_region[6]):
                for gene in core_core_region[6]:
                    row = [core_core_region[0],
                           core_core_region[1],
                           core_core_region[2],
                           gene,
                           'intermediate_frequency']
                    writer.writerow(row)

    out_file.close()


def write_consensus_core_gene_synteny(core_gene_synteny):
    with open('consensus_core_gene_synteny.txt', 'w', newline='', encoding='utf-8') as out_file:
        for gene in core_gene_synteny:
            out_file.write(f'{gene}\n')
    out_file.close()


def write_core_gene_coverage(core_path_coverage):
    """ Function to write """
    with open('core_gene_coverage.tsv', 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file, delimiter='\t')

        header = ['Core_gene_1', 'Core_gene_2', 'Connections']
        writer.writerow(header)

        for connection in core_path_coverage:
            writer.writerow(connection)
    out_file.close()


def write_alternative_core_gene_counts(alternative_core_gene_counts):
    with open('alternative_core_pairs_count.tsv', 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file, delimiter='\t')

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


