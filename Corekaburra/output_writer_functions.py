import csv
import os
import time


def master_info_writer(master_info, out_path, prefix, quiet):
    if not quiet:
        print("Printing master output")

    # Write general content
    out_file_name = 'low_frequency_gene_placement.tsv'
    if prefix is not None:
        out_file_name = prefix + '_' + out_file_name
    with open(os.path.join(out_path, out_file_name), 'w', newline='', encoding='utf-8') as out_file:
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
    out_file_name = 'core_core_accessory_gene_content.tsv'
    if prefix is not None:
        out_file_name = prefix + '_' + out_file_name

    with open(os.path.join(out_path, out_file_name), 'w', newline='', encoding='utf-8') as out_file:
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


def summary_info_writer(master_summary_info, out_path, prefix, quiet):
    if not quiet:
        print("Printing master output")

    # Generate file name
    out_file_name = 'core_pair_summary.csv'
    if prefix is not None:
        out_file_name = prefix + '_' + out_file_name

    # Write general content
    with open(os.path.join(out_path, out_file_name), 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file)

        # Create header
        header = ['Core_pair', 'n',
                  'min_dist', 'max_dist', 'mean_dist', 'median_dist',
                  'min_acc', 'max_acc', 'mean_acc', 'median_acc']
        writer.writerow(header)

        # Write remaining rows:
        for key in sorted(master_summary_info.keys()):
            info = master_summary_info[key]

            writer.writerow(info)
    out_file.close()

def segment_writer(segments, out_path, prefix, quiet):
    if not quiet:
        print("Printing core segments")

    # Generate file name
    out_file_name = 'core_segments.csv'
    if prefix is not None:
        out_file_name = prefix + '_' + out_file_name

    # Write general content
    with open(os.path.join(out_path, out_file_name), 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file)

        # Create header
        header = ['Segment_name', 'segment_position', 'core_gene']
        writer.writerow(header)

        # Write remaining rows:
        for key in sorted(segments.keys()):
            for index, gene in enumerate(segments[key]):
                info = [key, index+1, gene]

                writer.writerow(info)
    out_file.close()


def no_acc_segment_writer(no_acc_segments, out_path, prefix, quiet):
    if not quiet:
        print("Printing core segments without accessory content")

    # Generate file name
    out_file_name = 'no_accessory_core_segments.csv'
    if prefix is not None:
        out_file_name = prefix + '_' + out_file_name

    # Write general content
    with open(os.path.join(out_path, out_file_name), 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file)

        # Create header
        header = ['Parent_Segment_name', 'Sub_segment_name', 'Parent_segment_position', 'Sub_segment_position', 'core_gene']
        writer.writerow(header)

        # Write remaining rows:
        for key in sorted(no_acc_segments.keys()):
            for sub_index, subsegment in enumerate(no_acc_segments[key]):
                sub_name = f'{subsegment[0]}--{subsegment[-1]}'
                for index, gene in enumerate(subsegment):
                    info = [key, sub_name, sub_index + 1, index + 1, gene]

                    writer.writerow(info)
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


