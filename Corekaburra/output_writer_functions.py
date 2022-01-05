import csv
import os
import time


def master_info_writer(master_info, out_path, prefix, quiet):
    """
    Function to write two output .tsv files related to regions content and size for each genome
    :param master_info: Dict of info for each core gene pair across all genomes
    :param out_path: Path to the output folder
    :param prefix: A possible prefix for the output files.
    :param quiet:
    :return: Nothing
    """
    if not quiet:
        print("Printing master output")

    # Write general content
    out_file_name = 'low_frequency_gene_placement.tsv'
    if prefix is not None:
        out_file_name = f'{prefix}_{out_file_name}'
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

    # Write gene content in long format
    out_file_name = 'core_core_accessory_gene_content.tsv'
    if prefix is not None:
        out_file_name = f'{prefix}_{out_file_name}'

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


def summary_info_writer(master_summary_info, out_path, prefix, quiet):
    """
    Function for writing the summary table for regions identified across genomes
    :param master_summary_info: Dict holding summary statistics for core pair region identified
    :param out_path: Path to the output folder
    :param prefix: Prefix for any output files
    :param quiet: # TODO - log instead
    :return: Nothing
    """
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
                  'occurrence_core_1', 'occurrence_core_2', 'co_occurrence',
                  'min_dist', 'max_dist', 'mean_dist', 'median_dist',
                  'min_acc', 'max_acc', 'mean_acc', 'median_acc']
        writer.writerow(header)

        # Write remaining rows:
        for key in sorted(master_summary_info.keys()):
            info = master_summary_info[key]

            writer.writerow(info)


def segment_writer(segments, out_path, prefix, quiet):
    """
    Function to write segments of core genes identified across the pan-genome
    :param segments: Dict of segments (lists) in values, under name of segments as keys.
    :param out_path: Path to output folder
    :param prefix: Prefix for any output files
    :param quiet: # TODO - logger
    :return: Nothing
    """
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
        header = ['Segment_name', 'Segment_position', 'Core_gene']
        writer.writerow(header)

        # Write remaining rows:
        for key in sorted(segments.keys()):
            for index, gene in enumerate(segments[key]):
                info = [key.replace('--', '-'), index+1, gene]

                writer.writerow(info)


def no_acc_segment_writer(no_acc_segments, out_path, prefix, quiet):
    """
    Function for writing segments of core genes with no accessory between them.
    :param no_acc_segments: Dict of segments with (lists) in values with sub-lists being segments with no accessory genes between them, under name of segments as keys.
    :param out_path: Path to output folder
    :param prefix: Prefix for any output files
    :param quiet: # TODO - logger
    :return: Nothing
    """
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
        header = ['Parent_segment_name', 'Sub_segment_name', 'Parent_segment_position', 'Sub_segment_position', 'Core_gene']
        writer.writerow(header)

        # Write remaining rows:
        for key in sorted(no_acc_segments.keys()):
            for sub_index, subsegment in enumerate(no_acc_segments[key]):
                sub_name = f'{subsegment[0]}-{subsegment[-1]}'
                for index, gene in enumerate(subsegment):
                    info = [key.replace('--', '-'), sub_name, sub_index + 1, index + 1, gene]

                    writer.writerow(info)


if __name__ == "__main__":
    pass