import csv
import time

def master_info_writer(master_info, verbose=False):
    if verbose:
        print("Printing master output")


    with open('low_frequency_gene_placement.tsv', 'w', newline='', encoding='utf-8', ) as out_file:
        writer = csv.writer(out_file, delimiter="\t")

        # Create header
        header = ['Gff', 'Core_gene_1', 'Core_gene_2', 'Core_region_size',
                  'Core_region_accessory_count', 'Core_region_low_frequency_count']
        writer.writerow(header)

        # Write remaining rows:
        for key in master_info.keys():
            row_info = master_info[key][0:3] + master_info[key][3] + [master_info[key][4]] + [len(master_info[key][5])]
            writer.writerow(row_info)
    out_file.close()
