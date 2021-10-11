# General function
from parse_gene_presence_absence import read_gene_presence_absence
from gff_parser import segment_genome_content
from correct_gffs import correct_gffs
from merge_dicts import merge_dicts_counts, merge_dicts_lists, merge_first_genes
from output_writer_functions import master_info_writer, summary_info_writer, segment_writer, no_acc_segment_writer, \
    write_consensus_core_gene_synteny, \
    write_core_gene_coverage, \
    write_alternative_core_gene_counts, \
    write_core_gene_types
from consesus_core_genome import determine_genome_segments, determine_core_gene_consesus, identify_rearrangements, \
    characterise_rearrangements, core_pair_matrix
from check_inputs import define_input_source, check_gene_data, check_gff_files, check_gene_alignments, check_gff_in_pan
from time_calculator import time_calculator
from commandline_interface import get_commandline_arguments
from read_complete_genome_file import parse_complete_genome_file
from summary_table import calculate_n_create_summaries
from construct_multi_fasta_genome import construct_consensus_alignment
# from plots import consesus_genome_coverage
import concurrent.futures
from os import listdir, mkdir, rmdir
from os.path import join, isdir
import time
import sys


def main():
    total_time_start = time.time()
    # get arguments from the commandline
    args = get_commandline_arguments(sys.argv[1:])


    # Check presence of gff files. # TODO - is this necessary?
    if check_gff_files(args.input_gffs):
        print("All .gff files were found!")


    # Parse complete genome file and check that genomes are present
    if args.comp_genomes is not None:
        comp_genomes = parse_complete_genome_file(args.comp_genomes, args.input_gffs)
    else:
        comp_genomes = None

    # Check input
    if not args.quiet:
        print("\n----Checking presence of input files in pan genome folder----\n")

    # Check if Panaroo or Roary input folder is given
    source_program, input_pres_abs_file_path = define_input_source(args.input_pan)

    # Check if gene_data file is present if Panaroo input is given an gffs should be annotated
    if args.annotate and source_program is not 'Rorary':
        gene_data_path = check_gene_data(args.input_pan)

    if not args.quiet:
        print(f"Pan genome determined to come from {source_program}")
        print("All files found, let's move on!\n")
        print("--------------------------------------------------------------\n")

    # TODO - Make the program work with less than all files in the pangenome. Just make sure that all gff files supplied can be found in the pan genome. This will make is possible to look at hotspots and segments in different lineages
    check_gff_in_pan(args.input_gffs, input_pres_abs_file_path)

    # Construct output folder
    try:
        mkdir(args.output_path)
        if not args.quiet:
            print("Output folder constructed")
    except FileExistsError:
        if not args.quiet:
            print("Output folder exists")

    # Construct temporary folder:
    # TODO - check that the temporary folder does not exist and that the user does not have a folder with same name already. (Maybe use a time stamp for the start to make it unique.)
    temp_folder_path = join(args.output_path, 'genome_corer_tmp')
    mkdir(temp_folder_path)

    ## Read in gene presence absence file
    time_start = time.time()
    # TODO - Add the user specified thresholds for core and low frequency genes.
    core_dict, low_freq_dict, acc_gene_dict, attribute_dict = read_gene_presence_absence(input_pres_abs_file_path,
                                                                                         1, 0.05, source_program,
                                                                                         args.input_gffs, temp_folder_path)
    if not args.quiet:
        time_calculator(time_start, time.time(), "reading in gene presence/absence file")

    # If source program is Panaroo, see if alignments folder is available, and core genes are contained in it.
    # if source_program == 'Panaroo':
    #     alignment_folder = check_gene_alignments(args.input_pan, core_dict)

    # Add in the refound genes into the gff files and print the corrected GFF files.
    if source_program == "Panaroo" and args.annotate:
        time_start = time.time()
        print(f"\n----------Adding in refound annotations for gff files---------")

        corrected_folder = correct_gffs(args.input_gffs, gene_data_path, args.output_path, attribute_dict, temp_folder_path)

        args.input_gffs = [join(corrected_folder, file) for file in listdir(corrected_folder) if '.gff' in file]
        if not args.quiet:
            time_calculator(time_start, time.time(), "add refound annotations to gff files")


    # Loop over all gffs and extract info from each of them.
    time_start = time.time()
    # Initialise dictionaries to contain results from all gff files
    core_neighbour_pairs = {}
    core_neighbour_distance = {}
    core_neighbour_accessory_count = {}
    core_neighbour_low_freq = {}
    master_info_total = {}
    non_core_contig_info = {}
    merged_start_gene_clusters = []
    merged_second_gene_clusters = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=15) as executor:
        print(f"\n------Start core region identification of given gff files-----")
        print(f'{len(args.input_gffs)} GFF files to process')

        results = [executor.submit(segment_genome_content, gff, core_dict, low_freq_dict, acc_gene_dict, i, comp_genomes)
                   for i, gff in enumerate(args.input_gffs)]

        for output in concurrent.futures.as_completed(results):
            # Split the outputs
            core_pairs, distance, acc_count, \
            low_freq, master_info_return, \
            core_less_contigs_return, start_gene_cluster = output.result()

            # Merge results into single/master dictionaries
            core_neighbour_pairs = merge_dicts_counts(core_neighbour_pairs, core_pairs)
            core_neighbour_distance = merge_dicts_lists(core_neighbour_distance, distance)
            core_neighbour_accessory_count = merge_dicts_lists(core_neighbour_accessory_count, acc_count)
            core_neighbour_low_freq = merge_dicts_lists(core_neighbour_low_freq, low_freq)
            master_info_total.update(master_info_return)
            non_core_contig_info.update(core_less_contigs_return)
            merged_start_gene_clusters, merged_second_gene_clusters = merge_first_genes(start_gene_cluster,
                                                                                        merged_start_gene_clusters,
                                                                                        merged_second_gene_clusters,
                                                                                        core_pairs[0])
    time_calculator(time_start, time.time(), "searching gff files for core genes")


    # TODO - Try to identify variable regions
    #   * Check if there are any sequence breaks in the dataset
    #   * if sequence breaks are present, examine which core genes could substitute sequence breaks
    #   If sequence break has single substitute and variance is observed in length or accessory content then flag as variable
    #   If sequence break has multipl substitutes then flag as variable anyway, as this eludes to rearrangements


    ### FUNCTION ###
    # Determine the most common core gene synteny.
    # time_start = time.time()
    # Find the core gene synteny and possible core genes with alternative neighbours
    # consensus_core_genome, \
    #     possible_rearrangement_genes, \
    #     core_path_coverage = determine_core_gene_consesus(core_neighbour_pairs,
    #                                                       merged_start_gene_clusters,
    #                                                       merged_second_gene_clusters, args.output_path)

    ### Determine segments of core genome ###
    print(f"\n--------------Identifying segments in pan genome--------------")
    time_start = time.time()
    # Count number of unique accessory genes inserted into a core-core region across the genomes
    acc_region_count = {key: len(set(core_neighbour_low_freq[key])) for key in core_neighbour_low_freq}
    # Count number of unique low frequency genes inserted into a core-core region across the genomes
    low_frew_region_count = {key: len(set(core_neighbour_accessory_count[key])) for key in core_neighbour_accessory_count}

    # Combine the accessory and low frequency counts:
    combined_acc_gene_count = {key: low_frew_region_count[key] + acc_region_count[key] for key in low_frew_region_count}

    double_edge_segements, no_acc_segments = determine_genome_segments(core_neighbour_pairs, combined_acc_gene_count, len(args.input_gffs))

    time_calculator(time_start, time.time(), "identifying segments in pan genome")
    # Assign core-gene synteny types:
    # Identify alternative connections and their occurrence
    # alt_core_pairs, alt_core_pair_count, \
    #     core_genome_types, alt_core_comp_types = identify_rearrangements(consensus_core_genome,
    #                                                                      possible_rearrangement_genes,
    #                                                                      master_info_total,
    #                                                                      args.input_gffs)

    # rearrangement_predictions = characterise_rearrangements(alt_core_pairs, consensus_core_genome)

    # time_calculator(time_start, time.time(), "determining best core gene synteny")
    ### DO CALCULATIONS ###
    # TODO mean number length between core genes
    # for neighbours in core_neighbour_distance:
    #     print(mean(core_neighbour_distance[neighbours]))
    #     print(std(core_neighbour_distance[neighbours]))

    #######################
    # TODO - make function that relates every core gene to a given reference genomes' locus_tags, if given such a reference. - On request from Andrew


    # Make a Summary table like the one produced in R.
    master_summary_info = calculate_n_create_summaries(master_info_total)



    ### WRITE OUTPUTS ###
    print(f"\n-----------------------Printing outputs-----------------------")
    # Write master information to output file
    time_start = time.time()
    master_info_writer(master_info_total, args.output_path, args.output_prefix, args.quiet)
    summary_info_writer(master_summary_info, args.output_path, args.output_prefix, args.quiet)
    # if return_of_segments is not None
    # TODO - Contruct output for segments - parent column.
    segment_writer(double_edge_segements, args.output_path, args.output_prefix, args.quiet)
    no_acc_segment_writer(no_acc_segments, args.output_path, args.output_prefix, args.quiet)
    # print(non_core_contig_info) TODO - Print core less contigs.
    # TODO print a list of accessory genes that have not been related to any region?
    # TODO - Possibly output core gene graph. with segment annotations?

    time_calculator(time_start, time.time(), "writing output files")

    # Finish up running
    time_calculator(total_time_start, time.time(), "running the entire program")

    # Remove temporary database holding gff databases
    # TODO - Implement a nice crash function where the temporary folder is removed not to cause unessecary frustration for the user when trying to rerun the program.
    print(isdir(temp_folder_path))
    if isdir(temp_folder_path):
        rmdir(temp_folder_path)

if __name__ == "__main__":
    main()
