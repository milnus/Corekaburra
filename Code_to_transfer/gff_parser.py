import numpy as np
import os


def parse_gff(input_file):
    """ Read a gff file and return it as a generator object that return all line containing CDS """
    with open(input_file, 'r') as gff_file:
        for line in gff_file:
            if "##FASTA" in line:
                break
            if "#" not in line and 'CDS' in line:
                # Strip line for newline and split columns into list
                line = line.strip()
                line = line.split("\t")

                # See if refound gene or Prokka annotated and isolate ID in gene_presence_absence.csv accordingly
                if "old_locus_tag=" in line[8]:
                    gene_id = line[8][line[8].find('old_locus_tag'):]
                else:
                    gene_id = line[8][line[8].find('ID'):line[8].find(';')]

                # Remove equal sign from id and add as identifyer for the returned gff line
                gene_id = gene_id[gene_id.find('=') + 1:]
                line[8] = gene_id
                yield line


def get_genome_size_from_gff(input_file): #TODO - it should be possible to remove this and just use the directory given by get_contig_lengths function
    """ Get the genome size from a GFF3 file by counting characters in fasta appendix"""
    genome_size = 0
    fasta_reached = False
    with open(input_file, 'r', ) as gff_file:
        for line in gff_file:
            if fasta_reached and '>' not in line:
                genome_size += len(line.rstrip())
            if "##FASTA" in line:
                fasta_reached = True

    return genome_size


def get_contig_lengths(input_file):
    fasta_reached = False
    contig_size = 0
    contig_size_dir = {}

    # Open the given gff file, find the fasta section of the file and count the length of each contig.
    with open(input_file, 'r', ) as gff_file:
        for line in gff_file:
            if fasta_reached and '>' not in line:
                contig_size += len(line.rstrip())
            if fasta_reached and '>' in line:
                if contig_size > 0:
                    # Record the previous contig
                    if contig_name not in contig_size_dir:
                        contig_size_dir[contig_name] = contig_size
                    else:
                        raise ValueError(f"contig name: {contig_name}, in file {input_file} is duplicated! Please fix this")

                    # Set the contig name to the next contig
                    contig_name = line.strip().split(' ')[0].replace('>', '')
                    contig_size = 0
                else:
                    contig_name = line.strip().split(' ')[0].replace('>', '')

            if "##FASTA" in line:
                fasta_reached = True

        # Record last contig
        contig_size_dir[contig_name] = contig_size

    return contig_size_dir


def record_core_core_region(core_genes, gff_name, gff_line, contig_end, previous_core_gene_id,
                            previous_core_gene_end_coor, acc_genes_in_region, low_freq_genes_in_region,
                            core_gene_pair_distance, accessory_gene_content,
                            low_freq_gene_content, core_gene_pairs, num_acc_genes_in_region, master_info):

    """ Function to record information about a core gene pair or a core gene and a sequence break,
    along with accessory information between the two features"""

    # Examine if a fragmented core genome is given, if then skip the recording, reset the accessory counters and return the new end coordinate of the second fragment
    if gff_line is not None and previous_core_gene_id is not "Sequence_break":
        if core_genes[gff_name][previous_core_gene_id] == core_genes[gff_name][gff_line[8]]:
            previous_core_gene_id = gff_line[8]
            previous_core_gene_end_coor = int(gff_line[4])
            acc_genes_in_region = []
            low_freq_genes_in_region = []
            return (previous_core_gene_id, previous_core_gene_end_coor, acc_genes_in_region, low_freq_genes_in_region,
                    core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, core_gene_pairs,
                    num_acc_genes_in_region, master_info)


    # Set core cluster names
    # If no line from gff is given there is a sequence break,
    # if it is given then set current cluster and try to find previous if not found it is a sequence break
    if gff_line is not None:
        current_core_gene_cluster = core_genes[gff_name][gff_line[8]]
        try:
            previous_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
        except KeyError:
            previous_core_gene_cluster = previous_core_gene_id

        core_gene_neighbours = sorted([previous_core_gene_cluster, current_core_gene_cluster])

    else:
        current_core_gene_cluster = "Sequence_break"
        previous_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
        core_gene_neighbours = [previous_core_gene_cluster, current_core_gene_cluster]

    # Merge core neighbours and add to pairs
    core_gene_neighbours_str = f'{core_gene_neighbours[0]}--{core_gene_neighbours[1]}'
    core_gene_pairs.append(core_gene_neighbours_str)

    # Set core neighbour distance
    # Check if measuring between two genes on same contig and not measuring from a sequence break
    if gff_line is not None and previous_core_gene_cluster is not "Sequence_break":
        core_core_distance = int(gff_line[3]) - previous_core_gene_end_coor - 1
    else:
        # Check if measuring between sequence break and first core gene on contig, if then set start coordinate to zero
        if previous_core_gene_id == "Sequence_break":
            core_core_distance = int(gff_line[3]) - 1
            # If gene starts at contig end, make distance zero
            if core_core_distance < 0:
                core_core_distance = 0
        elif current_core_gene_cluster == 'Sequence_break':
            core_core_distance = abs(contig_end - previous_core_gene_end_coor - 1)
        else:
            NotImplementedError(
                'An error occured when measuring the distance between core gene and contig end. Something went wrong!')

    # Add core neighbour distance
    core_gene_pair_distance[core_gene_neighbours_str] = core_core_distance
    # Add number of accessory genes in region
    #num_acc_genes_in_region[core_gene_neighbours_str] = len(acc_genes_in_region) + len(low_freq_genes_in_region)

    # Add counts and annotation for accessory and low frequency genes
    accessory_gene_content[core_gene_neighbours_str] = acc_genes_in_region.copy() # Copies has been added here
    low_freq_gene_content[core_gene_neighbours_str] = low_freq_genes_in_region.copy()

    # Add info to master dict
    master_info[f'{core_gene_neighbours_str}--{gff_name}'] = [gff_name,
                                                              core_gene_neighbours[0],
                                                              core_gene_neighbours[1],
                                                              core_core_distance,
                                                              len(acc_genes_in_region) + len(low_freq_genes_in_region),
                                                              acc_genes_in_region.copy(),
                                                              low_freq_genes_in_region.copy()]

    # Update previous core gene id and end of core gene
    if gff_line is not None:
        # if previous_core_gene_id is not "Sequence_break":
        previous_core_gene_id = gff_line[8]
        previous_core_gene_end_coor = int(gff_line[4])
    else:
        previous_core_gene_id = "Sequence_break"
        # previous_core_gene_end_coor = 0

    # Reset values for accessory genes next core-core region
    acc_genes_in_region = []
    low_freq_genes_in_region = []

    return (previous_core_gene_id, previous_core_gene_end_coor, acc_genes_in_region, low_freq_genes_in_region,
            core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, core_gene_pairs,
            num_acc_genes_in_region, master_info)


def connect_first_n_last_gene_on_contig(core_genes, gff_name, previous_core_gene_id, previous_core_gene_end_coor,
                                        first_core_gene_gff_line, acc_genes_in_region, first_core_accessory_content,
                                        low_freq_genes_in_region, first_core_low_freq_genes, contig_size,
                                        core_gene_pairs, core_gene_pair_distance,
                                        accessory_gene_content, low_freq_gene_content, master_info):

    # TODO - Take first and last gene on chromosome and record them as connected.

    if previous_core_gene_id == 'Complete_genome_end_fail':
        print('Complete_genome_end_fail given!')
        previous_core_gene_id = first_core_gene_gff_line[8]
        # raise ValueError(f'Complete_genome_end_fail was given as input!\n'
        #                  f'gff: {gff_name}\n'
        #                  f'last_core_gene_cluster: {core_genes[gff_name][first_core_gene_gff_line[8]]}')

    last_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
    first_core_gene_cluster = core_genes[gff_name][first_core_gene_gff_line[8]]

    # TODO - Handle core gene that is on contig alone?
    if first_core_gene_cluster == last_core_gene_cluster:
        print('Same gene')
    #     raise ValueError(f'First and last core gene is the same!\n'
    #                      f'gff: {gff_name}\n'
    #                      f'first_core_gene_cluster: {first_core_gene_cluster}\n'
    #                      f'last_core_gene_cluster: {last_core_gene_cluster}')

    # Add core neighbours
    core_gene_neighbours = sorted([last_core_gene_cluster, first_core_gene_cluster])
    core_gene_neighbours_str = f'{core_gene_neighbours[0]}--{core_gene_neighbours[1]}'
    core_gene_pairs.append(core_gene_neighbours_str)

    # Add core neighbour distance
    core_core_distance = contig_size - previous_core_gene_end_coor \
                         + int(first_core_gene_gff_line[3])  # TODO - make sure this -1 is correct

    core_gene_pair_distance[core_gene_neighbours_str] = core_core_distance

    # Add accessory information from between last and first core gene
    last_first_accessory_content = acc_genes_in_region.copy() + first_core_accessory_content.copy()
    last_first_low_freq_count = low_freq_genes_in_region.copy() + first_core_low_freq_genes.copy()

    accessory_gene_content[core_gene_neighbours_str] = last_first_accessory_content.copy()  # Copies has been added here
    low_freq_gene_content[core_gene_neighbours_str] = last_first_low_freq_count.copy()

    # Add to master info dict
    master_info[f'{core_gene_neighbours_str}--{gff_name}'] = [gff_name,
                                                              core_gene_neighbours[0],
                                                              core_gene_neighbours[1],
                                                              core_core_distance,
                                                              len(last_first_accessory_content) +
                                                              len(last_first_low_freq_count),
                                                              last_first_accessory_content,
                                                              last_first_low_freq_count]

    previous_core_gene_id = 'Complete_genome_end_fail'
    previous_core_gene_end_coor = int(first_core_gene_gff_line[4])
    acc_genes_in_region = []
    low_freq_genes_in_region = []

    return (previous_core_gene_id, previous_core_gene_end_coor, acc_genes_in_region, low_freq_genes_in_region,
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info)


def segment_gff_content(gff_generator, core_genes, low_freq_genes, gff_path, acc_genes, complete_genomes):
    """ Function that takes a gff generator, core and low frequency genes and identify each core-core gene region
    counts the number of accessory genes in the region, records the number of low frequency genes in the region,
    and records the distance from one core gene to the next"""

    # Initialize data structures to be returned
    # List of sorted core gene neighbours separated by '--'. (Sorted as the orientation of contigs is sometimes random)
    core_gene_pairs = []
    # Dict of distances between core gene neighbours. Key is the core_gene_pairs, and values is the distance in basepairs from end to start for neighbouring genes.
    core_gene_pair_distance = {}
    # Dict of number of accessory genes between core gene neighbours, Keys are core_gene_pairs and values are the number of accessory genes in the region
    num_acc_genes_in_region = {}
    # Dict of low fequency genes between core gene pairs. Key is value from core_gene_pairs, and value is the names of low frequency genes found in region
    low_freq_gene_content = {}
    # Dict of lists with keys being the core_gene_pairs value with the value being the number of accessory genes found between the core genes
    accessory_gene_content = {}
    # Dict containing master information - used to write comprehensive output files
    master_info = {}
    # Dict to store info on accessory genes from contigs where no core gene is present.
    coreless_contigs = {}
    # initiate variable that holds the first gene for file is genome is complete. - # How do you determine from a gff file if a genome is complete?
    start_gene_cluster = False
    # Examine if a complete genomes has been given
    if complete_genomes is None:
        complete_genome = False
    else:
        if os.path.basename(gff_path).replace('.gff', '').replace('_corrected', '') in complete_genomes: # TODO - Use remove instead of replace?
            complete_genome = True
        else:
            complete_genome = False

    # Get size of contigs for given gff
    contig_sizes = get_contig_lengths(gff_path)

    # Split input path to gff to get genome name
    gff_name = gff_path.split('/')[-1]
    gff_name = gff_name.split('.')[0]
    if 'corrected' in gff_name:
        gff_name = gff_name.split('_corrected')[0]

    # Set that first core gene has not been found
    first_core_gene = True
    # Set Initialise variable to be used globally
    previous_core_gene_id = ""
    # Set variable to determine of genome is complete or single contig, or if multi contigs are present
    #single_contig = True
    # Set variable to get first contig
    first_contig = True
    # Initialise the accessory gene counter and low frequency gene list
    low_freq_genes_in_region = []
    acc_genes_in_region = []

    # Go through each line of GFF file
    for line in gff_generator:
        # Set first contig fofund in file
        if first_contig:
            previous_contig = line[0]
            first_contig = False

        # Check if contig has changed - if then finish contig, if not examine next gene on contig gene
        if line[0] == previous_contig:

            # Check if gene is core - if then test if fist or not, else assume to be accessory gene
            if line[8] in core_genes[gff_name]:

                # Check if core gene is the first observed in file - if then set information else record information
                if first_core_gene:
                    # Set information on first core gene to be used when finishing search
                    first_core_gene_gff_line = line
                    first_core_accessory_content = acc_genes_in_region.copy()
                    first_core_low_freq_genes = low_freq_genes_in_region.copy()

                    # Set information to be used with next core gene neighbour
                    previous_core_gene_end_coor = int(line[4])
                    previous_core_gene_id = line[8]

                    # Set that first core gene has been observed
                    first_core_gene = False

                    # Reset accessory and low frequency gene counters
                    acc_genes_in_region = []
                    low_freq_genes_in_region = []

                else:
                    # Record core gene pair information
                    (previous_core_gene_id,
                     previous_core_gene_end_coor,
                     acc_genes_in_region,
                     low_freq_genes_in_region,
                     core_gene_pair_distance,
                     accessory_gene_content,
                     low_freq_gene_content,
                     core_gene_pairs,
                     num_acc_genes_in_region,
                     master_info) = record_core_core_region(core_genes, gff_name, line, None, previous_core_gene_id,
                                                            previous_core_gene_end_coor, acc_genes_in_region,
                                                            low_freq_genes_in_region, core_gene_pair_distance,
                                                            accessory_gene_content, low_freq_gene_content,
                                                            core_gene_pairs, num_acc_genes_in_region, master_info)

            else:
                # Check if accessory is low frequency - else just regular accessory
                if line[8] in low_freq_genes[gff_name]:
                    low_freq_genes_in_region.append(low_freq_genes[gff_name][line[8]])
                else:
                    try:
                        acc_genes_in_region.append(acc_genes[gff_name][line[8]])
                    except KeyError:
                        gene_key = [key for key in acc_genes[gff_name].keys() if line[8] in key]
                        if len(gene_key) > 1:
                            acc_genes_in_region.append(acc_genes[gff_name][gene_key][0])

        else:
            # Note that gff has multiple contigs
            #single_contig = False
            # Check if there is a core gene on traversed contig or if a core gene is present on the first contig -
            # if then record it, if not record the accessory and low frequency genes found on contig and reset.
            if previous_core_gene_id is not "Sequence_break" and previous_core_gene_id is not "":
                if complete_genome: # TODO - Write a unit check for this good danm thing! - May be wrap in function to be used further down too, when last line in file has been read? ***
                    (previous_core_gene_id,
                     previous_core_gene_end_coor,
                     acc_genes_in_region,
                     low_freq_genes_in_region,
                     core_gene_pairs,
                     core_gene_pair_distance,
                     accessory_gene_content,
                     low_freq_gene_content,
                     master_info) = connect_first_n_last_gene_on_contig(core_genes, gff_name, previous_core_gene_id,
                                                                        previous_core_gene_end_coor,
                                                                        first_core_gene_gff_line, acc_genes_in_region,
                                                                        first_core_accessory_content,
                                                                        low_freq_genes_in_region,
                                                                        first_core_low_freq_genes, contig_sizes[line[0]],
                                                                        core_gene_pairs, core_gene_pair_distance,
                                                                        accessory_gene_content, low_freq_gene_content,
                                                                        master_info)

                    first_core_gene = True

                    # last_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
                    # first_core_gene_cluster = core_genes[gff_name][first_core_gene_gff_line[8]]
                    #
                    # # Add core neighbours
                    # core_gene_neighbours = sorted([last_core_gene_cluster, first_core_gene_cluster])
                    # core_gene_neighbours_str = f'{core_gene_neighbours[0]}--{core_gene_neighbours[1]}'
                    # core_gene_pairs.append(core_gene_neighbours_str)
                    #
                    # # Add core neighbour distance
                    # core_core_distance = contig_sizes[first_core_gene_gff_line[0]] - previous_core_gene_end_coor \
                    #                      + int(first_core_gene_gff_line[3])  # TODO - make sure this -1 is correct
                    #
                    # core_gene_pair_distance[core_gene_neighbours_str] = core_core_distance
                    #
                    # # Add accessory information from between last and first core gene
                    # last_first_accessory_content = acc_genes_in_region + first_core_accessory_content
                    # accessory_gene_content[core_gene_neighbours_str] = len(last_first_accessory_content)
                    # last_first_low_freq_count = low_freq_genes_in_region + first_core_low_freq_genes
                    # low_freq_gene_content[core_gene_neighbours_str] = last_first_low_freq_count
                    #
                    # # Add to master info dict
                    # master_info[f'{core_gene_neighbours_str}--{gff_name}'] = [gff_name,
                    #                                                           core_gene_neighbours[0],
                    #                                                           core_gene_neighbours[1],
                    #                                                           core_core_distance,
                    #                                                           len(last_first_accessory_content) +
                    #                                                           len(last_first_low_freq_count),
                    #                                                           last_first_accessory_content,
                    #                                                           last_first_low_freq_count]

                # TODO - check if genome is given as complete, if then record the fist and last gene as being neighbours and set first core gene to True to record a new first gene.
                #        Else add the core gene as being next to a seqeunce break
                else:
                    # Record the core gene neighbouring a sequence break
                    (previous_core_gene_id,
                     previous_core_gene_end_coor,
                     acc_genes_in_region,
                     low_freq_genes_in_region,
                     core_gene_pair_distance,
                     accessory_gene_content,
                     low_freq_gene_content,
                     core_gene_pairs,
                     num_acc_genes_in_region,
                     master_info) = record_core_core_region(core_genes, gff_name, None, contig_sizes[previous_contig], previous_core_gene_id,
                                                            previous_core_gene_end_coor, acc_genes_in_region,
                                                            low_freq_genes_in_region, core_gene_pair_distance,
                                                            accessory_gene_content, low_freq_gene_content,
                                                            core_gene_pairs, num_acc_genes_in_region,
                                                            master_info)

                    # Check if first gene on contig is a core gene, if the record it.
                    if line[8] in core_genes[gff_name]:
                        previous_core_gene_id = "Sequence_break"

                        # Get the starting position of the first core gene on contig to record the gene.
                        # Make it negative to fit the calculation of the distance between genes. Add one to further adjust
                        cur_core_gene_start = -int(line[3])

                        (previous_core_gene_id,
                         previous_core_gene_end_coor,
                         acc_genes_in_region,
                         low_freq_genes_in_region,
                         core_gene_pair_distance,
                         accessory_gene_content,
                         low_freq_gene_content,
                         core_gene_pairs,
                         num_acc_genes_in_region,
                         master_info) = record_core_core_region(core_genes, gff_name, line, 0, previous_core_gene_id,
                                                                cur_core_gene_start, acc_genes_in_region,
                                                                low_freq_genes_in_region, core_gene_pair_distance,
                                                                accessory_gene_content, low_freq_gene_content,
                                                                core_gene_pairs, num_acc_genes_in_region,
                                                                master_info)

                    # Add as accessory
                    else:
                        # Check if accessory is low frequency - else just regular accessory
                        if line[8] in low_freq_genes[gff_name]:
                            low_freq_genes_in_region.append(low_freq_genes[gff_name][line[8]])
                        else:
                            try:
                                acc_genes_in_region.append(acc_genes[gff_name][line[8]])
                            except KeyError:
                                gene_key = [key for key in acc_genes[gff_name].keys() if line[8] in key]
                                if len(gene_key) > 1:
                                    acc_genes_in_region.append(acc_genes[gff_name][gene_key][0])

                # TDOD - Add an else that count the accessory/low frequency gene.
                # Set new contig
                previous_contig = line[0]

            else:
                # Record info on accessory genes on core-less contig, if any accessory genes are present
                if len(acc_genes_in_region) + len(low_freq_genes_in_region) > 0:
                    coreless_contigs[f'{gff_name}--{previous_contig}'] = [acc_genes_in_region, low_freq_genes_in_region]
                acc_genes_in_region = []
                low_freq_genes_in_region = []

                # Set new contig
                previous_contig = line[0]

    # Check if genome is complete or a single contig. If then add information for last and first core gene, if not
    # then add the first and last core gene as being neighbours to sequence breaks.
    if complete_genome:
        if previous_core_gene_id is not 'Complete_genome_end_fail':
            (previous_core_gene_id,
             previous_core_gene_end_coor,
             acc_genes_in_region,
             low_freq_genes_in_region,
             core_gene_pairs,
             core_gene_pair_distance,
             accessory_gene_content,
             low_freq_gene_content,
             master_info) = connect_first_n_last_gene_on_contig(core_genes, gff_name, previous_core_gene_id,
                                                                previous_core_gene_end_coor,
                                                                first_core_gene_gff_line, acc_genes_in_region,
                                                                first_core_accessory_content,
                                                                low_freq_genes_in_region,
                                                                first_core_low_freq_genes, contig_sizes[line[0]],
                                                                core_gene_pairs, core_gene_pair_distance,
                                                                accessory_gene_content, low_freq_gene_content,
                                                                master_info)
        else:
            # Record info on accessory genes on core-less contig, if any accessory genes are present
            if len(acc_genes_in_region) + len(low_freq_genes_in_region) > 0:
                coreless_contigs[f'{gff_name}--{previous_contig}'] = [acc_genes_in_region, low_freq_genes_in_region]


        # last_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
        # first_core_gene_cluster = core_genes[gff_name][first_core_gene_gff_line[8]]
        #
        # # Add core neighbours
        # core_gene_neighbours = sorted([last_core_gene_cluster, first_core_gene_cluster])
        # core_gene_neighbours_str = f'{core_gene_neighbours[0]}--{core_gene_neighbours[1]}'
        # core_gene_pairs.append(core_gene_neighbours_str)
        #
        # # Add core neighbour distance
        # core_core_distance = contig_sizes[first_core_gene_gff_line[0]] - previous_core_gene_end_coor \
        #                     + int(first_core_gene_gff_line[3]) # TODO - make sure this -1 is correct
        #
        # core_gene_pair_distance[core_gene_neighbours_str] = core_core_distance
        #
        # # Add accessory information from between last and first core gene
        # last_first_accessory_content = acc_genes_in_region + first_core_accessory_content
        # accessory_gene_content[core_gene_neighbours_str] = len(last_first_accessory_content)
        # last_first_low_freq_count = low_freq_genes_in_region + first_core_low_freq_genes
        # low_freq_gene_content[core_gene_neighbours_str] = last_first_low_freq_count
        #
        # # Add to master info dict
        # master_info[f'{core_gene_neighbours_str}--{gff_name}'] = [gff_name,
        #                                                      core_gene_neighbours[0],
        #                                                      core_gene_neighbours[1],
        #                                                      core_core_distance,
        #                                                       len(last_first_accessory_content) +
        #                                                       len(last_first_low_freq_count),
        #                                                      last_first_accessory_content,
        #                                                      last_first_low_freq_count]

        # Add first gene in file to list, to aid as start point for consensus core gene graph
        # start_gene_cluster = first_core_gene_cluster

    else:
        # Add first core gene as being neighbour to a sequence break
        (_,
         _,
         _,
         _,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content, core_gene_pairs,
         num_acc_genes_in_region, master_info) = record_core_core_region(core_genes, gff_name, first_core_gene_gff_line, 0,
                                                                         'Sequence_break',
                                                                         9999999999999999,
                                                                         first_core_accessory_content,
                                                                         first_core_low_freq_genes,
                                                                         core_gene_pair_distance,
                                                                         accessory_gene_content,
                                                                         low_freq_gene_content,
                                                                         core_gene_pairs,
                                                                         num_acc_genes_in_region,
                                                                         master_info)
        # Add last core gene as being neighbour to a sequence break
        # if the last core gene has not been recorded already
        if previous_core_gene_id is not "Sequence_break":
            (previous_core_gene_id,
             previous_core_gene_end_coor,
             acc_genes_in_region,
             low_freq_genes_in_region,
             core_gene_pair_distance,
             accessory_gene_content,
             low_freq_gene_content, core_gene_pairs,
             num_acc_genes_in_region, master_info) = record_core_core_region(core_genes, gff_name, None, contig_sizes[previous_contig],
                                                                             previous_core_gene_id,
                                                                             previous_core_gene_end_coor,
                                                                             acc_genes_in_region,
                                                                             low_freq_genes_in_region,
                                                                             core_gene_pair_distance,
                                                                             accessory_gene_content,
                                                                             low_freq_gene_content,
                                                                             core_gene_pairs,
                                                                             num_acc_genes_in_region,
                                                                             master_info)

    return core_gene_pairs, core_gene_pair_distance, accessory_gene_content, \
           low_freq_gene_content, master_info, coreless_contigs, start_gene_cluster


def segment_genome_content(input_file, core_genes, low_freq_genes, acc_gene_dict, i, complete_genomes):
    """ Single function segmenting the gff into core gene regions to be used for simple multi processing"""
    if (i+1) % 25 == 0 or i == 0:
        print(f"Determining core-core synteny for GFF file #{i+1}")
    gff_generator = parse_gff(input_file)
    return_data = segment_gff_content(gff_generator=gff_generator,
                                      gff_path=input_file,
                                      core_genes=core_genes,
                                      low_freq_genes=low_freq_genes,
                                      acc_genes=acc_gene_dict,
                                      complete_genomes=complete_genomes)

    return return_data


# TODO - Add in so that the length between a given core gene and a sequnece break is conserved and given in master output.