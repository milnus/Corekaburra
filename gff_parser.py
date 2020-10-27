from parse_gene_presence_absence import read_gene_presence_absence


def parse_gff(input_file):
    """ Read a gff file and return it as a generator object that return all line containing CDS """
    with open(input_file, 'r', ) as gff_file:
        for line in gff_file:
            if "##FASTA" in line:
                break
            if "#" not in line and 'CDS' in line:
                line = line.strip()
                line = line.split("\t")
                gene_id = line[8][line[8].find('ID'):line[8].find(';')]
                gene_id = gene_id[gene_id.find('=')+1:]
                line[8] = gene_id
                yield line


def get_genome_size_from_gff(input_file):
    """ Get the genome size from a GFF3 file by counting characters in fasta appendix"""
    genome_size = 0
    fasta_reached = False
    with open(input_file, 'r', ) as gff_file:
        for line in gff_file:
            if "##FASTA" in line:
                fasta_reached = True
            if fasta_reached:
                genome_size += len(line.strip())

    return genome_size


def record_core_core_region(core_genes, gff_name, gff_line, previous_core_gene_id,
                            previous_core_gene_end_coor, accessory_gene_count, low_freq_genes_in_region,
                            core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, core_gene_pairs,
                            master_info):

    """ Function to record information about a core gene pair or a core gene and a sequence break,
    along with accessory information between the two features"""

    # Set core cluster names
    if gff_line is not None:
        current_core_gene_cluster = core_genes[gff_name][gff_line[8]]
        try:
            previous_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
        except KeyError:
            previous_core_gene_cluster = previous_core_gene_id

    else:
        current_core_gene_cluster = "Sequence_break"
        previous_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]

    # Add core neighbours
    core_gene_neighbours = sorted([previous_core_gene_cluster, current_core_gene_cluster])
    core_gene_neighbours_str = f'{core_gene_neighbours[0]}-{core_gene_neighbours[1]}'
    core_gene_pairs.append(core_gene_neighbours_str)

    # Add core neighbour distance
    # TODO - Check possible ovrelap and standedness - Maybe
    if gff_line is not None:
        core_core_distance = [int(gff_line[3]) - previous_core_gene_end_coor]
        core_gene_pair_distance[core_gene_neighbours_str] = core_core_distance
    else:
        core_core_distance = ['']



    # Add counts and annotation for accessory and low frequency genes
    accessory_gene_content[core_gene_neighbours_str] = accessory_gene_count
    low_freq_gene_content[core_gene_neighbours_str] = low_freq_genes_in_region

    # Add info for master dict
    master_info[f'{core_gene_pairs[-1]}-{gff_name}'] = [gff_name,
                                                            core_gene_neighbours[0],
                                                            core_gene_neighbours[1],
                                                            core_core_distance,
                                                            accessory_gene_count,
                                                            low_freq_genes_in_region]

    # Update previous core gene id and end of core gene
    if gff_line is not None:
        # if previous_core_gene_id is not "Sequence_break":
        previous_core_gene_id = gff_line[8]
        previous_core_gene_end_coor = int(gff_line[4])
    else:
        previous_core_gene_id = "Sequence_break"

    # Reset values for accessory genes next core-core region
    accessory_gene_count = 0
    low_freq_genes_in_region = []

    return (previous_core_gene_id, previous_core_gene_end_coor, accessory_gene_count, low_freq_genes_in_region,
            core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, core_gene_pairs, master_info)


def segment_gff_content(gff_generator, core_genes, low_freq_genes, gff_path):
    """ Function that takes a gff generator, core and low frequency genes and identify each core-core gene region
    counts the number of accessory genes in the region, records the number of low frequency genes in the region,
    and records the distance from one core gene to the next"""

    # TODO Add in a dict that has an entery for every core-core region in all genomes (key = core_1-core_2-gff_file), where the value is a list with gff_name, core_1, core_2, distance, accessory gene count, low_frequency genes in list
    # TODO - Integrate refound genes from somewhere!
    # TODO - Make it so that all accessory genes are recorded by their name in the core-core region. This would allow to
    # To observe movement in the pan-genome of individual accessory genes.

    # Initialize data structures to be returned
    # List of sorted core gene neighbours separated by '-'. (Sorted as the orientation of contigs is sometimes random)
    core_gene_pairs = []
    # Dict of distances between core gene neighbours. Key is the core_gene_pairs value, and values is the distance in basepairs from end to start.
    core_gene_pair_distance = {}
    # Dict of low fequency genes between core gene pairs. Key is value from core_gene_pairs, and value is the names of low frequency genes found in region
    low_freq_gene_content = {}
    # Dict of lists with keys being the core_gene_pairs value with the value being the number of accessory genes found between the core genes
    accessory_gene_content = {}
    # Dict containing master information - used to write comprehensive output files
    master_info = {}

    # Split input path to gff to get genome name
    gff_name = gff_path.split('/')[-1]
    gff_name = gff_name.split('.')[0]

    # Set that first core gene has not been found
    first_core_gene = True
    # Set Initialise variable to be used globally
    previous_core_gene_id = ""
    # previous_core_gene_end_coor = 0
    # first_core_gene_start_coor = 0
    # first_core_gene_id = "dummy_core_id"
    # first_core_accessory_count = 0
    # first_core_low_freq_genes = []
    # Set variable to determine of genome is complete or single contig, or if multi contigs are present
    single_contig = True
    # Set variable to get first contig
    first_contig = True
    # Initialise the accessory gene counter and low frequency gene list
    accessory_gene_count = 0
    low_freq_genes_in_region = []

    # Go through each line of GFF
    for line in gff_generator:
        # Set first contig
        if first_contig:
            previous_contig = line[0]
            first_contig = False

        # Check if contig has changed - if then finish contig, if not examine gene
        if line[0] == previous_contig:

            # Check if gene is core
            if line[8] in core_genes[gff_name]:
                # Check if core gene is the first
                if first_core_gene:
                    # Set information on first core gene to be used if complete genome
                    first_core_gene_start_coor = int(line[3])
                    first_core_gene_id = line[8]
                    first_core_accessory_count = accessory_gene_count
                    first_core_low_freq_genes = low_freq_genes_in_region
                    first_core_gene = False

                    # Set information to be used with next core gene neighbour
                    previous_core_gene_end_coor = int(line[4])
                    previous_core_gene_id = line[8]

                else:
                    # Record core gene pair information
                    (previous_core_gene_id,
                     previous_core_gene_end_coor,
                     accessory_gene_count,
                     low_freq_genes_in_region,
                     core_gene_pair_distance,
                     accessory_gene_content,
                     low_freq_gene_content,
                     core_gene_pairs,
                     master_info) = record_core_core_region(core_genes, gff_name, line, previous_core_gene_id,
                                                            previous_core_gene_end_coor, accessory_gene_count,
                                                            low_freq_genes_in_region, core_gene_pair_distance,
                                                            accessory_gene_content, low_freq_gene_content,
                                                            core_gene_pairs, master_info)

            else:
                # Check if accessory is low frequency - else just regular 'regular' accessory
                if line[8] in low_freq_genes[gff_name]:
                    low_freq_genes_in_region.append(low_freq_genes[gff_name][line[8]])
                    accessory_gene_count = + 1
                else:
                    accessory_gene_count += 1

        else: # TODO - FIX iF there is no core gene on first contig!
            # Set new contig
            previous_contig = line[0]

            # Note that gff has multiple contigs
            single_contig = False

            # Check if there is core gene on traversed contig or if first contig and there is no core gene,
            # if then skip.
            if previous_core_gene_id is not "Sequence_break" and previous_core_gene_id is not "":
                # Record the core gene neighbouring a sequence break
                (previous_core_gene_id,
                 previous_core_gene_end_coor,
                 accessory_gene_count,
                 low_freq_genes_in_region,
                 core_gene_pair_distance,
                 accessory_gene_content,
                 low_freq_gene_content,
                 core_gene_pairs,
                 master_info) = record_core_core_region(core_genes, gff_name, None, previous_core_gene_id,
                                                        previous_core_gene_end_coor, accessory_gene_count,
                                                        low_freq_genes_in_region, core_gene_pair_distance,
                                                        accessory_gene_content, low_freq_gene_content,
                                                        core_gene_pairs, master_info)

    # Check if genome is complete or a single contig. If then add information for last and first core gene, if not
    # then add the first and last core gene as being neighbours to sequence breaks.
    if single_contig:
        last_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
        first_core_gene_cluster = core_genes[gff_name][first_core_gene_id]

        # Add core neighbours
        core_gene_neighbours = sorted([last_core_gene_cluster, first_core_gene_cluster])
        core_gene_neighbours = f'{core_gene_neighbours[0]}-{core_gene_neighbours[1]}'
        core_gene_pairs.append(core_gene_neighbours)

        # Add core neighbour distance
        core_gene_pair_distance[core_gene_neighbours] = [get_genome_size_from_gff(gff_path) \
                                                        - previous_core_gene_end_coor + first_core_gene_start_coor]

        # Add accessory information from between last and first core gene
        accessory_gene_content[core_gene_neighbours] = accessory_gene_count + first_core_accessory_count
        low_freq_gene_content[core_gene_neighbours] = low_freq_genes_in_region + first_core_low_freq_genes

    else:
        # Add first core gene as being neighbour to a sequence break
        (previous_core_gene_id,
         previous_core_gene_end_coor,
         accessory_gene_count,
         low_freq_genes_in_region,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content, core_gene_pairs, master_info) = record_core_core_region(core_genes, gff_name, None,
                                                                                        first_core_gene_id,
                                                                                        first_core_gene_start_coor,
                                                                                        first_core_accessory_count,
                                                                                        first_core_low_freq_genes,
                                                                                        core_gene_pair_distance,
                                                                                        accessory_gene_content,
                                                                                        low_freq_gene_content,
                                                                                        core_gene_pairs, master_info)
        # Add last core gene as being neighbour to a sequence break
        # if the last core gene has not been recorded already
        if previous_core_gene_id is not "Sequence_break":
            (previous_core_gene_id,
             previous_core_gene_end_coor,
             accessory_gene_count,
             low_freq_genes_in_region,
             core_gene_pair_distance,
             accessory_gene_content,
             low_freq_gene_content, core_gene_pairs, master_info) = record_core_core_region(core_genes, gff_name, None,
                                                                               previous_core_gene_id,
                                                                               previous_core_gene_end_coor,
                                                                               accessory_gene_count,
                                                                               low_freq_genes_in_region,
                                                                               core_gene_pair_distance,
                                                                               accessory_gene_content,
                                                                               low_freq_gene_content,
                                                                               core_gene_pairs, master_info)

    return core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info


def segment_genome_content(input_file, core_genes, low_freq_genes, i):
    """ Single function segmenting the gff into core gene regions to be used for simple multi processing"""
    if (i+1) % 25 == 0 or i == 0:
        print(f"Processing GFF file #{i+1}")

    gff_generator = parse_gff(input_file)
    return_data = segment_gff_content(gff_generator=gff_generator, gff_path=input_file, core_genes=core_genes, low_freq_genes=low_freq_genes)

    return return_data
