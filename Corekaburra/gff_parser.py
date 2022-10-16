import os
import gzip


def open_file_generator(input_file_path):
    """
    Function to read a gff file as either gzip or text
    :param input_file_path: String path to file
    :return: Generator with line by line
    """
    # Open input file as if it was gzipped
    try:
        with gzip.open(input_file_path, 'rt') as open_file:
            # Test if gzipped by reading line
            open_file.readline()

            for line in open_file:
                yield line

    except (OSError, gzip.BadGzipFile):
        # Open input as if normal
        with open(input_file_path, 'r') as open_file:
            for line in open_file:
                yield line


def parse_gff(input_file):
    """
    Try to read a GFF file as gzipped then normal
    pass it to the parser and return it as a generator object that return all line containing CDS.
    Break whenever the fasta sequence is reached.
    :param input_file: File-path to a given gff file to be processed
    :return: line from generator object returning CDS from a gff file
    """
    file_generator = open_file_generator(input_file)

    for line in file_generator:
        if "##FASTA" in line:
            # FASTA found - end loop
            file_generator.close()
            break
        if "#" not in line and ('CDS' in line or 'candidate_gene' in line):
            # Strip line for newline and split columns into list
            line = line.strip()
            line = line.split("\t")

            # See if refound gene or Prokka annotated and isolate ID in gene_presence_absence.csv accordingly
            if "old_locus_tag=" in line[8]:
                gene_id = line[8][line[8].find('old_locus_tag'):]
                if ';' in gene_id:
                    gene_id = gene_id[:gene_id.find(';')]
            else:
                gene_id = line[8][line[8].find('ID'):line[8].find(';')]

            # Remove equal sign from id and add as identifier for the returned gff line
            gene_id = gene_id[gene_id.find('=') + 1:]
            line[8] = gene_id
            yield line


def get_contig_lengths(input_file):
    """
    Function that takes an input gff file path and records the length of each contig in the file
    :param input_file: File path for a given gff file
    :return: directory with key being the contig name (before first white spae)
    and value the size of the contig in base pairs
    """
    fasta_reached = False
    contig_size = 0
    contig_size_dir = {}

    # Open the given gff file, find the fasta section of the file and count the length of each contig.
    for line in open_file_generator(input_file):
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
                            low_freq_gene_content, core_gene_pairs, master_info):
    """
    Function to record information about a core gene pair or a core gene and a sequence break,
    along with accessory information between the two features
    :param core_genes: Dict of core genes passed to genomes and the pan-genome clusters.
    :param gff_name: Name of the gff file currently being examined
    :param gff_line: List for the gene currently being recorded. (Can be None, when previous gene is next to sequence break)
    :param contig_end: Length of the contig in question, used to calculate the distance from last gene to end of contig
    :param previous_core_gene_id: ID of the previous gene recorded in pair, can be Sequence_break
    :param previous_core_gene_end_coor: End coordinate of the previously recorded genome (Can be None when new contig is initiated)
    :param acc_genes_in_region: List of recorded accessory-frequency genes in the region being recorded
    :param low_freq_genes_in_region: List of recorded low-frequency genes in the region being recorded
    :param core_gene_pair_distance: Dict of distances between core pairs recorded. key being the pair and value the distance in base-pairs
    :param accessory_gene_content: Dict of accessory frequency genes and their mapping to genomes
    :param low_freq_gene_content: Dict of low frequency genes and their mapping to genomes
    :param core_gene_pairs: List of core pairs recorded
    :param master_info: Large dict holding key for each pair recorded with a list of info around the pair as value
    :return: A tuple of dictionaries and lists, most of which are also given as input - See descriptions above
    """

    # Check that a line from gff is provided and previous gene is not a sequence break
    if gff_line is not None and previous_core_gene_id != "Sequence_break":
        # Check if core gene is fragmented, if then change coordinates to the last part of the fragment.
        if core_genes[gff_name][previous_core_gene_id] == core_genes[gff_name][gff_line[8]]:
            previous_core_gene_id = gff_line[8]
            previous_core_gene_end_coor = int(gff_line[4])
            acc_genes_in_region = []
            low_freq_genes_in_region = []
            return (previous_core_gene_id, previous_core_gene_end_coor, acc_genes_in_region, low_freq_genes_in_region,
                    core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, core_gene_pairs, master_info)

    # Set core cluster names
    # If no line from gff is given there is a sequence-break,
    # if it is given then set current cluster and try to find previous if not found it is a sequence break
    if gff_line is not None:
        current_core_gene_cluster = core_genes[gff_name][gff_line[8]]
        try:
            previous_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
            core_gene_neighbours = sorted([previous_core_gene_cluster, current_core_gene_cluster])
        # Catch is previous gene was a sequence break.
        except KeyError:
            previous_core_gene_cluster = previous_core_gene_id
            core_gene_neighbours = [previous_core_gene_cluster, current_core_gene_cluster]

        # core_gene_neighbours = sorted([previous_core_gene_cluster, current_core_gene_cluster])

    else:
        current_core_gene_cluster = "Sequence_break"
        previous_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
        core_gene_neighbours = [previous_core_gene_cluster, current_core_gene_cluster]

    # Merge core neighbours and add to pairs
    core_gene_neighbours_str = f'{core_gene_neighbours[0]}--{core_gene_neighbours[1]}'
    core_gene_pairs.append(core_gene_neighbours_str)

    # Set core neighbour distance
    # Check if measuring between two genes on same contig and not measuring from a sequence break
    if gff_line is not None and previous_core_gene_cluster != "Sequence_break":
        core_core_distance = int(gff_line[3]) - previous_core_gene_end_coor - 1
    else:
        # Check if measuring between sequence break and first core gene on contig, if then set start coordinate to zero
        if previous_core_gene_id == "Sequence_break":
            core_core_distance = int(gff_line[3]) - 1
            # If gene starts at contig end, make distance zero
            if core_core_distance < 0:
                core_core_distance = 0
        elif current_core_gene_cluster == 'Sequence_break':
            core_core_distance = abs(contig_end - previous_core_gene_end_coor)
        else:
            NotImplementedError(
                'An error occured when measuring the distance between core gene and contig end. Something went wrong!')

    # Add core neighbour distance
    core_gene_pair_distance[core_gene_neighbours_str] = core_core_distance

    # Add counts and annotation for accessory and low frequency genes
    acc_genes_in_region = sorted(list(set(acc_genes_in_region)))
    low_freq_genes_in_region = sorted(list(set(low_freq_genes_in_region)))

    accessory_gene_content[core_gene_neighbours_str] = acc_genes_in_region.copy()
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
            core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, core_gene_pairs, master_info)


def connect_first_n_last_gene_on_contig(core_genes, gff_name, previous_core_gene_id, previous_core_gene_end_coor,
                                        first_core_gene_gff_line, acc_genes_in_region, first_core_accessory_content,
                                        low_freq_genes_in_region, first_core_low_freq_genes, contig_size,
                                        core_gene_pairs, core_gene_pair_distance,
                                        accessory_gene_content, low_freq_gene_content, master_info):
    """
    Function to record the connection between the first and the last genome on a closed contig

    :param core_genes: A dict of dicts mapping genomes to gene IDs to pan-genome clusters of core genes.
    :param gff_name: Name of the gff file currently being examined
    :param previous_core_gene_id: ID of the previous gene recorded in pair, can be Sequence_break.
    :param previous_core_gene_end_coor: End coordinate of the previously recorded genome.
    :param first_core_gene_gff_line: List for the gene first encountered on the contig.
    :param acc_genes_in_region: List of recorded accessory-frequency genes in the region after last core genes
    :param first_core_accessory_content: List of recorded accessory-frequency genes in the region before first core genes
    :param low_freq_genes_in_region: List of recorded low-frequency genes in the region after last core genes
    :param first_core_low_freq_genes: List of recorded low-frequency genes in the region before first core genes
    :param contig_size: Size of the contig currently being looked at.
    :param core_gene_pairs: List of core pairs recorded
    :param core_gene_pair_distance: Dict of distances between core pairs recorded. key being the pair and value the distance in base-pairs
    :param accessory_gene_content: Dict of accessory frequency genes and their mapping to genomes
    :param low_freq_gene_content: Dict of low frequency genes and their mapping to genomes
    :param master_info: Large dict holding key for each pair recorded with a list of info around the pair as value

    :return previous_core_gene_id: String of the last core genes ID/locus_tag
    :return previous_core_gene_end_coor: Int of the end coordinate for the latest core gene
    :return acc_genes_in_region: Empty list to store accessory genes
    :return low_freq_genes_in_region: Empty list to store low-frequency genes
    :return core_gene_pairs: List of core gene pairs recorded
    :return core_gene_pair_distance: A dict of distances between a specific pair of core genes.
    :return accessory_gene_content: A dict of the accessory genes found between a pair of core genes.
    :return low_freq_gene_content: A dict of the low-frequency genes found between a pair of core genes.
    :return master_info: A dict of multiple pieces of info for each core gene pair.
    """
    last_core_gene_cluster = core_genes[gff_name][previous_core_gene_id]
    first_core_gene_cluster = core_genes[gff_name][first_core_gene_gff_line[8]]


    # Add core neighbours
    core_gene_neighbours = sorted([last_core_gene_cluster, first_core_gene_cluster])
    core_gene_neighbours_str = f'{core_gene_neighbours[0]}--{core_gene_neighbours[1]}'
    core_gene_pairs.append(core_gene_neighbours_str)

    # Add core neighbour distance
    core_core_distance = contig_size - previous_core_gene_end_coor + int(first_core_gene_gff_line[3]) - 1

    core_gene_pair_distance[core_gene_neighbours_str] = core_core_distance

    # Add accessory information from between last and first core gene
    last_first_accessory_content = acc_genes_in_region.copy() + first_core_accessory_content.copy()
    last_first_low_freq_count = low_freq_genes_in_region.copy() + first_core_low_freq_genes.copy()

    last_first_accessory_content = sorted(list(set(last_first_accessory_content)))
    last_first_low_freq_count = sorted(list(set(last_first_low_freq_count)))

    accessory_gene_content[core_gene_neighbours_str] = last_first_accessory_content.copy()
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

    previous_core_gene_id = ""
    previous_core_gene_end_coor = int(first_core_gene_gff_line[4])
    acc_genes_in_region = []
    low_freq_genes_in_region = []

    return (previous_core_gene_id, previous_core_gene_end_coor, acc_genes_in_region, low_freq_genes_in_region,
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info)


def record_coreless_contig(coreless_contigs, acc_genes_in_region, low_freq_genes_in_region, gff_name, contig_name):
    """
    Function to record the presence and info of a contig without a core gene
    :param coreless_contigs: List with info of other core less contigs
    :param acc_genes_in_region: List of accessory genes on contig
    :param low_freq_genes_in_region: List of low frequency genes on contig
    :param gff_name: String name of input gff file.
    :param contig_name: String name of the contg

    :return: Dict of information round coreless contigs as list.
    """
    acc_genes_in_region = sorted(list(set(acc_genes_in_region)))
    low_freq_genes_in_region = sorted(list(set(low_freq_genes_in_region)))
    if len(acc_genes_in_region) + len(low_freq_genes_in_region) > 0:
        coreless_contigs[f'{gff_name}--{contig_name}'] = [acc_genes_in_region, low_freq_genes_in_region]

    return coreless_contigs


def segment_gff_content(gff_generator, core_genes, low_freq_genes, gff_path, acc_genes, complete_genomes):
    """
    Function that takes a gff generator, core, accessory, and low frequency genes and identify each core-core gene region
    counts the number of accessory genes in the region, records the number of low frequency genes in the region,
    and records the distance from one core gene to the next.
    :param gff_generator: A generator object providing each CDS line from a gff file to be segmented.
    :param core_genes: A dict of dicts mapping genomes to gene IDs to pan-genome clusters of core genes.
    :param low_freq_genes: Same structure as core_genes, but for low-frequency genes.
    :param gff_path: List of file paths to gff files.
    :param acc_genes: Same structure as core_genes, but for accessory genes.
    :param complete_genomes: List of gff names that are to be handled as complete.

    :return core_gene_pairs: A list of core genes (or sequence breaks) that are found to be neighbouring each other in a given gff file.
    :return core_gene_pair_distance: A dict of distances between a specific pair of core genes.
    :return accessory_gene_content: A dict of the accessory genes found between a pair of core genes.
    :return low_freq_gene_content: A dict of the low-frequency genes found between a pair of core genes.
    :return master_info: A dict of multiple pieces of info for each core gene pair (Gff file, core gene 1, core gene 2, distnace between them, genes between them, list of accesspry genes, list of low-frequency genes)
    :return coreless_contigs: Dict of contigs found to not encode any core genes on them. The accessory and low-frequency genes are recorded.
    """
    # Initialize data structures to be returned
    core_gene_pairs = []
    core_gene_pair_distance = {}
    low_freq_gene_content = {}
    accessory_gene_content = {}
    master_info = {}
    coreless_contigs = {}

    # Examine if a complete genome has been given
    if complete_genomes is None:
        complete_genome = False
    else:
        if os.path.basename(gff_path).replace('.gz', '').replace('.gff', '') in complete_genomes:
            complete_genome = True
        else:
            complete_genome = False

    # Get size of contigs for given gff
    contig_sizes = get_contig_lengths(gff_path)

    # Split input path of gff to get genome name
    gff_name = gff_path.split('/')[-1]
    gff_name = gff_name.replace('.gz', '').rsplit('.', 1)[0]
    if 'corrected' in gff_name:
        gff_name = gff_name.split('_corrected')[0]

    # Set that first core gene has not been found
    first_core_gene = True
    # Set Initialise variable for the ID of previous encountered core gene
    previous_core_gene_id = ""
    # Set variable to get first contig
    first_contig = True
    # Initialise the accessory and low-frequency gene lists
    low_freq_genes_in_region = []
    acc_genes_in_region = []

    # Go through each line of GFF file
    for line in gff_generator:
        # Set first contig to fund in file
        if first_contig:
            previous_contig = line[0]
            first_contig = False

        # Check if contig has changed - if then finish contig, if not examine next gene on contig
        if line[0] == previous_contig:
            # Check if gene is core, if then test if it is the first, else assume to be accessory gene
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
                    low_freq_genes_in_region = []
                    acc_genes_in_region = []

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
                     master_info) = record_core_core_region(core_genes, gff_name, line, None, previous_core_gene_id,
                                                            previous_core_gene_end_coor, acc_genes_in_region,
                                                            low_freq_genes_in_region, core_gene_pair_distance,
                                                            accessory_gene_content, low_freq_gene_content,
                                                            core_gene_pairs, master_info)

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
            # Check if there is a core gene on traversed contig or if a core gene is present on the first contig -
            # if then record it, if not record the accessory and low frequency genes found on contig and reset.
            if previous_core_gene_id != "Sequence_break" and previous_core_gene_id != "":
                if complete_genome:
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
                                                                        first_core_low_freq_genes, contig_sizes[previous_contig],
                                                                        core_gene_pairs, core_gene_pair_distance,
                                                                        accessory_gene_content, low_freq_gene_content,
                                                                        master_info)

                    first_core_gene = True

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
                     master_info) = record_core_core_region(core_genes, gff_name, None, contig_sizes[previous_contig],
                                                            previous_core_gene_id, previous_core_gene_end_coor,
                                                            acc_genes_in_region, low_freq_genes_in_region,
                                                            core_gene_pair_distance, accessory_gene_content,
                                                            low_freq_gene_content, core_gene_pairs, master_info)
            else:
                # Record info on accessory genes on core-less contig, if any accessory genes are present
                coreless_contigs = record_coreless_contig(coreless_contigs, acc_genes_in_region,
                                                          low_freq_genes_in_region, gff_name, previous_contig)

                # Reset accessory and low-frequency gene lists
                acc_genes_in_region = []
                low_freq_genes_in_region = []

                # Set new contig
                previous_contig = line[0]

            # Check if fist gene is core on complete genome, if then record details.
            if line[8] in core_genes[gff_name] and complete_genome:
                # Set information on first core gene to be used when finishing search
                first_core_gene_gff_line = line

                # Set information to be used with next core gene neighbour
                previous_core_gene_end_coor = int(line[4])
                previous_core_gene_id = line[8]

                # Set that first core gene has been observed
                first_core_gene = False

            # Check if first gene on new contig is a core gene, if then record it.
            elif line[8] in core_genes[gff_name]:
                previous_core_gene_id = "Sequence_break"

                # Get the starting position of the first core gene on contig to record the gene.
                # Make it negative to fit the calculation of the distance between genes.
                cur_core_gene_start = -int(line[3])

                (previous_core_gene_id,
                 previous_core_gene_end_coor,
                 acc_genes_in_region,
                 low_freq_genes_in_region,
                 core_gene_pair_distance,
                 accessory_gene_content,
                 low_freq_gene_content,
                 core_gene_pairs,
                 master_info) = record_core_core_region(core_genes, gff_name, line, 0, previous_core_gene_id,
                                                        cur_core_gene_start, acc_genes_in_region,
                                                        low_freq_genes_in_region, core_gene_pair_distance,
                                                        accessory_gene_content, low_freq_gene_content,
                                                        core_gene_pairs, master_info)

            # Add as accessory - if first gene is not core
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

            # Set new contig
            previous_contig = line[0]

    # Check if genome is complete or a single contig. If then add information for last and first core gene, if not
    # then add the first and last core gene as being neighbours to sequence breaks.
    if complete_genome:
        if previous_core_gene_id != "":
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
            record_coreless_contig(coreless_contigs, acc_genes_in_region, low_freq_genes_in_region, gff_name, line[0])

    else:
        # Add first core gene as being neighbour to a sequence break
        (_,
         _,
         _,
         _,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content, core_gene_pairs,
         master_info) = record_core_core_region(core_genes, gff_name, first_core_gene_gff_line, 0,
                                                                         'Sequence_break',
                                                                         9999999999999999,
                                                                         first_core_accessory_content,
                                                                         first_core_low_freq_genes,
                                                                         core_gene_pair_distance,
                                                                         accessory_gene_content,
                                                                         low_freq_gene_content,
                                                                         core_gene_pairs,
                                                                         master_info)
        # Add last core gene as being neighbour to a sequence break
        # if the last core gene has not been recorded already
        if previous_core_gene_id != "Sequence_break":
            (previous_core_gene_id,
             previous_core_gene_end_coor,
             acc_genes_in_region,
             low_freq_genes_in_region,
             core_gene_pair_distance,
             accessory_gene_content,
             low_freq_gene_content, core_gene_pairs,
            master_info) = record_core_core_region(core_genes, gff_name, None, contig_sizes[previous_contig],
                                                   previous_core_gene_id,
                                                   previous_core_gene_end_coor,
                                                   acc_genes_in_region,
                                                   low_freq_genes_in_region,
                                                   core_gene_pair_distance,
                                                   accessory_gene_content,
                                                   low_freq_gene_content,
                                                   core_gene_pairs,
                                                   master_info)
        else:
            # Add a core-less contig if there has been accessory genes:
            coreless_contigs = record_coreless_contig(coreless_contigs, acc_genes_in_region,
                                                      low_freq_genes_in_region, gff_name, line[0])

    return core_gene_pairs, core_gene_pair_distance, accessory_gene_content, \
           low_freq_gene_content, master_info, coreless_contigs


def segment_genome_content(input_gff_file, core_genes, low_freq_genes, acc_gene_dict, complete_genomes):
    """
    Single function segmenting the gff into core gene regions to be used for simple multi processing
    :param input_gff_file: File-path to the given gff file to be segmented
    :param core_genes: Dictionary over core genes
    :param low_freq_genes: Dictionary over low-frequency genes
    :param acc_gene_dict: Dictionary over accessory genes
    :param complete_genomes: Bool indicating if this genome should be considered as a complete genome

    :return input_gff_file: File path to the gff being searched
    :return core_genes: Dict of core genes passed to genomes and the pan-genome clusters.
    :return low_freq_genes: Same structure as core_genes, but for low-frequency genes.
    :return acc_gene_dict: Same structure as core_genes, but for accessory genes.
    :return i: The index of the gff in the larger scheme of the analysis
    :return complete_genomes: List of genomes given as complete by the user.
    """

    gff_generator = parse_gff(input_gff_file)
    return_data = segment_gff_content(gff_generator=gff_generator,
                                      gff_path=input_gff_file,
                                      core_genes=core_genes,
                                      low_freq_genes=low_freq_genes,
                                      acc_genes=acc_gene_dict,
                                      complete_genomes=complete_genomes)

    return return_data
