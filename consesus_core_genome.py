import networkx as nx
import numpy as np
from os.path import basename


def construct_core_graph(core_neighbour_pairs):
    # Initiate core gene graph
    G = nx.Graph()

    # Add all core gene pairs and their edges
    for core_set in core_neighbour_pairs.keys():
        # split core genes
        core_genes = core_set.split('--')

        # Check that sequence break is not present:
        if 'Sequence_break' != core_genes[0] and 'Sequence_break' != core_genes[1]:
            # Construct edge list
            core_genes = [(core_genes[0], core_genes[1], {'weight': core_neighbour_pairs[core_set]})]
            # Add edge
            G.add_edges_from(ebunch_to_add=core_genes)

    # TODO - Return graph with nodes of core genes and weight being the number of times observed to be connected
    return G


def get_most_frequent_gene(gene_list):
    genes, frequency = np.unique(gene_list, return_counts=True)
    max_index = np.where(frequency == np.amax(frequency))
    most_frequent_gene = genes[max_index][0]

    return most_frequent_gene


def get_new_neighbours(current_gene, visited_genes, core_gene_graph):
    """ Function to get new neighbours and their weights and filter away already visited neighbours """
    neighbours = core_gene_graph.neighbors(current_gene)

    new_neighbours = [neighbour for neighbour in neighbours if neighbour not in visited_genes]

    edge_weights = [core_gene_graph.get_edge_data(current_gene, neighbour)['weight'] for neighbour in new_neighbours]

    return new_neighbours, edge_weights


def find_core_gene_synteny(core_gene_graph, start_gene_clusters, second_gene_cluster):
    """ Function to identify core gene synteny from core gene graph
     Also returns a set of genes that possibly are associated with rearrangements"""
    # Check that at least one start gene is found, if then get most frequent gene in first and second position.
    if len(start_gene_clusters) > 0:
        start_gene = get_most_frequent_gene(start_gene_clusters)
        second_gene = get_most_frequent_gene(second_gene_cluster)

        print(f"Most frequent first gene in complete genomes and "
              f"chosen as starting gene for consensus core gene synteny is: {start_gene}")
        print(f"Most frequent second gene in complete genes and used to inform core gene synteny is: {second_gene}")

        # Get place to start core gene synteny
        visited_list = [start_gene, second_gene]
        new_neighbours, edge_weights = get_new_neighbours(second_gene, visited_list, core_gene_graph)

        # Initialise list to hold the weights for core-core edges.
        # First core gene
        _, edge_weights_first = get_new_neighbours(second_gene, [], core_gene_graph)
        _, edge_weights_second = get_new_neighbours(second_gene, [start_gene], core_gene_graph)
        core_path_coverage = [[edge_weights_first[0]], edge_weights_second]


    # TODO - Else, when only contigged genomes are pressent (Search for dnaA? first naively using regex in list then more sofisticated)

    # Initialise list to store genes of interest and narrow the search for rearrangements
    possible_rearrangement_genes = []

    # Search Graph for path
    while len(new_neighbours) > 0:
        # Check if more than one neighbour is present in, then chose the one with largest edge weight,
        # Else walk along only path available
        if len(new_neighbours) > 1:
            max_weight = np.max(edge_weights)
            largest_weight_index = np.where([index == max_weight for index in edge_weights])[0][0]
            best_neighbour = new_neighbours[largest_weight_index]

        else:
            best_neighbour = new_neighbours[0]

        # Add next node in core genome synteny
        visited_list.append(best_neighbour)
        new_neighbours, edge_weights = get_new_neighbours(best_neighbour, visited_list, core_gene_graph)
        core_path_coverage = core_path_coverage + [edge_weights]
        # TODO - Use core_path_coverage to produce a plot that show the how much coverage each edge has.

        # Check if more than one neighbour is available due to possible rearrangements
        if len(new_neighbours) > 1:
            # Map possible combinations of branch in core-core genome path
            possible_pairs = list(map(lambda e: [best_neighbour, e], new_neighbours))
            possible_rearrangement_genes = possible_rearrangement_genes + possible_pairs

    # Insert connections from last gene to first gene
    _ , core_path_coverage[-1] = get_new_neighbours(visited_list[-1], [visited_list[-2]], core_gene_graph)
    return visited_list, possible_rearrangement_genes, core_path_coverage


def determine_core_gene_consesus(core_neighbour_pairs, start_gene_cluster, second_gene_cluster):
    core_gene_graph = construct_core_graph(core_neighbour_pairs)

    consensus_core_genome, \
    possible_rearrangement_genes, \
    core_path_coverage = find_core_gene_synteny(core_gene_graph,
                                                start_gene_cluster,
                                                second_gene_cluster)

    # TODO - Guess size of rearrangement using distance from genome and fill in with average, if missing.

    return consensus_core_genome, possible_rearrangement_genes, core_path_coverage


def assign_core_gene_synteny_types(alternative_core_pairs, gff_names):
    """ Function that assigns a type to each genome based on its core genome synteny and the consensus.
     returns a dict with each genome as a key and the value being the consensus core genome type.
     a type = 1 is consensus everything else is not consensus."""
    # Get name of genomes
    gff_names = [basename(gff).split('.')[0] for gff in gff_names]

    # construct dict with each genome as key and value as core gene synteny type
    type_dict = dict.fromkeys(gff_names)

    # Get list of all genomes with an alternative core pair
    alternative_genomes = alternative_core_pairs.keys()

    # Construct dict to hold all alternative core-pair and their associated type
    alt_core_comp_types = {}

    # Initialise alternative type counter
    type_counter = 2
    # Go through all genomes and assign core gene synteny type dependant on their composition
    for key in type_dict.keys():
        # Check if alternative core gene neighbours are present,
        # If they give search for type, if not then give consencus type, 1.
        if key in alternative_genomes:
            # extract alternative core neighbours and sort and make string to be possible key in dict
            current_pairs = alternative_core_pairs[key]
            current_pairs.sort()
            current_pairs = str(current_pairs)

            # Check if combination of alternative core gene neighbours has been assigned a type,
            # if not assign new type.
            if current_pairs in list(alt_core_comp_types.keys()):
                type_dict[key] = alt_core_comp_types[current_pairs]
            else:
                alt_core_comp_types[str(alternative_core_pairs[key])] = type_counter
                type_dict[key] = type_counter
                type_counter = type_counter + 1
        else:
            type_dict[key] = 1

    return type_dict, alt_core_comp_types


def identify_rearrangements(consensus_core_genome, possible_rearrangement_genes, master_info_dict, gff_names):
    # Pair all neighbouring genes in the consensus core genome
    core_genome_pairs = []
    for i in range(len(consensus_core_genome)):
        sorted_neighbours = sorted([consensus_core_genome[i], consensus_core_genome[i-1]])
        core_genome_pairs.append(f'{sorted_neighbours[0]}--{sorted_neighbours[1]}')

    # Make searchable set of all core gene pairs
    core_genome_pairs = set(core_genome_pairs)

    # TODO!
    # Search each key that contain a gene with possible rearrangements.
    # If two core gene pairs is not present then predict rearrangement and size. - store genome and rearrangements
    # If only one is present predict possible rearrangements, but unsure due to only one core pair identified
    # If three or more then it is complex and no further prediction can be made. - Possible contamination?

    # initialize dict of alternative core gene neighbours
    alt_core_pairs = {}

    # Identify the alternative core gene neighbours for each genome:
    for key in master_info_dict.keys():
        split_key = key.split("--")
        if "Sequence_break" != split_key[0] and "Sequence_break" != split_key[1]:
            core_pair = {f"{split_key[0]}--{split_key[1]}"}

            difference = core_pair.difference(core_genome_pairs)
            if len(difference):
                try:
                    alt_core_pairs[split_key[2]] = alt_core_pairs[split_key[2]] + list(difference)
                except KeyError:
                    alt_core_pairs[split_key[2]] = list(difference)


    print(f'Number of genomes with alternative core neighbours {len(alt_core_pairs.keys())} - '
          f'{round(len(alt_core_pairs.keys())/len(gff_names) * 100, ndigits=1)}%\n')

    # Record the number of times each alternative core pair occur
    alt_core_pair_count = {}
    for key in alt_core_pairs.keys():
        for pair in alt_core_pairs[key]:
            if pair in list(alt_core_pair_count.keys()):
                alt_core_pair_count[pair] += 1
            else:
                alt_core_pair_count[pair] = 1

    # Assign core gene synteny types to genomes
    core_genome_types, alt_core_comp_types = assign_core_gene_synteny_types(alt_core_pairs, gff_names)

    return alt_core_pairs, alt_core_pair_count, core_genome_types, alt_core_comp_types


def simple_rearrangement_prediction(gene_pairs, consensus_core_genome):
    print("predicting rearrangements")
    # Pair all neighbouring genes in consensus sequence
    core_genome_pairs = []
    for i in range(len(consensus_core_genome)):
        sorted_neighbours = sorted([consensus_core_genome[i], consensus_core_genome[i - 1]])
        core_genome_pairs.append(f'{sorted_neighbours[0]}--{sorted_neighbours[1]}')


    # TODO - Check if alternative neighbours are neighbours in the concensus if paired differently
    individual_genes = gene_pairs[0].split('--') + gene_pairs[1].split('--')

    new_gene_pairs = [f'{individual_genes[0]}--{individual_genes[2]}',
                      f'{individual_genes[1]}--{individual_genes[3]}']

    if new_gene_pairs[0] in core_genome_pairs and new_gene_pairs[1] in core_genome_pairs:
        print("Solution found")
    else:
        new_gene_pairs = [f'{individual_genes[1]}--{individual_genes[2]}',
                          f'{individual_genes[0]}--{individual_genes[3]}']

        if new_gene_pairs[0] in consensus_core_genome and new_gene_pairs[1] in consensus_core_genome:
            print("Solution found")


def characterise_rearrangements(alt_core_pairs, consensus_core_genome):
    print("characterising rearrangements")
    rearrangement_predictions = [[], []]

    # Go through each genome with alternative core pairs
    for genome in alt_core_pairs.keys():

        # If only one alternative pair is found, too little info is available. (May be due to contig breaks)
        if len(alt_core_pairs[genome]) == 1:
            rearrangement_predictions[0].append(genome)
            rearrangement_predictions[1].append('Too little information')

        # If exactly two alternative core pairs are present it is possible to make a prediction, with some uncertainty.
        if len(alt_core_pairs[genome]) == 2:
            rearrangement_predictions[0].append(genome)
            rearrangement_predictions[1].append('Possible prediction')

            simple_rearrangement_prediction(alt_core_pairs[genome], consensus_core_genome)

        # If more than two are present then more than one solution is possible and no prediction can be made.
        if len(alt_core_pairs[genome]) > 2:
            rearrangement_predictions[0].append(genome)
            rearrangement_predictions[1].append('Possible prediction')
            # TODO - Look at predicting more complex rearrangements.

    '''Use consensus core synteny to work out how rearrangements may have occured. 
    If A and B are connected in the concesus and C and D
    are connected, then if A and C are connected, and B and D are found next to 
    sequence breaks then it may be a possible recombination'''


    return rearrangement_predictions # TODO - Use this output


def determine_partners_neighbours(alt_core_gene, core_genome_pairs, cur_genome_pairs):
    # Find consensus neighbours for genes in the alternative pair
    consensus_pairs = [gene for gene in core_genome_pairs if gene[0] is alt_core_gene or gene[1] is alt_core_gene]

    # Remove consensus pairs found in the current genome
    consensus_pairs = [pair for pair in consensus_pairs if pair not in cur_genome_pairs]

    # Flatten list of list to list
    consensus_neighbours = [gene for pair in consensus_pairs for gene in pair]

    # Find remove the gene from the alternative pair and keep consensus neighbours
    consensus_neighbours = [gene for gene in consensus_neighbours if gene is not alt_core_gene]

    # Find neighbours of the consesus neighbours in the current genome
    cur_consensus_neighbours = [pair for pair in cur_genome_pairs if pair[1] in consensus_neighbours
                               or pair[0] in consensus_neighbours]

    # Flatten the list of neighbours to the consensus neighbours
    consensus_neighbours = [gene for pair in cur_consensus_neighbours for gene in pair]
    # print(consensus_neighbours)

    # Count the number of sequence breaks
    partner_sequence_breaks = consensus_neighbours.count('Sequence_break')

    return partner_sequence_breaks


def core_pair_matrix(core_genome_types, alt_core_comp_types, alt_core_pair_count, master_info_total, consensus_genome):
    """ Function to produce a matrix of presence or absence of alternative core neighbours to be viewed in Phandango"""
    # Construct header for output file
    header = ['genome']
    for pair in alt_core_pair_count.keys():
        header.append(pair)

    # Initialise the dict that will become the matrix with gene pairs as coluns and genomes as rows.
    alt_core_pair_matrix = {}
    for genome in core_genome_types:
        alt_core_pair_matrix[genome] = dict.fromkeys(header[1:], 0)
        alt_core_pair_matrix[genome].update({'genome': genome})

    # Construct dict that contain the genome type as key and alternative core neghbours associated with it as values.
    genome_type_gene_dict = {}
    for gene_pairs in alt_core_comp_types:
        core_gene_type = alt_core_comp_types[gene_pairs]
        gene_pairs = gene_pairs.replace('[', '')
        gene_pairs = gene_pairs.replace(']', '')
        gene_pairs = gene_pairs.split(', ')

        gene_pairs = [eval(gene_pair) for gene_pair in gene_pairs]

        genome_type_gene_dict[core_gene_type] = gene_pairs

    # Fill the matrix with presence (3) or leave as absent (0)
    for key in core_genome_types:
        genome_type = core_genome_types[key]
        if genome_type is not 1:
            for core_pair in genome_type_gene_dict[genome_type]:
                alt_core_pair_matrix[key][core_pair] = 3

    # TODO - Go though and see of an isolate contain a sequence break by both genes in an alternative pair.
    #  if so then mark it as 0.5 or probable.
    # Find core genes where both genes are next to a sequence break and can be found in a alternative pair
    sequence_break_pairs = [pair.split('--') for pair in master_info_total.keys() if 'Sequence_break' in pair]
    # Remove the Sequence_break enteries. Keep gene neighbouring sequnece break and genome name
    sequence_break_pairs = [[element for element in pair if element != 'Sequence_break'] for pair in sequence_break_pairs]

    # Pair all neighbouring genes in the consensus core genome
    core_genome_pairs = []
    for i in range(len(consensus_genome)):
        sorted_neighbours = sorted([consensus_genome[i], consensus_genome[i - 1]])
        core_genome_pairs.append(sorted_neighbours)


    # TODO - predict if a consensus - core-core neighbourship is possible.
    # Initiate counters to keep track of core-core variants searched and found
    core_variants_predicted_weak = 1
    core_variants_predicted_strong = 1
    core_variants_predicted = 1
    genomes_searched = 1
    # Go through each header and find genomes that contain both core genes of a alternative pair next to a sequence break
    for cur_genome in alt_core_pair_matrix.keys():
        # Find breaks for current genome
        breaks_cur_genome = [pair for pair in sequence_break_pairs if cur_genome in pair]
        if len(breaks_cur_genome):
            # Isolate core gene neighbour pairs for genome
            cur_genome_core_pairs = [pair.split('--')[0:2] for pair in master_info_total.keys() if cur_genome in pair]

            # Go through all alternative core pairs and examine the ones not present in the current genome.
            for alt_core_pair in alt_core_pair_matrix[cur_genome]:
                # Check that genome entry is not checked
                if alt_core_pair is not 'genome':
                    # Check that the alternative core neighbours are not present
                    if alt_core_pair_matrix[cur_genome][alt_core_pair] != 2:
                        # Split alternative variant into two elements in list
                        core_pair = alt_core_pair.split('--')

                        # Search if a combination of core genes near sequnce breaks can match the alternative pair,
                        # If then set evidence level to 1.
                        if len([gene for gene in breaks_cur_genome if gene[0] in core_pair[0]]) and \
                                len([gene for gene in breaks_cur_genome if gene[0] in core_pair[1]]) \
                                or \
                                len([gene for gene in breaks_cur_genome if gene[0] in core_pair[1]]) and \
                                len([gene for gene in breaks_cur_genome if gene[0] in core_pair[0]]):

                            # Check that the possible alternative pair is not because of many breaks in genome
                            alt_core_pair_possible = False

                            seq_breaks_gene_1 = determine_partners_neighbours(core_pair[0],
                                                                               core_genome_pairs,
                                                                               cur_genome_core_pairs)
                            seq_breaks_gene_2 = determine_partners_neighbours(core_pair[1],
                                                                              core_genome_pairs,
                                                                              cur_genome_core_pairs)
                            if seq_breaks_gene_1 > 0 and seq_breaks_gene_2 > 0:
                                alt_core_pair_possible = True

                            if alt_core_pair_possible:
                                alt_core_pair_matrix[cur_genome][alt_core_pair] = 1
                                core_variants_predicted_weak += 1
                            else:
                                alt_core_pair_matrix[cur_genome][alt_core_pair] = 2
                                core_variants_predicted_strong += 1

                            core_variants_predicted += 1


            # Increment the counter for variants searched
            genomes_searched += 1

    # TODO - set as verbose oberated:
    print(f'A total of {genomes_searched} genomes were searched for alternative'
          f' core-core neighbours separated by a sequence break')
    print(f'{core_variants_predicted} alternative core-core neighbours were predicted, when examining genomes with '
          f'sequnce breaks.')
    print(f'{core_variants_predicted_weak} of the predicted alternative core neighbours had weak evidence')
    print(f'{core_variants_predicted_strong} of the predicted alternative core neighbours had strong evidence')

    # Convert to list to be writen as output file
    matrix_list = [dict_content for dict_content in alt_core_pair_matrix.values()]

    return [matrix_list, header]
