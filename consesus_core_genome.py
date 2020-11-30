import networkx as nx
import numpy as np


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
            largest_weight_index = np.where(np.max(edge_weights))[0][0]
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


def identify_rearrangements(consensus_core_genome, possible_rearrangement_genes, master_info_dict):
    # Pair all neighbouring genes in the consensus core genome
    core_genome_pairs = []
    for i in range(len(consensus_core_genome)):
        sorted_neighbours = sorted([consensus_core_genome[i], consensus_core_genome[i-1]])
        core_genome_pairs.append(f'{sorted_neighbours[0]}--{sorted_neighbours[1]}')

    core_genome_pairs = set(core_genome_pairs)

    # Search each key that contain a gene with possible rearrangements.
    # If two core gene pairs is not present then predict rearrangement and size. - store genome and rearrangements
    # If only one is present predict possible rearrangements, but unsure due to only one core pair identified
    # If three or more then it is complex and no further prediction can be made. - Possible contamination?

    # initialize dict of alternative core neighbours
    alt_core_pairs = {}

    # Identify the alternative core neighbours for each genome:
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
          f'{round(len(alt_core_pairs.keys())/332 * 100, ndigits=1)}%')

    # Record the number of times each alternative core pair occur
    alt_core_pair_count = {}
    for key in alt_core_pairs.keys():
        for pair in alt_core_pairs[key]:
            if pair in alt_core_pair_count.keys():
                alt_core_pair_count[pair] += 1
            else:
                alt_core_pair_count[pair] = 1

    return alt_core_pairs, alt_core_pair_count


def characterise_rearrangements(alt_core_pairs):
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

        # If more than two are present then more than one solution is possible and no prediction can be made.
        if len(alt_core_pairs[genome]) > 2:
            rearrangement_predictions[0].append(genome)
            rearrangement_predictions[1].append('Possible prediction')

    return rearrangement_predictions # TODO - Use this output
