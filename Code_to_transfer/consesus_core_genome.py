import networkx as nx
import numpy as np
from os.path import basename, join
import concurrent.futures
from itertools import repeat
from time import time


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

    return G


def clean_core_graph(core_gene_graph, output_path):
    """ Function to clean the core genome for any nodes that has more than two degrees. This also resolves any cliques.
    Additionally all genes with more than two dregrees are noted as being possible rearrangement spots."""

    # Get all nodes in graph
    graph_nodes = list(core_gene_graph.nodes())

    # Initialise gene and pas counters
    j = 0
    n_pas = 1

    # Initialise list to hold possible rearrangements
    possible_rearrangement_genes = []

    # Get the number of cliques and nodes that do not have 2 degrees
    n_cliques = len([clique for clique in nx.enumerate_all_cliques(core_gene_graph) if len(clique) >= 3])
    n_degrees = sum([1 for node in core_gene_graph.degree if node[1] != 2])

    # Print initial info
    print(f'Cleaning core genome graph\n') #TODO - MAKE MORE CLEAR

    print(f'Pas number: {n_pas}')
    print(f'Number of cliques present in graph: {n_cliques}')
    print(f'Number of nodes not having degree of 2 present in graph: {n_degrees}\n')

    # Start cleaning core gene graph until done
    while True:
        neighbours = 0

        # Search until a node with at least three degrees is found
        while neighbours < 3:
            # Check if the last node is reached and if there are no more cliques and if all nodes has 2 degrees,
            # if then return the graph - all is done
            if j == len(graph_nodes) and n_cliques == 0 and n_degrees == 0:
                print('Finished cleaning core genome graph. No more cliques and all nodes has two edges')
                nx.write_gml(core_gene_graph, join(output_path, 'cleaned_core_graph.gml'))
                return core_gene_graph, possible_rearrangement_genes

            # Check if last node is reached, then restart as not all nodes has been resolved for additional degrees
            elif j == len(graph_nodes):
                j = 0
                n_pas += 1

                nx.write_gml(core_gene_graph, join(output_path, f'prelim_core_graph_{n_pas}.gml'))
                print(f'Pas number: {n_pas}')
                print(f'Number of cliques present in graph: {n_cliques}')
                print(f'Number of nodes not having degree of 2 present in graph: {n_degrees}\n')

                if n_pas == 50:
                    cluster_2_degrees = [node[0] for node in core_gene_graph.degree if node[1] != 2]

                    nx.write_gml(core_gene_graph, join(output_path, 'prelim_core_graph.gml'))

                    raise NotImplementedError(f'The core graph cleaning process reached 50 rounds. Is it stuck? '
                                              f'Please check that the number of cliques and degrees not having 2 edges '
                                              f'has gone down during the couple of passes!\n'
                                              f'Core genes with multiple edges: {cluster_2_degrees}')

            # If the final node is not reached fetch the next node in line and examine its number of degrees
            else:
                edges_oi = core_gene_graph.edges(graph_nodes[j])

                neighbours = len(edges_oi)

                j += 1

        # Get the weight of each edge leading for a given node
        edge_weights = [core_gene_graph.get_edge_data(edge[0], edge[1])['weight'] for edge in edges_oi]

        # Find the least weighted edge
        min_edge_index = list(np.where(edge_weights == np.amin(edge_weights))[0])

        # If a single minimum weight can be found or
        # if the number of degrees is at least three larger than the number of degrees with the minimum weight,
        # then pressed remove the degrees
        if len(min_edge_index) == 1 or len(edges_oi) > len(min_edge_index) + 3:
            edges_oi = list(edges_oi)

            minimum_edge_pairs = [edges_oi[index] for index in min_edge_index]
            # print(f"minimum_edge_pairs {minimum_edge_pairs}")
            # Check if node at the end of edge has at least 3 or more edges, if then remove if not then keep edge
            # print(f'Node in question: {graph_nodes[j - 1]}')
            # neighbour_nodes = set([node for tup in minimum_edge_pairs for node in tup])

            # neighbour_nodes.remove(graph_nodes[j - 1])
            # print(f'Nodes round with low edge weight: {neighbour_nodes}')
            # minimum_edge_pairs = [node for node in neighbour_nodes if core_gene_graph.degree(node) > 2]
            minimum_edge_pairs = [pair for pair in minimum_edge_pairs if core_gene_graph.degree(pair[1]) > 2]
            # print(f'Number of degrees for neighbours: {[core_gene_graph.degree(node) for node in neighbour_nodes]}')
            # print(f'Nodes to which edges should be removed: {minimum_edge_pairs}')

            # [core_gene_graph.remove_edge(*list(edges_oi)[index]) for index in minimum_edge_pairs]
            [core_gene_graph.remove_edge(*pair) for pair in minimum_edge_pairs]

            # possible_rearrangement_genes += [edges_oi[index] for index in minimum_edge_pairs]
            possible_rearrangement_genes += [pair for pair in minimum_edge_pairs]

            # Recalculate degree and clique numbers
            n_cliques = len([clique for clique in nx.enumerate_all_cliques(core_gene_graph) if len(clique) >= 3])
            n_degrees = sum([1 for node in core_gene_graph.degree if node[1] != 2])


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
        # Connection to first core gene
        _, edge_weights_first = get_new_neighbours(second_gene, [], core_gene_graph)
        # Connection to second core gene
        third_gene, edge_weights_second = get_new_neighbours(second_gene, [start_gene], core_gene_graph)
        # start core coverage output
        core_path_coverage = [[start_gene, second_gene, edge_weights_first[0]],
                              [second_gene, third_gene[0], edge_weights_second[0]]]

    # TODO - Else, when only contigged genomes are pressent (Search for dnaA? first naively using regex in list then more sofisticated)

    # Initialise list to store genes of interest and narrow the search for rearrangements
    possible_rearrangement_genes = []

    # Search Graph for path
    # TODO - The logic of running until new_neoghbours is zero can be dangerous if only incomplete genomes are available
    #  as a break in the same place in all of them may lead to a truncated consensus. Implement a check to see
    #  if all core genes are included.

    # TODO - Fix so that all nodes are visited in graph to constuct consensus core genome.
    while len(new_neighbours) > 0:
        # Check if more than one neighbour is present if, then choose the one with largest edge weight,
        # Else walk along only path available
        if len(new_neighbours) > 1:
            max_weight = np.max(edge_weights)
            largest_weight_index = np.where([index == max_weight for index in edge_weights])[0]#[0]
            if len(largest_weight_index) == 1:
                best_neighbour = new_neighbours[largest_weight_index[0]]
            else:
                raise ValueError(f'Multiple core gene neighbours were equally likely to neighbour to: {visited_list[-1]},\n'
                                 f'The new neighbours variable was: {new_neighbours}\n'
                                 f'The edge weights were: {edge_weights}.')
                # TODO - Handle better!

        else:
            best_neighbour = new_neighbours[0]

        # Add next node in core genome synteny
        visited_list.append(best_neighbour)
        new_neighbours, edge_weights = get_new_neighbours(best_neighbour, visited_list, core_gene_graph)
        coverage_pairs = map(lambda best, new, cov:
                             [best, new, cov],
                             [best_neighbour]*len(new_neighbours), new_neighbours, edge_weights)
        # for i in coverage_pairs:
        #     print(i)
        # print(list(coverage_pairs))
        # core_path_coverage = core_path_coverage + [edge_weights]
        core_path_coverage = core_path_coverage + list(coverage_pairs)
        # TODO - Use core_path_coverage to produce a plot that show the how much coverage each edge has.
        # TDOD - possibly produce a list and an output file that contain only the coverage of consensus core gene synteny.

        # Check if more than one neighbour is available due to possible rearrangements
        if len(new_neighbours) > 1:
            # Map possible combinations of branch in core-core genome path
            possible_pairs = list(map(lambda e: [best_neighbour, e], new_neighbours))
            possible_rearrangement_genes = possible_rearrangement_genes + possible_pairs

    # Insert connections from last gene to first gene
    _, end_to_start_coverage = get_new_neighbours(visited_list[-1], [visited_list[-2]], core_gene_graph)
    core_path_coverage = core_path_coverage + [[visited_list[-1], start_gene, end_to_start_coverage[0]]]

    print(f"visited_list[-1] {visited_list[-1]}")
    print(f"start_gene{start_gene}")

    print(f"Number of nodes in core graph: {core_gene_graph.number_of_nodes()}")
    print(f'Number of nodes visited: {len(visited_list)}')

    return visited_list, possible_rearrangement_genes, core_path_coverage


def identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count):
    # Go through each segment flanked by a core gene with >2 edges and identify segments with no accessory between them

    # Create dict of subsegments of the larger segments
    # DO not use dict.fromkeys() with value being a list, as this uses the same list in all keys.
    sub_segment_dict = {key: [] for key in double_edge_segements}

    # Go through segments
    for segment in double_edge_segements:
        empty_segment_genes = []

        cur_segment = double_edge_segements[segment]
        # Check each region of the segment for core genes
        for i in range(0, len(cur_segment)-1):
            core_neighbours = sorted([cur_segment[i], cur_segment[i+1]])
            core_region = f'{core_neighbours[0]}--{core_neighbours[1]}'
            # Get accessory genes in region
            core_region_acc_genes = combined_acc_gene_count[core_region]

            # If number of core genes is zero then add the core pair to current segment and increment the counter for the length of current pair
            # Else add the segment to the directory of sub segments and reset counters.

            # If core region does not contain accessory genes, add to current segment. Else add the segment and start a new
            if core_region_acc_genes == 0:
                # If  first pair in segment add both, if not first only add the last gene
                if len(empty_segment_genes) == 0:
                    empty_segment_genes += [cur_segment[i], cur_segment[i+1]]
                else:
                    empty_segment_genes += [cur_segment[i+1]]

            else:
                # Check if first pair in subsegment and add first gene as being 'lonely'
                if len(empty_segment_genes) == 0:
                    empty_segment_genes += [cur_segment[i]]

                # Record the segment and reset the subsegment to contain no core genes
                sub_segment_dict[segment].append(empty_segment_genes)
                empty_segment_genes = []

            # Check of segment end has been reached and more than two genes are in the segment, if then add the segment
            if i == len(cur_segment) - 2 and len(empty_segment_genes) >= 2:
                sub_segment_dict[segment].append(empty_segment_genes)
                empty_segment_genes = []
            # Check if the second gene in pair is last in segment, and accessory genes are present between second to last and last core gene,
            # if then add the last gene as being 'lonely'
            elif i == len(cur_segment) - 2:
                empty_segment_genes += [cur_segment[i + 1]]
                sub_segment_dict[segment].append(empty_segment_genes)

    return sub_segment_dict


def determine_genome_segments(core_neighbour_pairs, combined_acc_gene_count, num_gffs):
    # Construct a graph from core gene neighbours
    core_graph = construct_core_graph(core_neighbour_pairs)

    # Find segments in the genome between core genes with multiple neighbors
    double_edge_segements, connect_dict, multi_edge_nodes = identify_segments(core_graph, num_gffs)

    if double_edge_segements is not None:
        # Find segments of core genes with no accessory in between
        no_acc_segments = identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)
    else:
        no_acc_segments = None

    return double_edge_segements, no_acc_segments


def determine_core_gene_consesus(core_neighbour_pairs, start_gene_cluster, second_gene_cluster, output_path):
    core_gene_graph = construct_core_graph(core_neighbour_pairs, output_path)

    # TODO - Clean up graph and find alternative core pairs
    core_gene_graph, possible_rearrangement_genes = clean_core_graph(core_gene_graph, output_path)

    consensus_core_genome, \
    _, \
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
    if '_corrected' in gff_names[0]:
        gff_names = [gff.split('_corrected')[0] for gff in gff_names]

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

    # Correct gff file names
    # if '_corrected' in gff_names[0]:
    #     gff_names = [name.split('corrected') for name in gff_names]

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
          f'{round(len(alt_core_pairs.keys())/len(gff_names) * 100, ndigits=1)}%\n') # TODO - Change to say number out of total instead of a percentage

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
    # print("predicting rearrangements")
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

    # TODO - set as verbose operated:
    print(f'A total of {genomes_searched} genomes were searched for alternative'
          f' core-core neighbours separated by a sequence break')
    print(f'{core_variants_predicted} alternative core-core neighbours were predicted, when examining genomes with '
          f'sequnce breaks.')
    print(f'{core_variants_predicted_weak} of the predicted alternative core neighbours had weak evidence')
    print(f'{core_variants_predicted_strong} of the predicted alternative core neighbours had strong evidence')

    # Convert to list to be writen as output file
    matrix_list = [dict_content for dict_content in alt_core_pair_matrix.values()]

    return [matrix_list, header]


def construct_directed_graph(core_graph):
    # Find a node that contains only two edges
    node_degrees = list(core_graph.degree)
    degree_counter = 0

    n_degrees = node_degrees[degree_counter][1]
    while True:
        degree_counter += 1

        if degree_counter == len(node_degrees):
            raise NotImplementedError("No nodes in the network has two degrees!")

        n_degrees = node_degrees[degree_counter][1]

        if n_degrees == 2:
            neighbours = list(core_graph.neighbors(node_degrees[degree_counter][0]))
            neighbour_1_ndegree, neighbour_2_ndegree = [core_graph.degree(neighbour) for neighbour in neighbours]

            if neighbour_1_ndegree == 2:
                split_index = 0

            elif neighbour_2_ndegree == 2:
                split_index = 1

            split_neighbour = neighbours[split_index]
            break

    # Extract the node with only two degrees
    split_node = node_degrees[degree_counter][0]

    print(f'split_node: {split_node}')
    print(f'split_neighbour: {split_neighbour}')
    # Remove the edge between the two nodes:
    core_graph.remove_edge(split_node, split_neighbour)

    # Go through and construct all edges in directed graph
    di_core_graph = nx.DiGraph()

    previous_nodes = [split_node]
    all_visited_nodes = []
    last_visited_nodes = []
    all_edges_added = []
    alternative_edges = []

    print(f'Number of of nodes in input graph: {len(core_graph.nodes)}')
    # REMOVE
    counter = 0

    while previous_nodes != [split_neighbour] and len(previous_nodes) != 0:
        if counter == 200:
            print(counter)
            break

        next_pairs = [(node, gene) for node in previous_nodes for gene in get_new_neighbours(node, last_visited_nodes, core_graph)[0]]

        # filter edges that are already added with direction
        next_pairs = [pair for pair in next_pairs if (pair[1], pair[0]) not in all_edges_added]

        dup_next_pair = [pair for pair in next_pairs if (pair[1], pair[0]) in next_pairs]

        # print(f'Number of duplicated pairs: {len(set(dup_next_pair))}')
        dup_next_pair = list(set(dup_next_pair))
        # Go though all duplicated pairs until they have been resolved
        while len(dup_next_pair) > 0:
            cur_pair = dup_next_pair.pop()

            # Constuct the alternative
            alt_dup_pair = (cur_pair[1], cur_pair[0])

            # Remove the pair and the alternative
            next_pairs = [pair for pair in next_pairs if pair != alt_dup_pair]
            next_pairs = [pair for pair in next_pairs if pair != cur_pair]

            alternative_edges.append((cur_pair[1], cur_pair[0]))

            dup_next_pair.remove((cur_pair[1], cur_pair[0]))

            next_pairs.append(cur_pair)

        [di_core_graph.add_edge(*pair, weight=core_graph.get_edge_data(*pair)['weight']) for pair in next_pairs]

        all_visited_nodes += [pair[1] for pair in next_pairs]
        all_edges_added += next_pairs
        last_visited_nodes = [pair[0] for pair in next_pairs]
        previous_nodes = [pair[1] for pair in next_pairs]

        counter += 1

    print(f'Number of node in di_core_graph: {len(di_core_graph.nodes)}')
    print(f'Number of edges {len(di_core_graph.edges)}')

    # print([node if node in enumerate(di_core_graph.degree) if node[1] > ])
    degrees = core_graph.degree
    for node in di_core_graph.degree:
        for ori_node in degrees:
            if ori_node[0] == node[0] and ori_node[1] < node[1]:
                print(f'node {node}')
                print(f'ori_node {ori_node}')

    return di_core_graph, alternative_edges
    # Remove one of the edges from that node

    # Remember what the nodes' names were

    # Make all edges from one node to the other directed with weight


def identify_segments(core_graph, num_gffs):
    """ Identify all segments (paths) of the graph the goes from one node with >2 degrees to the next,
    where all nodes in between contain only two degrees and thus the path has no ambiguity """

    # Identify all nodes that contain more than two degrees.
    multi_edge_nodes = [node for node, connections in core_graph.degree if connections > 2]

    # Check if any multi node edges are present, if not then return.
    if len(multi_edge_nodes) == 0:
        return None, None, None
        # raise NotImplementedError("A core gene graph with no nodes having more than two degrees was constructed.")

    # identify the multi connected nodes that are connected to one another, as these do not need to be searched for a simple path between them.
    connect_dict = {}

    # Identify neighbouring nodes with >2 degrees for all the nodes with >2 degrees themself
    for node in multi_edge_nodes:
        connect_dict[node] = [neighbor for neighbor in core_graph.neighbors(node) if neighbor in multi_edge_nodes]

    # Turn the weight into a 'distance' or number of times not found together.
    for edge in core_graph.edges(data=True):
        core_graph[edge[0]][edge[1]]['weight'] = num_gffs - core_graph[edge[0]][edge[1]]['weight']

    # find all simple paths between nodes with >2 degrees
    double_edge_segements = {}

    # Go through all source and taget nodes and see if a path can be found where all nodes between them have only two degrees
    for source_node in multi_edge_nodes:
        for target_node in multi_edge_nodes:
            if target_node != source_node and target_node not in connect_dict[source_node]:
                # Get path (segment) segment from source to target
                segment = nx.shortest_path(core_graph, source_node, target_node, weight='weight')

                # Get length of path
                segment_length = len(segment)

                # Get length of segment with multi nodes removed
                two_degree_segment_length = len([node for node in segment if node not in multi_edge_nodes])

                # Check if no node between the source and target has more than two edges, if not then record the segment/path
                if segment_length - 2 == two_degree_segment_length:
                    # Construct name for path
                    source_target_name = sorted([source_node, target_node])
                    source_target_name = f'{source_target_name[0]}--{source_target_name[1]}'

                    # Check that path has not been recorded in the opposite direction, if not then record it
                    if source_target_name not in double_edge_segements:
                        double_edge_segements[source_target_name] = segment
                    else:
                        if double_edge_segements[source_target_name] != segment[::-1]:
                            raise NotImplementedError("Path from one node to another was found, but did not match previously found path!")

    # Calculate the expected number of paths
    total_edges_from_multi_edge_nodes = sum([connections for _, connections in core_graph.degree if connections > 2])
    num_edges_between_multi_edge_nodes = sum([len(connect_dict[key]) for key in connect_dict])
    expected_segment_number = int((total_edges_from_multi_edge_nodes / 2) - (num_edges_between_multi_edge_nodes / 2))

    # Check if less than the number of expected paths has been found, if then try to identify missing paths
    if expected_segment_number != len(double_edge_segements):
        # Get number of expected edges for each node that has more than two edges
        expected_edge_num_dict = {node: connections for node, connections in core_graph.degree if connections > 2}

        # Get the number of edges directly between multi connected nodes
        identified_edge_num_dict = {node: len(connect_dict[node]) for node in connect_dict}

        # Get the number of paths connecting multi connected nodes via segments found in previous loop
        for connection in double_edge_segements:
            connection_nodes = connection.split('--')
            for node in connection_nodes:
                identified_edge_num_dict[node] += 1

        # Compare the identified number of connections expected to the identified, to find nodes that are missing connections
        nodes_missing_connections = []
        for node in expected_edge_num_dict:
            if identified_edge_num_dict[node] != expected_edge_num_dict[node]:
                nodes_missing_connections.append(node)

        # Go through nodes that are missing at least one path and try to identify missing paths
        for node in nodes_missing_connections:

            for current_target_node in nodes_missing_connections:
                # Check that the current target node is not a neighbouring node or the current node itself
                if current_target_node not in connect_dict[node] and current_target_node != node:

                    # Copy the graph to manipulate it
                    core_graph_copy = core_graph.copy()

                    # Extract the nodes for source and target
                    source_node = node
                    target_node = current_target_node

                    # Construct a pair name
                    suspected_pair = sorted([source_node, target_node])
                    suspected_pair = f'{suspected_pair[0]}--{suspected_pair[1]}'

                    # Check that the pair has not been found in a previous run
                    if suspected_pair not in double_edge_segements:
                        # Counter to stop loop
                        counter = 0
                        # Identifier to see if path has been found, to stop loop
                        path_identified = False
                        while not path_identified:
                            counter += 1

                            # Get all shortest path between source and target.
                            all_shortest_paths = nx.all_shortest_paths(core_graph_copy, source_node, target_node)

                            # Go through each path to see if is satisfies the criteria
                            try:
                                for index, path in enumerate(all_shortest_paths):
                                    # Get length of path
                                    segment_length = len(path)

                                    # Get length of segment with multi nodes removed
                                    two_degree_segment_length = len([node for node in path if node not in multi_edge_nodes])

                                    # Check that the path does not contain nodes with >2 degrees outside of source and target,
                                    # if then add path if not then find nodes that has >2 edges and remove an edge that leads to the to break path for next run through loop
                                    if segment_length - 2 == two_degree_segment_length and two_degree_segment_length != 0:
                                        # Add in path
                                        double_edge_segements[suspected_pair] = path
                                        path_identified = True
                                        pass

                                    else:
                                        # Check if path is length >2, if then find >2 degree nodes and remove an edge to them, if not just remove edge found between nodes.
                                        if len(path) > 2:
                                            multi_node_in_path = [[path[index], path[index+1]] for index, node in enumerate(path) if node in multi_edge_nodes and node != source_node and node != target_node]
                                            # print(f'list(set(multi_node_in_path)) {list(set(multi_node_in_path))}')
                                            for multi_node_pair in multi_node_in_path:
                                                # Try to remove edge found to multi node, if already removed move on.
                                                try:
                                                    core_graph_copy.remove_edge(*multi_node_pair)
                                                except nx.exception.NetworkXError:
                                                    continue
                                                    # raise ValueError("NetworkX was not able to remove an edge in the core graph network during segment identification!")

                                        else:
                                            core_graph_copy.remove_edge(*path)

                                    if counter == 1000:
                                        raise IndexError("Counter reached limit! in detecting a new path for pair.")
                            except nx.NetworkXNoPath:
                                # No simple paths could be found for the source and target thus the while loop is terminated.
                                path_identified = True

    return double_edge_segements, connect_dict, multi_edge_nodes



    # TODO - Add in genes that contain more than two edges and are only connected to other genes with more than two edges. These should be added as segments!

    # This idea scales exponentially with the number of degrees that has >2 degrees.

    # Another idea could be to move out from a node, until another node is identifies and then save the path.

def find_min_path(simple_path, recquired_len, i):
    print(f'process {i}')
    print(simple_path)
    print(recquired_len)
    if len(simple_path) == recquired_len:
        print("nice length")
        return True
    else:
        print("too short")
        return False


def connect_segments(double_edge_segements, connect_dict, multi_edge_nodes, core_graph):
    # Construct some represenatation of edge weights between end nodes in segments
    # Calculate the best path through them. All simple paths, that fulfill the number of all nodes and then search for largest value?

    # Build simplified network
    segment_graph = nx.Graph()

    for node in connect_dict:
        # Add in direct connections
        for conencted_node in connect_dict[node]:
            conenction_weight = core_graph.get_edge_data(node, conencted_node)['weight']
            conenction_distance = 532 - conenction_weight + 1

            segment_graph.add_weighted_edges_from([(node, conencted_node, conenction_distance)])

        # Add in segment connections
    for segment in double_edge_segements:
        source, target = segment.split('--')
        conenction_distance = 0

        segment_graph.add_weighted_edges_from([(source, target, conenction_distance)])

    # Remove all edges that connect a core gene, between two segments, to another core gene elsewhere in the graph.

    nx.write_gml(G=segment_graph, path='/Users/mjespersen/Documents/Davies_scripts/segment_graph.gml')

    all_simple_paths = nx.all_simple_paths(segment_graph, source, target)

    num_nodes = len(segment_graph.nodes)
    num_segments = sum([1 for *_, info in list(segment_graph.edges(data=True)) if info['weight'] == 0])

    # TODO - TRY TO SPEED UP THIS PROCEESS BY MULTIPROCESSING
    # TODO - insert check if more than one core is given, if then use multiprocessing, if not then run in for-loop.
    # with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
    #     results = [executor.submit(find_min_path, path, num_nodes, i) for i, path in enumerate(all_simple_paths)]
    #
    #     print(results)
    #     for re_value in concurrent.futures.as_completed(results):
    #         print(re_value.result())

    # for index, path in enumerate(all_simple_paths):
    #     if index == 100:
    #         exit()
    #
    #     print(find_min_path(path, num_nodes))

    # counter = 1
    min_path = []
    min_path_weight = 999999999
    best_path_weights = []
    alternative_bests = []
    # TODO - speed this up by multi proccessing. Divide the list of all simple paths into smaller chunks and get the smallest value as the return from them.
    start_time = time()
    for index, path in enumerate(all_simple_paths):
        if len(path) == num_nodes:

            sum_cur_path_weight = 0
            cur_path_weights = []
            for pair in zip(path[1:], path[:-1]):
                path_weight = segment_graph.get_edge_data(*pair)['weight']
                sum_cur_path_weight += path_weight
                cur_path_weights.append(path_weight)

            # Check if the found path has a accumulated wight small than the previous best edge and if all segments has been visited
            # if min_path_weight > sum_cur_path_weight and cur_path_weights.count(0) == num_segments:
            if sum_cur_path_weight < min_path_weight: # cur_path_weights.count(0) == num_segments:
                print('new best path')
                min_path = path
                min_path_weight = sum_cur_path_weight
                best_path_weights = cur_path_weights
                num_segments = cur_path_weights.count(0)
                print(f"num_segments {num_segments}")

                # Rest alternative bests
                alternative_bests = []
                print(best_path_weights)

            elif min_path_weight > sum_cur_path_weight:
                alternative_bests.append(path)

        if index % 10000000 == 0:
            print(f'Round: {index} reached')

    print(f"num_segments {num_segments}")
    print(f'time: {time() - start_time}')

if __name__ == '__main__':
    G = nx.read_gml('/Users/mjespersen/Downloads/core_graph.gml')

    double_edge_segements, connect_dict, multi_edge_nodes = identify_segments(G)

    connect_segments(double_edge_segements, connect_dict, multi_edge_nodes, G)

# TODO - Find a way to give the segments in an output.
#   * Find all segments that contain no core genes with >2 degrees
#   * Find all segments that contain no accesory genes between core genes.