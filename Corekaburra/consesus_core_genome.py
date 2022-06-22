import networkx as nx
import concurrent.futures

try:
    from Corekaburra.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error
EXIT_SEGMENT_IDENTIFICATION_ERROR = 4

# pylint: disable=E1123, E1121


def construct_core_graph(core_neighbour_pairs):
    """
    Function to construct a graph from the core pairs and number of times each is observed.
    :param core_neighbour_pairs: Dict of core pairs and the number of times each is observed
    :return: A graph with nodes being core genes, edges being a connection between them with a weight of the number of times they are connected
    """
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


def count_gene_co_occurrence(core_gene_dict, two_gene_segment):
    """
    Function to find the number of genomes in which two genes co-occur across the input genomes.
    :param core_gene_dict: Dictionary over core genes mapped from genome, to locus_tag, to pan-genome cluster
    :param two_gene_segment: List of two genes forming a segment
    :return: Int - number of co-occurrences for the two genes in the input two_gene_segment
    """
    co_occurrence = 0
    gene_occurrence = dict.fromkeys(two_gene_segment, 0)

    # Get pan-genome clusters for all genomes in a list of lists
    core_gene_presences = [list(core_genes.values()) for core_genes in core_gene_dict.values()]

    # Go through all genomes and check if genes co-occur
    for core_gene_set in core_gene_presences:
        # count the co-occurrences
        if set(two_gene_segment).issubset(core_gene_set):
            co_occurrence += 1

        # Count the individual occurrences
        if two_gene_segment[0] in core_gene_set:
            gene_occurrence[two_gene_segment[0]] += 1
        if two_gene_segment[1] in core_gene_set:
            gene_occurrence[two_gene_segment[1]] += 1

    return co_occurrence, gene_occurrence


def identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count):
    """
    Function that takes segments between multi connected core genes, and a Dict of number of accessory genes in core-core regions.
    Divides segments into smaller subsegments, in which no accessory genes can be found between core pairs.
    :param double_edge_segements: Dict of segments of core genes identified. Keys are genes at edges of segments. Value is a List of genes in the segment from one side to the other.
    :param combined_acc_gene_count: Dict of the number of accessory genes (value) identified between a set core genes (Key)
    :return: Dict of subsegments. Same keys as for the segment dict, but keys are a list of lists. Each sub-list is a subsegment.
    """
    # TODO - ATM this does not occur if every core gene only has two connections. Should is still occur to let a complete static genome synteny be divided into no-accessory segments?
    # Create dict of subsegments of the larger segments
    sub_segment_dict = {key: [] for key in double_edge_segements}

    # Go through segments to identify subsegments
    for segment in double_edge_segements:
        empty_segment_genes = []

        cur_segment = double_edge_segements[segment]
        # Check each region of the segment for core genes
        for i in range(0, len(cur_segment)-1):
            core_neighbours = sorted([cur_segment[i], cur_segment[i+1]])
            core_region = f'{core_neighbours[0]}--{core_neighbours[1]}'
            # Get accessory genes in region
            core_region_acc_genes = combined_acc_gene_count[core_region]

            # If core region does not contain accessory genes, add to current segment.
            # Else add the segment and start a new
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


def search_for_path(core_graph_copy, source_node, target_node, multi_edge_nodes):
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
                # if then add path,
                # else then find nodes that has >2 edges and remove an edge that leads to the node, to break the path for next run through loop
                if segment_length - 2 == two_degree_segment_length and two_degree_segment_length != 0:
                    # double_edge_segements[suspected_pair] = path
                    return path

                else:
                    # Check if path is length >2,
                    # if then find >2 degree nodes and remove an edge to them,
                    # else just remove edge found between nodes.
                    if len(path) > 2:
                        multi_node_in_path = [[path[index], path[index + 1]] for index, node in enumerate(path) if
                                              node in multi_edge_nodes and node != source_node and node != target_node]
                        for multi_node_pair in multi_node_in_path:
                            # Try to remove edge found to multi node, if already removed move on.
                            try:
                                core_graph_copy.remove_edge(*multi_node_pair)
                            except nx.exception.NetworkXError:
                                continue
                    else:
                        core_graph_copy.remove_edge(*path)

                if counter == 1000:
                    raise IndexError(
                        f"Counter reached limit in detecting a new path for pair, with name {target_node = } and {source_node = }")
        except nx.NetworkXNoPath:
            # No simple paths could be found for the source and target thus the while loop is terminated.
            return


def identify_segments(core_graph, num_gffs, core_gene_dict, logger):
    """
    Function to identify stretches of core genes between core genes neighbouring multiple different genes
    :param core_graph: Graph over core genes with weights being the number of connections between the genes
    :param num_gffs: Number of gffs inputted
    :param core_gene_dict: Dict with keys being genomes, each genome is a dict with keys being genes and values the mapped pan-genome gene cluster.
    :param logger: Logger for the program

    :return: Dict over stretches of core genes found in the core gene graph.
    """

    # Identify all nodes that contain more than two degrees and only one degree.
    multi_edge_nodes = [node for node, connections in core_graph.degree if connections > 2]
    single_edge_nodes = [node for node, connections in core_graph.degree if connections == 1]

    # Check if any node have multiple edges, if not then return.
    if len(multi_edge_nodes+single_edge_nodes) == 0:
        return None

    # Dict to hold connections between >2 edge nodes
    connect_dict = {}

    # for all nodes with >2 degrees themself, identify neighbouring nodes with >2 degrees
    for node in multi_edge_nodes+single_edge_nodes:
        connect_dict[node] = [neighbor for neighbor in core_graph.neighbors(node)
                              if neighbor in multi_edge_nodes or neighbor in single_edge_nodes]

    # Turn the weight into a 'distance' or number of times not found together.
    for edge in core_graph.edges(data=True):
        core_graph[edge[0]][edge[1]]['weight'] = num_gffs - core_graph[edge[0]][edge[1]]['weight']

    # find all simple paths between nodes with >2 degrees
    double_edge_segements = {}
    multi_edge_connect_adjust = []

    # Go through all source and taget nodes,
    # see if a path can be found where all nodes between them have only two degrees
    for source_node in multi_edge_nodes+single_edge_nodes:
        for target_node in multi_edge_nodes+single_edge_nodes:
            if target_node != source_node:
                # Get path (segment) from source to target
                segment = nx.shortest_path(core_graph, source_node, target_node, weight='weight', method='dijkstra')

                # Get length of path
                segment_length = len(segment)

                # Get length of segment with multi nodes removed
                two_degree_segment_length = len([node for node in segment if node not in multi_edge_nodes+single_edge_nodes])

                # Check if no node between the source and target has more than two edges,
                # if then move to record the segment/path
                if segment_length - 2 == two_degree_segment_length:
                    # Check if two gene segment occur in every possible genome, if not then skip
                    if segment_length == 2:
                        gene_co_occurrences, _ = count_gene_co_occurrence(core_gene_dict, segment)
                        if num_gffs - core_graph[segment[0]][segment[1]]['weight'] < gene_co_occurrences:
                            continue
                        else:
                            # Check if segment has been added in opposite direction, if not they add it to be further examined
                            if all([x != segment[::-1] for x in multi_edge_connect_adjust]): multi_edge_connect_adjust.append(segment)

                    # Construct name for path
                    source_target_name = sorted([source_node, target_node])
                    source_target_name = f'{source_target_name[0]}--{source_target_name[1]}'

                    # Check that path has not been recorded in the opposite direction, if not then record it
                    if source_target_name not in double_edge_segements:
                        double_edge_segements[source_target_name] = segment
                    else:
                        if double_edge_segements[source_target_name] != segment[::-1]:
                            exit_with_error(EXIT_SEGMENT_IDENTIFICATION_ERROR,
                                            f"Path from one node to another ({source_target_name}) was found, but did not match previously found path!", logger)


    # Calculate the expected number of paths
    total_edges_from_non_two_edge_core_genes = sum([connections for _, connections in core_graph.degree if connections > 2 or connections < 2])
    num_edges_between_non_two_edge_core_genes = sum([len(connect_dict[key]) for key in connect_dict])
    expected_segment_number = int((total_edges_from_non_two_edge_core_genes / 2) - (num_edges_between_non_two_edge_core_genes / 2)) + len(multi_edge_connect_adjust)

    # Check if less than the number of expected paths has been found,
    # if then try to identify missing paths
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

        # Adjust for the number of segments identified between multi connected nodes
        for segment in multi_edge_connect_adjust:
            for node in segment:
                identified_edge_num_dict[node] -= 1

        # Compare the number of connections expected to the number identified, to find nodes that are miss connections
        nodes_missing_connections = []
        for node in expected_edge_num_dict:
            if identified_edge_num_dict[node] != expected_edge_num_dict[node]:
                nodes_missing_connections.append(node)

        # Go through nodes that are missing at least one path and try to identify missing paths
        for source_node in nodes_missing_connections:
            for target_node in nodes_missing_connections:

                # Check that the source and target are not the same node
                if target_node != source_node:
                    # Copy the graph to manipulate it
                    core_graph_copy = core_graph.copy()

                    # Construct a pair name
                    suspected_pair = sorted([source_node, target_node])
                    suspected_pair = f'{suspected_pair[0]}--{suspected_pair[1]}'

                    # Check that the current target node is not a neighbouring node
                    if target_node not in connect_dict[source_node]:
                        # Search for path
                        return_path = search_for_path(core_graph_copy, source_node, target_node, multi_edge_nodes)

                    else:
                        # Remove the link between the two core genes that are neighbours
                        core_graph_copy.remove_edge(*suspected_pair.split('--'))
                        # Search for path
                        return_path = search_for_path(core_graph_copy, source_node, target_node, multi_edge_nodes)

                    # Check if proper path is returned and insert it
                    if return_path is not None:
                        double_edge_segements[suspected_pair] = return_path

    return double_edge_segements


def determine_genome_segments(core_neighbour_pairs, combined_acc_gene_count, num_gffs, core_gene_dict, max_cpus, logger):

    """
    Function to be called from main that collects the functions for determining core segments in pan-genome

    :param core_neighbour_pairs: Dict of the number of times core pairs have been detected
    :param combined_acc_gene_count: Number of accessory and low-frequency genes detected between core gene pairs
    :param num_gffs: Number of inputted gff files
    :param core_gene_dict: A dictionary of core genes across genomes and their identifier
    :param max_cpus: Int for the maximum number of cpus allowed to be used during graph component search
    :param logger: Program logger

    :return double_edge_segements:
    :return no_acc_segments:
    """

    logger.debug(f"--------------Searching for segments in pan genome--------------")

    # Construct a graph from core gene neighbours
    core_graph = construct_core_graph(core_neighbour_pairs)
    num_core_graph_components = nx.number_connected_components(core_graph)

    logger.debug(f'Identified: {num_core_graph_components} components in core genome graph')

    double_edge_segements = {}
    # Identify all segments in components of core graph
    #for component in nx.connected_components(core_graph):

    logger.debug(f'Searching components of core gene graph')
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_cpus) as executor:
        return_object = [executor.submit(identify_segments,
                                         core_graph.subgraph(component).copy(), num_gffs,
                                         core_gene_dict, logger)
                         for component in nx.connected_components(core_graph)]
     #       identify_segments(core_graph.subgraph(component).copy(), num_gffs, core_gene_dict, num_core_graph_components, logger)
        for output in concurrent.futures.as_completed(return_object):
            return_segments = output.result()
            if return_segments is not None:
                double_edge_segements = double_edge_segements | return_segments


    # if double_edge_segements is not None:
    if double_edge_segements:
        logger.debug(f'A total of {len(double_edge_segements)} core genes were identified to have multiple neighbours.')
        logger.debug(f'Genes with multiple neighbours: {double_edge_segements}')

        logger.debug('Search for Segments with no accessory genes starts now')

        # Find segments of core genes with no accessory in between
        no_acc_segments = identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        logger.debug('Segments with no accessory genes is done')
    else:
        logger.debug(f'No segments can be identified in given pan-genome\n')
        no_acc_segments = None

    return double_edge_segements, no_acc_segments, core_graph


if __name__ == '__main__':
    print('Nothing was computed. This is not a main program. Run __main__.py')
    pass
