import numpy

def merge_dicts_counts(parent_dict, merge_dict):
    """ Function that can merge two dicts by keys and adding 1 to the value each time key is observed"""
    if isinstance(merge_dict, dict):
        keys = merge_dict.keys()
    elif isinstance(merge_dict, list):
        keys = merge_dict

    for key in keys:
        if key in parent_dict:
            parent_dict[key] += 1
        else:
            parent_dict[key] = 1

    return parent_dict


def merge_dicts_lists(parent_dict, merge_dict):
    """ Function to add two dictionaries by adding lists of matching keys """

    for key in merge_dict.keys():
        # Check if key is present, if then append the value to the key
        if key in parent_dict:
            if not isinstance(merge_dict[key], list):
                parent_dict[key] += [merge_dict[key]]
            else:
                parent_dict[key] += merge_dict[key]

        # If key is not present construct the key
        else:
            if not isinstance(merge_dict[key], list):
                parent_dict[key] = [merge_dict[key]]
            else:
                parent_dict[key] = merge_dict[key]

    return parent_dict


def merge_first_genes(start_gene_cluster, merged_start_gene_clusters, merged_second_gene_clusters, first_core_pair):
    """ Function that merge first genes and find second genes found in complete genomes,
    used to inform core-genome synteny direction and start gene """
    # Check if first gene is available
    if start_gene_cluster:
        # Record first gene
        merged_start_gene_clusters = merged_start_gene_clusters + [start_gene_cluster]

        # Record second gene in file
        second_gene = first_core_pair.split('--')
        second_gene.remove(start_gene_cluster)
        merged_second_gene_clusters = merged_second_gene_clusters + second_gene

    return merged_start_gene_clusters, merged_second_gene_clusters
