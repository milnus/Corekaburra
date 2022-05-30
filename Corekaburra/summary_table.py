import numpy as np
# pylint: disable=E4702

try:
    from Corekaburra.consesus_core_genome import count_gene_co_occurrence
except ModuleNotFoundError:
    from consesus_core_genome import count_gene_co_occurrence


def calculate_n_create_summaries(master_info, core_gene_dict):
    """
    Function used to calculate summary statistics for core gene pairs
    :param master_info: Dict over master info from core pairs identified from gff files
    :param core_gene_dict: Dict mapping genes to genomes and pan-genome cluster (genes)
    :return: Dict containing the summary statistics for each core pair.
    """
    # initialize the summary dict to be returned
    summary_dict = {}

    # Group the values of length of accessory gene content for the core pairs into the summary dict
    for pair in master_info:
        core_pair = pair.rsplit('--', 1)[0]

        # Try to add the values for a given core pair to the summary dict,
        # if not possible then add a two column array to they key and add the values
        try:
            summary_dict[core_pair].append(master_info[pair][3:5])
        except KeyError:
            summary_dict[core_pair] = [master_info[pair][3:5]]

    occurrence_dict = {key: {} for key in summary_dict}
    # dict.fromkeys(summary_dict.keys(), {})
    for core_pair in summary_dict:
        gene_list = core_pair.split('--')
        # Calculate the occurrence and co-occurrence of core genes

        # Calculate gene occurrence and co-occurrence
        occurrence_dict[core_pair]['co_occurrence'], individual_occurrences = count_gene_co_occurrence(core_gene_dict, gene_list)

        for gene in gene_list:
            occurrence_dict[core_pair][gene] = individual_occurrences[gene]

    # Take the mean over length of accessory genes for all core-pairs
    for core_pair in summary_dict:
        gene_list = core_pair.split('--')
        # cast the list into a numpy array, format it to have lengths and accessory gene count as a column each
        summary_dict[core_pair] = np.array(summary_dict[core_pair])
        summary_dict[core_pair] = summary_dict[core_pair].reshape((-1, 2))

        # Calculate summary statistics
        pair_occurrence = summary_dict[core_pair].shape[0]
        min_values = np.min(summary_dict[core_pair], axis=0)
        max_values = np.max(summary_dict[core_pair], axis=0)
        mean_values = np.mean(summary_dict[core_pair], axis=0)
        median_values = np.median(summary_dict[core_pair], axis=0)

        # Add the summary to the dict, round mean and median to one decimal
        summary_dict[core_pair] = [core_pair.replace('--', '-'),
                                   pair_occurrence,
                                   occurrence_dict[core_pair][gene_list[0]],
                                   occurrence_dict[core_pair][gene_list[1]],
                                   occurrence_dict[core_pair]['co_occurrence'],
                                   min_values[0], max_values[0],
                                   round(mean_values[0], 1), round(median_values[0], 1),
                                   min_values[1], max_values[1],
                                   round(mean_values[1], 1), round(median_values[1], 1)
                                   ]

    return summary_dict


if __name__ == '__main__':
    pass
