import numpy as np


def calculate_n_create_summaries(master_info): # TODO - Add in columns that gives that difference in length and acc genes from smallest to biggest.
    # Create dict to hold the information for each core pair with a key being the pair name and the value a 2 dimentional numpy array
    summary_dict = {}

    # Group the values of length of accessory gene content for the core pairs into the summary dict
    for pair in master_info.keys():
        core_pair = pair.rsplit('--', 1)[0]

        # Try to add the values for a given core pair to the summary dict,
        # if not possible then add a two column array to they key and add the values
        try:
            summary_dict[core_pair].append(master_info[pair][3:5])
        except KeyError:
            summary_dict[core_pair] = [master_info[pair][3:5]]
            # summary_dict[core_pair] = np.append(summary_dict[core_pair], master_info[pair][3:5])

    # Take the mean over length or accesory genes for all core-pairs
    for core_pair in summary_dict:
        # cast the list into a numpy array, format is to have lengths and accessory gene count as a column each
        summary_dict[core_pair] = np.array(summary_dict[core_pair])
        summary_dict[core_pair] = summary_dict[core_pair].reshape((-1, 2))

        # Calculate summary statistics
        pair_occurrence = summary_dict[core_pair].shape[0]
        min_values = np.min(summary_dict[core_pair], axis=0)
        max_values = np.max(summary_dict[core_pair], axis=0)
        mean_values = np.mean(summary_dict[core_pair], axis=0)
        median_values = np.median(summary_dict[core_pair], axis=0)

        # Add the summary to the dict, round mean and median to one decimal
        summary_dict[core_pair] = [core_pair, pair_occurrence,
                                   min_values[0], max_values[0],
                                   round(mean_values[0], 1), round(median_values[0], 1),
                                   min_values[1], max_values[1],
                                   round(mean_values[1], 1), round(median_values[1], 1)
                                   ]

    return summary_dict

if __name__ == '__main__':
    master_info = {
        'pan_cluster_1--pan_cluster_2--genome_1': ['genome_1', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                   ['Acc_1', 'Acc_2'], ['low_1']],
        'pan_cluster_1--pan_cluster_2--genome_2': ['genome_2', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                   ['Acc_1', 'Acc_2'], ['low_1']],
        'pan_cluster_1--pan_cluster_2--genome_3': ['genome_3', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                   ['Acc_1', 'Acc_2'], ['low_1']],
        'pan_cluster_2--pan_cluster_3--genome_1': ['genome_1', 'pan_cluster_2', 'pan_cluster_3', 100, 2,
                                                   ['Acc_1', 'Acc_2'], []],
        'pan_cluster_2--pan_cluster_3--genome_2': ['genome_2', 'pan_cluster_2', 'pan_cluster_3', 150, 1,
                                                   ['Acc_1', ], []],
        'pan_cluster_2--pan_cluster_3--genome_3': ['genome_3', 'pan_cluster_2', 'pan_cluster_3', 200, 0,
                                                   [], []]
                   }

    calculate_n_create_summaries(master_info)