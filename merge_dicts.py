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


# Seems depricated and uneccesary when .update is available
# def merge_master_dict(master_parent, master_merge):
#     for key in master_merge.keys():
#         master_parent[key] = master_merge[key]
#
#     return master_parent
