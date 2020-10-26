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
        if key in parent_dict:
            parent_dict[key] += merge_dict[key]
        else:
            parent_dict[key] = merge_dict[key]

    return parent_dict
