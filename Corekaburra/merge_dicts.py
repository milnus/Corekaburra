def merge_dicts_counts(parent_dict, merge_object):
    """
    Function that can merge two dicts by keys and adding 1 to the value each time key is observed
    :param parent_dict: Dict to which the second with should be merged into
    :param merge_object: Dict or List to be merged into the first.

    :return: Resulting dict following merge
    """

    if isinstance(merge_object, dict):
        keys = merge_object.keys()
    elif isinstance(merge_object, list):
        keys = merge_object

    for key in keys:
        if key in parent_dict:
            parent_dict[key] += 1
        else:
            parent_dict[key] = 1

    return parent_dict


def merge_dicts_lists(parent_dict, merge_dict):
    """
    Function to add two dictionaries by adding lists of matching keys
    :param parent_dict: The Dict to which the second dict should be merged with
    :param merge_dict: Dict to be merge with the parent

    :return: Dict having the two inputs merged
    """

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
