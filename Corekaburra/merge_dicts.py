def merge_dicts_counts(parent_dict, merge_object):
    """
    Function that can merge two dicts by keys and adding 1 to the value each time key is observed
    :param parent_dict: Dict to which the second with should be merged into
    :param merge_object: Dict or List to be merged into the first.

    :return: Resulting dict following merge
    """

    keys = merge_object.keys() if isinstance(merge_object, dict) else merge_object

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

    for key, value in merge_dict.items():
        value_as_list = value if isinstance(value, list) else [value]
        if key in parent_dict:
            parent_dict[key] += value_as_list
        else:
            parent_dict[key] = value_as_list

    return parent_dict
