import numpy as np

from oncodrivefml.stats import STATISTIC_TESTS


def flatten_partitions(results):

    for name, result in results.items():
        for partition in result['partitions']:
            yield (name, partition, result)


def partitions_list(total_size, chunk_size):
    """
    Create a list of values less or equal to chunk_size that sum total_size

    :param total_size: Total size
    :param chunk_size: Chunk size
    :return: list of integers
    """
    partitions = [chunk_size for _ in range(total_size // chunk_size)]

    res = total_size % chunk_size
    if res != 0:
        partitions += [res]

    return partitions


def compute_sampling(value):
    name, samples, result = value

    muts_count = result['muts_count']
    observed = result['observed']
    items_to_simulate = result['simulation_items']
    items_to_simulate_prob = result['simulation_probs']

    r = _compute_sampling(samples, muts_count, observed, items_to_simulate, items_to_simulate_prob)

    # TODO update result

    return


def _compute_sampling(samples, muts, observed, items, probs):
    # TODO add test
    raise NotImplementedError
