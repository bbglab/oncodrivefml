import numpy as np

from oncodrivefml.stats import STATISTIC_TESTS
from oncodrivefml.walker_cython import walker_sampling


def flatten_partitions(results):

    for name, result in results.items():
        for partition in result['partitions']:
            yield (name, partition, result, np.random.randint(0, 2**32-1))


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
    name, samples, result, seed = value

    scores = result['simulation_scores']
    muts_count = result['muts_count']
    probs = result['simulation_probs']
    observed = result['observed']
    statistic_name = result['statistic_name']

    np.random.seed(seed)

    if statistic_name == "amean":
        obs, neg_obs = compute_sampling_cython(samples, muts_count, np.mean(observed), np.array(scores), np.array(probs))
    else:
        obs, neg_obs = compute_sampling_python(samples, muts_count, observed, scores, probs, statistic_name)

    return name, obs, neg_obs


def compute_sampling_python(samples, muts, observed, scores, probs, statistic_name):
    statistic_test = STATISTIC_TESTS.get(statistic_name)
    background = np.random.choice(scores, size=(samples, muts), p=probs, replace=True)
    return statistic_test.calc_observed(background, np.array(observed))


def compute_sampling_cython(samples, muts, obs_val, scores, probs):

    # Walker alias initialization
    size = len(scores)
    probs = probs * size
    inx = -np.ones(size, dtype=int)
    short = np.where(probs < 1)[0].tolist()
    long = np.where(probs > 1)[0].tolist()
    while short and long:
        j = short.pop()
        k = long[-1]

        inx[j] = k
        probs[k] -= (1 - probs[j])
        if probs[k] < 1:
            short.append(k)
            long.pop()

    # Check maximum to avoid long overflow
    if samples < 2000000000:
        seed = np.random.randint(0, 2**31-1)
        return walker_sampling(samples, muts, obs_val, scores, probs, inx, seed)
    else:
        obs, neg_obs = 0, 0
        for p in partitions_list(samples, 2000000000):
            seed = np.random.randint(0, 2**31-1)
            o, no = walker_sampling(p, muts, obs_val, scores, probs, inx, seed)
            obs += o
            neg_obs += no

    return obs, neg_obs


if __name__ == "__main__":
    import time

    size = 1000
    samples = 1000000
    muts = 50

    np.random.seed(0)

    scores = np.random.rand(size)
    probs = np.random.rand(size)
    probs = probs / sum(probs)
    obs_val = np.mean(np.random.rand(muts))

    s = time.time()
    obs, neg_obs = compute_sampling_python(samples, muts, obs_val, scores, probs, 'amean')
    print("Python time: {} Obs:{} Neg_obs:{}".format(time.time() - s, obs, neg_obs))

    s = time.time()
    obs, neg_obs = compute_sampling_cython(samples, muts, obs_val, scores, probs)
    print("Cython time: {} Obs:{} Neg_obs:{}".format(time.time() - s, obs, neg_obs))
