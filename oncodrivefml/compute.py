import numpy as np
from scipy import stats



def gmean(a):
    return stats.gmean(np.array(a) + 1.0) - 1.0


def gmean_weighted(vectors, weights):
    v_a = [np.array(i) + 1 for i in vectors]
    v_b = [np.power(i, wi) for i, wi in zip(v_a,weights)]
    total_weight = np.sum(weights)
    v_c = [(np.prod([j[i] for j in v_b])**(1/total_weight)) - 1.0 for i in range(len(vectors[0]))]
    return np.array(v_c)



def random_scores(num_samples, sampling_size, background, signature, statistic_name):

    result = None
    if num_samples > 0:

        if len(background) == 0:
            return None

        if signature is None:
            # Subs sampling without signature
            p_normalized = None

        else:
            # Subs sampling with signature
            p_normalized = np.array(signature[:len(background)]) / sum(signature)

        if num_samples == 1:
            result = np.array(
                np.random.choice(background, size=sampling_size, p=p_normalized, replace=True),
                dtype='float32'
            )
        else:

            # Select mean
            if statistic_name == 'gmean':
                statistic_test = gmean
            elif statistic_name == 'rmean':
                statistic_test = rmean
            elif statistic_name == 'max':
                statistic_test = np.max
            else:
                statistic_test = np.mean

            result = np.array(
                [statistic_test(np.random.choice(background, size=num_samples, p=p_normalized, replace=True)) for a in range(sampling_size)],
                dtype='float32'
            )

    return result[:sampling_size]


# def by_segment():
#
#     scores_by_segment = defaultdict(list)
#     signature_by_segment = defaultdict(defaultdict_list)
#
#     # Signature
#     if alt is None or alt == '.':
#         scores_by_segment[region['segment']].append(value)
#         scores_by_segment[region['segment']].append(value)
#         scores_by_segment[region['segment']].append(value)
#     else:
#         scores_by_segment[region['segment']].append(value)
