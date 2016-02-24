import numpy as np
from scipy import stats


class ArithmeticMean(object):

    @staticmethod
    def calc(values):
        return np.mean(values)

    @staticmethod
    def weight(vectors, weights):
        return np.average(vectors, weights=weights, axis=0)


class GeometricMean(object):

    @staticmethod
    def calc(values):
        return stats.gmean(np.array(values) + 1.0) - 1.0

    @staticmethod
    def weight(vectors, weights):
        v_a = [np.array(i) + 1 for i in vectors]
        v_b = [np.power(i, wi) for i, wi in zip(v_a,weights)]
        total_weight = np.sum(weights)
        v_c = [(np.prod([j[i] for j in v_b])**(1/total_weight)) - 1.0 for i in range(len(vectors[0]))]
        return np.array(v_c)

STATISTIC_TESTS = {
    'amean': ArithmeticMean(),
    'gmean': GeometricMean()
}

