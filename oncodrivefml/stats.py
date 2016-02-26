
import numpy as np
from scipy import stats


class ArithmeticMean(object):

    @staticmethod
    def calc(values):
        return np.mean(values)

    @staticmethod
    def calc_observed(values, observed_value):
        values = np.mean(np.array(list(values)), axis=1)
        obs = len(values[values >= observed_value])
        neg_obs = len(values[values <= observed_value])
        return obs, neg_obs

    @staticmethod
    def calc_observed_slow(values, observed_value):
        obs, neg_obs = 0, 0
        for values in values:
            background_test = np.mean(values)
            if background_test >= observed_value:
                obs += 1
            if background_test <= observed_value:
                neg_obs += 1
        return obs, neg_obs

    @staticmethod
    def weight(vectors, weights):
        return np.average(vectors, weights=weights, axis=0)


class GeometricMean(object):

    @staticmethod
    def calc(values):
        return stats.gmean(np.array(values) + 1.0) - 1.0

    @staticmethod
    def calc_observed(values, observed_value):
        obs, neg_obs = 0, 0
        for values in values:
            background_test = GeometricMean.calc(values)
            if background_test >= observed_value:
                obs += 1
            if background_test <= observed_value:
                neg_obs += 1
        return obs, neg_obs

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

