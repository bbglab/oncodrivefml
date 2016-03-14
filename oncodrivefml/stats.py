
import numpy as np
from scipy import stats


class ArithmeticMean(object):

    @staticmethod
    def calc(values):
        return np.mean(values)

    @staticmethod
    def calc_observed(values, observed):
        observed_value = np.mean(observed)
        values = np.mean(np.array(list(values)), axis=1)
        obs = len(values[values >= observed_value])
        neg_obs = len(values[values <= observed_value])
        return obs, neg_obs


class MaximumAndArithmeticMean(ArithmeticMean):

    @staticmethod
    def calc(values):
        return np.max(values)


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


STATISTIC_TESTS = {
    'amean': ArithmeticMean(),
    'maxmean': MaximumAndArithmeticMean(),
    'gmean': GeometricMean()
}

