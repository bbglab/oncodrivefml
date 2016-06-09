
import numpy as np
from scipy import stats


class ArithmeticMean(object):

    @staticmethod
    def calc(values):
        return np.mean(values)

    @staticmethod
    def calc_observed(values, observed):
        """

        :param values:
        :param observed: [ scores of observed mutations ]
        :return:
        """
        observed_value = np.mean(observed)
        values = np.mean(values, axis=1)
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
        return stats.gmean(values + 1.0) - 1.0

    @staticmethod
    def calc_observed(values, observed):
        observed_value = stats.gmean(observed + 1.0) -1.0
        values = stats.gmean(values +1.0, axis=1) -1.0
        obs = len(values[values >= observed_value])
        neg_obs = len(values[values <= observed_value])
        return obs, neg_obs


STATISTIC_TESTS = {
    'amean': ArithmeticMean(),
    'maxmean': MaximumAndArithmeticMean(),
    'gmean': GeometricMean()
}

