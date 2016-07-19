"""
This modules contains different statistical methods used to compare
the observed and the simulated scores
"""


import numpy as np
from scipy import stats


class ArithmeticMean(object):

    @staticmethod
    def calc(values):
        """
        Computes the arithmetic mean

        Args:
            values (:obj:`list`, :obj:`numpy.array`): array of values

        Returns:
            float: mean value

        """
        return np.mean(values)

    @staticmethod
    def calc_observed(values, observed):
        """
        Measure how many times the mean of the values is higher than the mean of the observed values

        Args:
            values (:obj:`numpy.array`):  m x n matrix with scores (m: number of randomizations; n: number of mutations)
            observed (list, :obj:`numpy.array`): n size vector with the observed scores (n: number of mutations)

        Returns:
            tuple: the number of times that the mean value of a randomization is greater or equal than the mean observed value
             (as :obj:`int`) and the number of times that the mean value of a randomization is equal or lower than the mean
             observed value (as :obj:`int`).

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
    """
    The geometric mean used is not the standard.

    .. math::

        (\prod \limits_{i=1}^n (x_i+1))^{1/n}-1 &= \sqrt[n]{(x_1+1)(x_2+1) \cdots (x_n+1)} -1
    """

    @staticmethod
    def calc(values):
        """
        Computes the geometric mean of a set of values.

        Args:
            values (:obj:`list`, :obj:`numpy.array`): set of values

        Returns:
            (float): geometric mean
            (array): geometric mean by columns (if the input is a matrix)

        """
        return stats.gmean(values + 1.0) - 1.0

    @staticmethod
    def calc_observed(values, observed):
        """
         Measure how many times the geometric mean of the values is higher than the geometric mean of the observed values

        Args:
            values (:obj:`numpy.array`):  m x n matrix with scores (m: number of randomizations; n: number of mutations)
            observed (list, :obj:`numpy.array`): n size vector with the observed scores (n: number of mutations)

        Returns:
            tuple: the number of times that the mean value of a randomization is greater or equal than the mean observed value
             (as :obj:`int`) and the number of times that the mean value of a randomization is equal or lower than the mean
             observed value (as :obj:`int`).

        """
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

