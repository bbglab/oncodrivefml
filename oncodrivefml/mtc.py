"""
Module containing functions related to
multiple test correction
"""

import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt


def multiple_test_correction(results, num_significant_samples=2):
    """
    Performs a multiple test correction on the analysis results

    Args:
        results (dict): dictionary with the results
        num_significant_samples (int): mininum samples that a gene must have in order to perform the correction

    Returns:
        :obj:`~pandas.DataFrame`. DataFrame with the q-values obtained from a multiple test correction

    """

    results_all = pd.DataFrame.from_dict(results, orient='index')

    # Filter minimum samples
    try:
        results_good = results_all[(results_all['samples_mut'] >= num_significant_samples) & (~results_all['pvalue'].isnull())].copy()
        results_masked = results_all[(results_all['samples_mut'] < num_significant_samples) | (results_all['pvalue'].isnull())].copy()
    except KeyError as e:
        raise e

    # Multiple test correction
    if len(results_good) > 1:
        results_good['qvalue'] = mlpt(results_good['pvalue'], alpha=0.05, method='fdr_bh')[1]
        results_good['qvalue_neg'] = mlpt(results_good['pvalue_neg'], alpha=0.05, method='fdr_bh')[1]
    else:
        results_good['qvalue'] = np.nan
        results_good['qvalue_neg'] = np.nan

    # Concat results
    results_concat = pd.concat([results_good, results_masked], sort=False)
    return results_concat
