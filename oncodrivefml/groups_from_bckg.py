"""
Module containing functions related to
multiple test correction
"""
import re
import sys
import json
import typer
import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt

app = typer.Typer()

def multiple_test_correction(results_all, num_significant_samples=2):
    """
    Performs a multiple test correction on the analysis results

    Args:
        results (dict): dictionary with the results
        num_significant_samples (int): mininum samples that a gene must have in order to perform the correction

    Returns:
        :obj:`~pandas.DataFrame`. DataFrame with the q-values obtained from a multiple test correction

    """

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
    results_concat['scores'] = results_concat['scores'].apply(np.mean)
    results_concat.rename(columns={'gene_id':'GENE_ID', 'muts': 'MUTS', 'muts_recurrence': 'MUTS_RECURRENCE',
                        'scores' : 'AVG_SCORE_OBS',
                        'samples_mut': 'SAMPLES',
                        'pvalue': 'P_VALUE', 'qvalue': 'Q_VALUE','pvalue_neg': 'P_VALUE_NEG',
                        'qvalue_neg': 'Q_VALUE_NEG', 'snps': 'SNP', 'mnps':'MNP', 'indels': 'INDELS',
                        'population_mean' : 'POPULATION_MEAN', 'population_std' : 'POPULATION_STD', 'std_distribution_of_means': 'STD_OF_MEANS',
                        'z-score' : 'Z-SCORE',
                        'symbol': 'SYMBOL',
                        'genes_in_group': 'GENES_IN_GROUP'}, inplace=True)
    results_concat = results_concat[['GENE_ID', 'MUTS', 'MUTS_RECURRENCE', 'SAMPLES', 'AVG_SCORE_OBS', 
                                        'P_VALUE', 'Q_VALUE', 'P_VALUE_NEG', 'Q_VALUE_NEG',
                                        'SNP', 'MNP', 'INDELS',
                                        'POPULATION_MEAN', 'POPULATION_STD', 'STD_OF_MEANS',
                                        'Z-SCORE', 'SYMBOL', 'GENES_IN_GROUP']]
    return results_concat

def custom_replace(x):
        if x:
            if type(x) == str:
                return np.array(eval(x)) if ',' in x[:100] else np.array(eval(re.sub(r'(?<=\d|\.)\s+', ', ', x)))
            else:
                return np.array(x)
        return np.array()

def add_groups(results_all_file, groups_dict):
    """
    This function takes:
    - results: results from running Oncodrivefml dictionary with an entry for each successfully run gene.
    - groups_dict : a dictionary of groups of genes with name of the group as key and list of genes as value.

    Returns:
    - OncodriveFML results for each of the groups, when putting together the results of the genes belonging to each of the groups.
    """
    results_groups = dict()
    results_all = pd.read_csv(results_all_file, header = 0, sep = "\t")
    results_all.rename(columns={'GENE_ID':'gene_id', 'MUTS': 'muts', 'MUTS_RECURRENCE': 'muts_recurrence',
                                'SAMPLES': 'samples_mut', 'SNP': 'snps', 'MNP': 'mnps', 'INDELS': 'indels',
                                'POPULATION_MEAN': 'population_mean', 'POPULATION_STD': 'population_std',
                                'STD_OF_MEANS': 'std_distribution_of_means', 'SCORE_OBS': 'scores',
                                'BACK_MEANS': 'back_means', 'SYMBOL': 'symbol'}, inplace=True)

    # results_all["scores"] = results_all["scores"].apply(lambda x: np.array(eval(x)))
    # results_all["back_means"] = results_all["back_means"].apply(lambda x: np.array(eval(x)))
    results_all["scores"] = results_all["scores"].apply(lambda x: custom_replace(x))
    results_all["back_means"] = results_all["back_means"].apply(lambda x: custom_replace(x))
    floatts = results_all[results_all["back_means"].apply(lambda x: type(x) == float)]
    if floatts.shape[0] > 0:
        for index, row in floatts.iterrows():
            print(row[['gene_id', 'muts', 'muts_recurrence', 'snps', 'mnps', 'indels', 'scores']], row["back_means"])
        results_all = results_all[results_all["back_means"].apply(lambda x: type(x) != float)]

    diff_len = results_all[results_all["back_means"].apply(lambda x: len(x) != 100000)]
    if diff_len.shape[0] > 0:
        for index, row in diff_len.iterrows():
            print(row[['gene_id', 'muts', 'muts_recurrence', 'snps', 'mnps', 'indels', 'scores']], len(row["back_means"]))
        results_all = results_all[results_all["back_means"].apply(lambda x: len(x) == 100000)]


    for group in groups_dict:
        gene_list = groups_dict[group]
        results_group = results_all[results_all["gene_id"].isin(gene_list)].reset_index(drop = True)
        
        if results_group.shape[0] == 0:
            continue

        group_item = dict(results_group[['muts', 'muts_recurrence', 'snps', 'mnps', 'indels']].sum().items())
        group_item['samples_mut'] = results_group[['samples_mut']].max()[0]
        group_item['scores'] = np.concatenate(results_group["scores"])
        # group_item['scores'] = np.array( [x for x in y for y in results_group["scores"] ] )
        group_item['genes_in_group'] = ','.join(sorted(results_group["gene_id"].values))

        mean_score_obs = (results_group["scores"].apply(sum).sum()) / (group_item["muts"])
        mean_scores_bckg = ((results_group["muts"] * results_group["back_means"]).sum()) / (group_item["muts"])

        weighted_mean_from_pop_bckg = ((results_group["muts"] * results_group["population_mean"]).sum()) / (group_item["muts"])

        weighted_std_of_means = np.average(results_group["population_std"], weights=results_group["muts"])
        sem = weighted_std_of_means / np.sqrt( group_item["muts"] )

        group_item['z-score'] = ( mean_score_obs - weighted_mean_from_pop_bckg ) / sem
        
        group_item['population_mean'] = weighted_mean_from_pop_bckg
        group_item['population_std'] = weighted_std_of_means
        group_item['std_distribution_of_means'] = sem

        group_item['pvalue'] = np.sum(mean_scores_bckg > mean_score_obs) / len(mean_scores_bckg)
        group_item['pvalue_neg'] = np.sum(mean_scores_bckg < mean_score_obs) / len(mean_scores_bckg)

        group_item['symbol'] = group
        group_item['gene_id'] = group

        results_groups[group] = group_item

    results_groups_df = pd.DataFrame.from_dict(results_groups, orient='index').reset_index(drop = False)

    return results_groups_df.reset_index(drop = False)


if __name__ == '__main__':
    bckg_file = sys.argv[1]
    groups_json_file = sys.argv[2]
    store_results = sys.argv[3]
    number_of_decimals = 6

    with open(groups_json_file, 'r') as file_groups:
        groups_json = json.load(file_groups)
    group_results_df = add_groups(bckg_file, groups_json)

    if group_results_df.shape[0] > 0:
        group_results_corrected_df = multiple_test_correction(group_results_df)
        group_results_corrected_df.sort_values(by='P_VALUE', inplace=True)
        
        float_columns = group_results_corrected_df.select_dtypes(include=['float64']).columns
        group_results_corrected_df[float_columns] = group_results_corrected_df[float_columns].round(number_of_decimals)

        group_results_corrected_df.to_csv(store_results, sep="\t", header=True, index=False, compression="gzip")
