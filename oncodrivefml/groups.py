import numpy as np
import pandas as pd


def add_groups(results, groups_dict):
    """
    This function takes:
    - results: results from running Oncodrivefml dictionary with an entry for each successfully run gene.
    - groups_dict : a dictionary of groups of genes with name of the group as key and list of genes as value.

    Returns:
    - OncodriveFML results for each of the groups, when putting together the results of the genes belonging to each of the groups.
    """
    results_groups = dict()
    results_copy = results.copy()
    results_all = pd.DataFrame.from_dict(results, orient='index')
    results_all["scores"] = results_all["scores"].apply(lambda x: np.array(x))
    results_all["back_means"] = results_all["back_means"].apply(lambda x: np.array(x))

    for group in groups_dict:
        gene_list = groups_dict[group]
        results_group = results_all[results_all.index.isin(gene_list)]

        if results_group.shape[0] == 0:
            continue

        group_item = dict(results_group[['muts', 'muts_recurrence', 'snps', 'mnps', 'indels']].sum().items())
        group_item['samples_mut'] = results_group[['samples_mut']].max()[0]
        group_item['scores'] = np.concatenate(results_group["scores"])
        group_item['genes_in_group'] = ','.join(sorted(results_group.index.values))

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

        results_groups[group] = group_item

    
    # update the results dict to contain both results from individual genes and entire groups
    results_copy.update(results_groups)

    return results_groups, results_copy