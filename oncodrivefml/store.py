"""
This module contains the methods used to store the results.
"""

import gzip
import os


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def add_symbol(df):
    ensemble_file = os.path.join(__location__, "ensembl_genes_75.txt.gz")
    gene_conversion = {line.split("\t")[0]: line.strip().split("\t")[-1]
                       for line in gzip.open(ensemble_file, 'rt').readlines()}
    gene_symbols = df.GENE_ID.apply(lambda e: gene_conversion.get(e, e))
    df.SYMBOL.fillna(gene_symbols, inplace=True)
    return df


def store_tsv(results, result_file):
    """
    Saves the results in a tsv file sorted by pvalue

    Args:
        results (:obj:`~pandas.DataFrame`): results of the analysis
        result_file: file where to store the results

    """
    results.index.names = ['GENE_ID']
    results.sort_values(by='pvalue', inplace=True)
    fields = ['muts', 'muts_recurrence', 'samples_mut', 'pvalue', 'qvalue', 'pvalue_neg', 'qvalue_neg', 'snps', 'mnps', 'indels', 'symbol']
    df = results[fields].copy()
    df.reset_index(inplace=True)
    df.rename(columns={'muts': 'MUTS', 'muts_recurrence': 'MUTS_RECURRENCE', 'samples_mut': 'SAMPLES',
                       'pvalue': 'P_VALUE', 'qvalue': 'Q_VALUE','pvalue_neg': 'P_VALUE_NEG',
                       'qvalue_neg': 'Q_VALUE_NEG', 'snps': 'SNP', 'mnps':'MNP', 'indels': 'INDELS',
                       'symbol': 'SYMBOL'}, inplace=True)
    df = add_symbol(df)

    df.to_csv(result_file, sep="\t", header=True, index=False, compression="gzip")
