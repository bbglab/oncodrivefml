import argparse
import json
import logging
import datetime
import signaturesampling as sampling
import os
import pandas as pd
import numpy as np
import gzip

from statsmodels.sandbox.stats.multicomp import multipletests as mlpt

STEP1_VERSION = "v03"
STEP2_VERSION = "v01"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Constants
NC_PATH = "/projects_bg/bg/shared/datasets/"
EXPRESSION_DICTIONARY = {'astrocytoma': 'BRAIN_rnaseq_medians_EXP_STATED_ENS',
                         'breast': 'BREAST_rnaseq_medians_EXP_STATED_ENS',
                         'cll': 'BLOOD_rnaseq_medians_EXP_STATED_ENS',
                         'liver': 'PANCAN19_rnaseq_medians_EXP_STATED_ENS',
                         'lung': 'LUNG_rnaseq_medians_EXP_STATED_ENS',
                         'lymphoma': 'BLOOD_rnaseq_medians_EXP_STATED_ENS',
                         'medulloblastoma': 'BRAIN_rnaseq_medians_EXP_STATED_ENS',
                         'stad': 'STOMACH_rnaseq_medians_EXP_STATED_ENS',
                         'pan': 'PANCAN19_rnaseq_medians_EXP_STATED_ENS',
                         'blca': 'BLADDER_rnaseq_medians_EXP_STATED_ENS',
                         'brca': 'BREAST_rnaseq_medians_EXP_STATED_ENS',
                         'crc': 'COLON_AND_RECTUM_rnaseq_medians_EXP_STATED_ENS',
                         'gbm': 'BRAIN_rnaseq_medians_EXP_STATED_ENS',
                         'hnsc': 'OROPHARYNX_rnaseq_medians_EXP_STATED_ENS',
                         'kich': 'kidney_rnaseq_medians_EXP_STATED_ENS',
                         'kirc': 'kidney_rnaseq_medians_EXP_STATED_ENS',
                         'lgg': 'BRAIN_rnaseq_medians_EXP_STATED_ENS',
                         'luad': 'LUNG_rnaseq_medians_EXP_STATED_ENS',
                         'lusc': 'LUNG_rnaseq_medians_EXP_STATED_ENS',
                         'prad': 'PROSTATE_rnaseq_medians_EXP_STATED_ENS',
                         'skcm': 'PANCAN19_rnaseq_medians_EXP_STATED_ENS',
                         'thca': 'THYROID_rnaseq_medians_EXP_STATED_ENS',
                         'ucec': 'corpus_uteri_rnaseq_medians_EXP_STATED_ENS'}

SIGNATURE_MAP_PROJECT = { 'pilot_simulation_1kb': 'pilot',
                          'pilot_simulation_10kb': 'pilot',
                          'pilot_simulation_50kb': 'pilot',
                          'tcga_simulation_1kb': 'tcga',
                          'tcga_simulation_10kb': 'tcga',
                          'tcga_simulation_50kb': 'tcga',
                          'targeted': None}

SIGNATURE_MAP_TISSUE = { 'intogen':
    { 
    'ALL_Leuven': ('intogen', 'all'),
     'ALL_Memphis': ('intogen', 'all'),
     'ALL_NewSamples': ('intogen', 'all'),
     'tcga_laml': ('intogen', 'aml'),
     'tcga_blca': ('intogen', 'blca'),
     'icgc_snp_hopkins_Breast': ('intogen', 'brca'),
     'icgc_snp_sanger_Breast': ('intogen', 'brca'),
     'literature_breast_BC': ('intogen', 'brca'),
     'literature_breast_WU': ('intogen', 'brca'),
     'literature_breast_carcinoma_BROAD': ('intogen', 'brca'),
     'tcga_breast': ('intogen', 'brca'),
     'tcga_cesc': ('intogen', 'cesc'),
     'Broad_CLL': ('intogen', 'cll'),
     'icgc_snp_es_CLL': ('intogen', 'cll'),
     'literature_CLL_DF': ('intogen', 'cll'),
     'Melanoma_BROAD': ('intogen', 'cm'),
     'tcga_skcm': ('intogen', 'cm'),
     'icgc_snp_hopkins_Colon': ('intogen', 'coread'),
     'tcga_coadread': ('intogen', 'coread'),
     'Bcell_lymphoma': ('intogen', 'dlbc'),
     'DLBCL_broad': ('intogen', 'dlbc'),
     'Esophageal_DFCI': ('intogen', 'esca'),
     'icgc_snp_hopkins_GBM': ('intogen', 'gbm'),
     'tcga_glioblastoma': ('intogen', 'gbm'),
     'Liver_ICGC': ('intogen', 'hc'),
     'icgc_snp_fr_Liver': ('intogen', 'hc'),
     'literature_HNSCC_BROAD': ('intogen', 'hnsc'),
     'tcga_hnsc': ('intogen', 'hnsc'),
     'tcga_lgg': ('intogen', 'lgg'),
     'icgc_snp_tsp_Lung': ('intogen', 'luad'),
     'tcga_luad': ('intogen', 'luad'),
     'tcga_lusc': ('intogen', 'lusc'),
     'Medulloblastoma_NewSamples': ('intogen', 'mb'),
     'icgc_snp_de_Pedbrain': ('intogen', 'mb'),
     'Myeloma_NewSamples': ('intogen', 'mm'),
     'mm_broad': ('intogen', 'mm'),
     'Neuroblastoma_BROAD': ('intogen', 'nb'),
     'literature_NSCLC_WISCONSIN': ('intogen', 'nsclc'),
     'tcga_ovary': ('intogen', 'ov'),
     'Astrocytoma_NEWSAMPLES': ('intogen', 'pa'),
     'Astrocytoma_StJude': ('intogen', 'pa'),
     'icgc_snp_hopkins_Panc': ('intogen', 'paad'),
     'icgc_snp_oicr_Panc': ('intogen', 'paad'),
     'icgc_snp_qcmg_Panc': ('intogen', 'paad'),
     'tcga_paad': ('intogen', 'paad'),
     'PRAD_broad': ('intogen', 'prad'),
     'tcga_prad': ('intogen', 'prad'),
     'tcga_kidney': ('intogen', 'rccc'),
     'tcga_kirp': ('intogen', 'rpcc'),
     'literature_SCLC_JHU': ('intogen', 'sclc'),
     'literature_SCLC_UCOLOGNE': ('intogen', 'sclc'),
     'literature_gastric_PFIZER': ('intogen', 'stad'),
     'tcga_stad': ('intogen', 'stad'),
     'tcga_thca': ('intogen', 'thca'),
     'tcga_uteri': ('intogen', 'ucec')
    }
}


def _signature(project, tissue):

    if project in SIGNATURE_MAP_PROJECT:
        project = SIGNATURE_MAP_PROJECT[project]

    if project in SIGNATURE_MAP_TISSUE:
        if tissue in SIGNATURE_MAP_TISSUE[project]:
            project, tissue = SIGNATURE_MAP_TISSUE[project][tissue]

    if project is None or tissue is None:
        return 'none'

    return "{}_{}".format(project, tissue)


def _multiple_test_correction(feature, results, tissue):
    results_all = pd.DataFrame.from_dict(results, orient='index')
    # Filter expression of utr
    if feature == 'utr_5and3':
        expression_file = os.path.join(NC_PATH, "somatic_mutations/expression_filters_backup/",
                                       EXPRESSION_DICTIONARY[tissue])
        expression = [line.strip() for line in open(expression_file, 'r')]
        results_not_expressed = results_all[~results_all.index.isin(expression)].copy()

        # Filter the dataframe for expression
        results_all = results_all[results_all.index.isin(expression)]
    else:
        results_not_expressed = pd.DataFrame(columns=results_all.columns)

    # Filter minimum samples
    num_significant_samples = 5 if tissue == 'pan' else 2
    results_good = results_all[results_all['samples_mut'] >= num_significant_samples].copy()
    results_masked = results_all[results_all['samples_mut'] < num_significant_samples].copy()
    # Multiple test correction
    if len(results_good) > 1:
        results_good['qvalue'] = mlpt(results_good['pvalue'], alpha=0.05, method='fdr_bh')[1]
    else:
        results_good['qvalue'] = np.nan

    # Concat results
    results_concat = pd.concat([results_good, results_masked, results_not_expressed])
    return results_concat


def _run(project, tissue, score, feature):

    # Initialize
    if score not in sampling.SCORES:
        raise RuntimeError("Unknown score '{}'".format(score))

    # Skip if done
    output_folder = os.path.join(sampling.OUTPUT_FOLDER, project, tissue, feature)
    pvalues_file = os.path.join(output_folder, 'pvalues_{}_{}.tsv.gz'.format(score, STEP2_VERSION))
    if os.path.exists(pvalues_file):
        logger.info("{} - {} - {} - {} [skip]".format(project, tissue, score, feature))
        return
    sampling._silent_mkdir(output_folder)

    # Loading results
    logger.info("Loding results")
    results_file = os.path.join(sampling.OUTPUT_FOLDER, project, tissue, feature, "results_{}_{}.json.gz".format(score, STEP1_VERSION))
    if not os.path.exists(results_file):
        raise RuntimeError("Results file '{}' not found".format(results_file))
    with gzip.open(results_file, 'rt') as fd:
        results = json.load(fd)

    # Calculate empirical p-values
    logger.info("Calculate empirical p-values")
    if len(results) == 0:
        logger.error("There are no results")
        return

    num_randomizations = 10000
    to_run = [(element, means) for element, means in results.items()]
    while len(to_run) > 0:

        # Preprocess samplings
        next_to_run = []
        next_num_randomizations = num_randomizations * 2
        to_prepare = []
        logger.info("Round {} permutations with {} elements [start]".format(num_randomizations, len(to_run)))
        for element, m in to_run:

            for m_tissue in m['muts_by_tissue']['subs'].keys():

                m_signature = _signature(project, m_tissue)
                m_scores = m['muts_by_tissue']['subs'][m_tissue]
                m_count = len(m_scores)

                if num_randomizations == 10000:
                    sampling_folder = os.path.join(sampling.SAMPLING_CACHE_FOLDER, score, feature, element[len(element) - 2:], element.upper())
                    sampling_file = os.path.join(sampling_folder, 'sampling-{}-{}.bin'.format(m_signature, m_count))
                    if os.path.exists(sampling_file):
                        # If file exists skip sampling
                        continue
                    else:
                        background_file = os.path.join(sampling_folder, 'background.tsv')
                        if os.path.exists(background_file) and os.stat(background_file).st_size == 0:
                            # If the background file is empty skip sampling
                            continue

                to_prepare.append((element, m_count, m_signature))

        if len(to_prepare) > 0:
            logger.info("Send {} samplings to prepare".format(len(to_prepare)))
            sampling.sampling_prepare(to_prepare, feature=feature, score=score, verbose=True, max_jobs=450, sampling_size=num_randomizations)

        logger.info("Start sampling")
        for i, (e, m) in enumerate(to_run):

            if i % 100 == 0:
                logger.info("[{} of {}]".format(i, len(to_run)))

            values_mean = None
            values_mean_count = 0
            all_scores = []
            for m_tissue in m['muts_by_tissue']['subs'].keys():

                m_signature = _signature(project, m_tissue)
                m_scores = m['muts_by_tissue']['subs'][m_tissue]
                m_count = len(m_scores)

                values = sampling.sampling(score=score,
                       signature=m_signature,
                       feature=feature,
                       element=e,
                       num_samples=m_count,
                       sampling_size=num_randomizations)

                if values is None:
                    logger.warning("There are no scores at {}-{}-{}-{}-{}-{}".format(score, m_signature, feature, e, m_count, num_randomizations))
                    continue

                values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count], axis=0) if values_mean is not None else values
                values_mean_count += m_count

                all_scores += m_scores

            #TODO Compute indels

            obs = len(values_mean[values_mean >= np.mean(all_scores)]) if len(all_scores) > 0 else float(num_randomizations)

            # Check if we need more resolution at this element
            if obs <= 3 and next_num_randomizations <= 180000:
                next_to_run.append((e, m))
                logger.warning("We need more permutations at {}-{}-{}-{} (obs: {})".format(score, feature, e, num_randomizations, obs))

            m['pvalue'] = max(1, obs) / float(num_randomizations)

        # Next iteration with the elements that need more permutations
        to_run = next_to_run
        num_randomizations = next_num_randomizations

    # Run multiple test correction
    logger.info("Multiple test correction")
    results_concat = _multiple_test_correction(feature, results, tissue)

    # Sort and store results
    results_concat.sort('pvalue', 0, inplace=True)
    fields = ['symbol', 'muts', 'muts_recurrence', 'samples_mut', 'samples_max_recurrence', 'subs', 'subs_score', 'indels', 'indels_score', 'pvalue', 'qvalue', 'all_mean', 'max_mean', 'scores', 'positions']
    with gzip.open(pvalues_file, 'wt') as fd:
        results_concat[fields].to_csv(fd, sep="\t", header=True, index=True)

    logger.info("{} - {} - {} - {} [done]".format(project, tissue, score, feature))

if __name__ == "__main__":

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')

    # Parse the arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("project")
    parser.add_argument("tissue")
    parser.add_argument("score")
    parser.add_argument("feature")

    args = parser.parse_args()

    _run(args.project, args.tissue, args.score, args.feature)
