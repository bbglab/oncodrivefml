import json
import logging
import datetime
from oncodrivefm2.configuration import Configuration, SCORES
import os
import pandas as pd
import numpy as np
import gzip

from statsmodels.sandbox.stats.multicomp import multipletests as mlpt
from oncodrivefm2.signaturesampling import SignatureSampling, silent_mkdir

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class OncodriveFM2Test(object):
    def __init__(self, config: Configuration, sampling: SignatureSampling):
        self._conf = config
        self._output_file = None
        self._sampling = sampling


    def get_output_file(self):
        return self._output_file

    def _multiple_test_correction(self, results):
        feature = self._conf.feature
        tissue = self._conf.tissue

        results_all = pd.DataFrame.from_dict(results, orient='index')
        # Filter expression of utr
        if feature == 'utr_5and3':
            expression = [line.strip() for line in open(self._conf.expression_file, 'r')]
            results_not_expressed = results_all[~results_all.index.isin(expression)].copy()

            # Filter the dataframe for expression
            results_all = results_all[results_all.index.isin(expression)]
        else:
            results_not_expressed = pd.DataFrame(columns=results_all.columns)

        # Filter minimum samples
        num_significant_samples = 5 if tissue == 'pan' else 2
        results_good = results_all[results_all['muts'] >= num_significant_samples].copy()
        results_masked = results_all[results_all['muts'] < num_significant_samples].copy()
        # Multiple test correction
        if len(results_good) > 1:
            results_good['max_qvalue'] = mlpt(results_good['max_pvalue'], alpha=0.05, method='fdr_bh')[1]
            results_good['all_qvalue'] = mlpt(results_good['all_pvalue'], alpha=0.05, method='fdr_bh')[1]
        else:
            results_good['max_qvalue'] = np.nan
            results_good['all_qvalue'] = np.nan

        # Concat results
        results_concat = pd.concat([results_good, results_masked, results_not_expressed])
        return results_concat


    def run(self):
        feature = self._conf.feature
        tissue = self._conf.tissue
        project = self._conf.project
        score = self._conf.score
        cache_dir = self._conf.cache_dir

        # Initialize
        signature = self._conf.signature
        if score not in SCORES:
            raise RuntimeError("Unknown score '{}'".format(score))

        # Loading results
        logger.info("Loading score sampling output")
        output_folder = self._conf.output_dir
        sampling_results_file = self._conf.input
        if not os.path.exists(sampling_results_file):
            raise RuntimeError("Results file '{}' not found".format(sampling_results_file))
        with gzip.open(sampling_results_file, 'rt') as fd:
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
            logger.info("Starting {} permutations with {} elements [start]".format(num_randomizations, len(to_run)))
            for element, means in to_run:

                if num_randomizations == 10000:
                    sampling_folder = os.path.join(output_folder, score, feature, element[len(element) - 2:], element.upper())
                    sampling_file = os.path.join(sampling_folder, 'sampling-{}-{}.bin'.format(signature, means['muts']))
                    if os.path.exists(sampling_file):
                        # If file exists skip sampling
                        continue
                    else:
                        background_file = os.path.join(sampling_folder, 'background.tsv')
                        if os.path.exists(background_file) and os.stat(background_file).st_size == 0:
                            # If the background file is empty skip sampling
                            continue

                to_prepare.append((element, means['muts']))

            if len(to_prepare) > 0:
                logger.info("Send samplings to prepare")
                self._sampling.sampling_prepare(to_prepare, signature=signature, verbose=True, max_jobs=450, sampling_size=num_randomizations)

            logger.info("Start sampling")
            for i, (e, v) in enumerate(to_run):

                if i % 100 == 0:
                    logger.info("[{} of {}]".format(i, len(to_run)))

                values = self._sampling._sampling(
                       signature=signature,
                       element=e,
                       num_samples=v['muts'],
                       sampling_size=num_randomizations)

                if values is None:
                    logger.warning("There are no scores at {}-{}-{}-{}-{}-{}".format(score, signature, feature, e, v['muts'], num_randomizations))
                    continue

                obs = len(values[values >= v['all_mean']])
                obs_by_sample = len(values[values >= v['max_mean']])

                # Check if we need more resolution at this element
                if obs <= 3 and next_num_randomizations <= 180000:
                    next_to_run.append((e, v))
                    logger.warning("We need more permutations at {}-{}-{}-{}-{}-{} (obs: {}, obs_by_sample: {})".format(score, signature, feature, e, v['muts'], num_randomizations, obs, obs_by_sample))

                v['all_pvalue'] = max(1, obs) / float(num_randomizations)
                v['max_pvalue'] = max(1, obs_by_sample) / float(num_randomizations)

            # Next iteration with the elements that need more permutations
            to_run = next_to_run
            num_randomizations = next_num_randomizations

        # Run multiple test correction
        logger.info("Multiple test correction")
        results_concat = self._multiple_test_correction(results)

        # Sort and store results
        results_concat.sort('all_pvalue', 0, inplace=True)
        silent_mkdir(output_folder)
        fields = ['symbol', 'muts', 'muts_recurrence', 'samples_mut', 'samples_max_recurrence', 'subs', 'subs_score', 'indels', 'indels_score', 'all_pvalue', 'all_qvalue', 'all_mean', 'max_pvalue', 'max_qvalue', 'max_mean', 'scores', 'positions']
        output_file = os.path.join(output_folder, self._conf.PVALUES_FILE_NAME.format(score, datetime.date.today()))
        logger.info("Writing to output file {} ".format(output_file))
        with gzip.open(output_file, 'wt') as fd:
            results_concat[fields].to_csv(fd, sep="\t", header=True, index=True)
        self._output_file = output_file

        logger.info("{} - {} - {} - {} [done]".format(project, tissue, score, feature))
