import argparse
import gzip
import json
import logging
import numpy as np
import pandas as pd
import os

from multiprocessing.pool import Pool
from oncodrivefm2.utils import _load_variants, _file_name, _load_regions, _silent_mkdir, _intersect_mutations, _compute_score_means, _multiple_test_correction, \
    _create_background_signature, _sampling


class OncodriveFM2(object):

    def __init__(self, variants_file, regions_file, signature_file, score_file, output_folder, project_name=None, cores=os.cpu_count()):

        # Configuration
        self.cores = cores
        self.variants_file = variants_file
        self.regions_file = regions_file
        self.signature_file = signature_file
        self.score_file = score_file
        self.output_folder = output_folder
        self.project_name = project_name if project_name is not None else _file_name(variants_file)
        self.signature_name = _file_name(signature_file)

        # Output files
        self.variants_dict_file = os.path.join(output_folder, 'variants.json.gz')
        self.scores_dict_file = os.path.join(output_folder, 'scores.json.gz')
        self.results_file = os.path.join(output_folder, self.project_name + '-oncodrivefm2.tsv')

        # Some initializations
        _silent_mkdir(output_folder)

    def run(self):

        variants_dict = self._load_variants_dict()
        scores_dict = self._load_scores_dict(variants_dict)
        self._run_sampling(scores_dict)

    def _load_variants_dict(self):

        # Create variants dictionary
        if not os.path.exists(self.variants_dict_file):

            # Initialize pool
            pool = Pool(self.cores)

            # Loading mutations
            logging.info("Loading mutations")
            muts = _load_variants(self.variants_file, signature=self.signature_name)

            # Load features
            logging.info("Loading regions")
            regions = _load_regions(self.regions_file)

            # Aggregate mutations
            logging.info("Grouping mutations by element")
            mutations = {}
            elements_all = [(e, segments) for e, segments in regions]
            elements_chunks = np.array_split(elements_all, self.cores)
            compute_arguments = ((chunk, muts, num) for num, chunk in enumerate(elements_chunks, start=1))
            c = 0
            for done in pool.starmap(_intersect_mutations, compute_arguments):
                c += 1
                logging.info("Chunk {} of {} [done]".format(c, self.cores))
                for element, variants in done:
                    mutations[element] = variants

            # Store results
            logging.info("Storing variants dictionary")
            with gzip.open(self.variants_dict_file, 'wt') as fd:
                json.dump(mutations, fd)

        # Load variants dict
        with gzip.open(self.variants_dict_file, 'rt') as fd:
            variants_dict = json.load(fd)

        return variants_dict

    def _load_scores_dict(self, variants_dict):

        # Skip if done
        if not os.path.exists(self.scores_dict_file):

            # Initialize
            pool = Pool(self.cores)

            if not os.path.exists(self.score_file):
                raise RuntimeError("Score file not found '{}'".format(self.score_file))

            # Compute elements statistics
            logging.info("Computing score statistics by element")
            results = {}
            elements_all = list(variants_dict.items())
            elements_chunks = np.array_split(elements_all, self.cores)
            compute_arguments = ((chunk, self.regions_file, self.score_file, self.output_folder, num) for num, chunk in enumerate(elements_chunks, start=1))
            c = 0
            for done in pool.starmap(_compute_score_means, compute_arguments):
                c += 1
                logging.info("Chunk {} of {} [done]".format(c, self.cores))
                for element, means in done:
                    results[element] = means

            # Store results
            logging.info("Store scores dictionary")

            with gzip.open(self.scores_dict_file, 'wt') as fd:
                json.dump(results, fd)

        # Load scores dict
        with gzip.open(self.scores_dict_file, 'rt') as fd:
            scores_dict = json.load(fd)

        return scores_dict

    def _run_sampling(self, scores_dict):

        # Skip if done
        if os.path.exists(self.results_file):
            logging.info("Already calculated at '{}'".format(self.results_file))
            return

        # Prepare signature
        logging.info("Prepare signature files")
        signature_probabilities = pd.read_csv(self.signature_file, sep='\t')
        signature_probabilities.set_index(['Signature_reference', 'Signature_alternate'], inplace=True)
        signature_dict = signature_probabilities.to_dict()['Probability_' + _file_name(self.variants_file)]
        elements_all = [e for e, _ in scores_dict.items()]
        elements_chunks = np.array_split(elements_all, self.cores)
        arguments = [(chunk, self.output_folder, signature_dict, num) for num, chunk in enumerate(elements_chunks, start=1)]

        pool = Pool(self.cores)
        for c, _ in enumerate(pool.starmap(_create_background_signature, arguments), start=1):
            logging.info("Chunk {} of {} [done]".format(c, self.cores))

        # Calculate empirical p-values
        logging.info("Calculate empirical p-values")
        if len(scores_dict) == 0:
            logging.error("The scores file is empty")
            return

        num_randomizations = 10000
        to_run = [(element, means) for element, means in scores_dict.items()]
        while len(to_run) > 0:

            # Preprocess samplings
            next_to_run = []
            next_num_randomizations = min(1000000, num_randomizations * 2) if num_randomizations < 1000000 else 20000000

            logging.info("Start sampling ({})".format(num_randomizations))
            torun_chunks = np.array_split(to_run, self.cores)
            arguments = [(chunk, num_randomizations, self.output_folder, num) for num, chunk in enumerate(torun_chunks, start=1)]
            for result in pool.starmap(_sampling, arguments):
                for (e, m) in result:
                    # Check if we need more resolution at this element
                    if m['obs'] <= 5 and next_num_randomizations <= 1000000:
                        next_to_run.append((e, m))
                        logging.warning("We need more permutations at {}-{} (obs: {})".format(e, num_randomizations, m['obs']))
                    else:
                        scores_dict[e]['pvalue'] = m['pvalue']

            # Next iteration with the elements that need more permutations
            to_run = np.array(next_to_run)
            num_randomizations = next_num_randomizations

        # Run multiple test correction
        logging.info("Multiple test correction")
        results_concat = _multiple_test_correction(scores_dict, num_significant_samples=2)

        # Sort and store results
        results_concat.sort('pvalue', 0, inplace=True)
        fields = ['muts', 'muts_recurrence', 'samples_mut', 'samples_max_recurrence', 'subs', 'subs_score', 'pvalue', 'qvalue', 'all_mean', 'scores', 'positions']
        with open(self.results_file, 'wt') as fd:
            results_concat[fields].to_csv(fd, sep="\t", header=True, index=True)
        logging.info("Result save at: {}".format(self.results_file))


def cmdline():

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)

    # Parse the arguments
    parser = argparse.ArgumentParser()

    # Mandatory
    parser.add_argument('-i', '--input', dest='input_file', help='Variants file (maf, vcf or tab formated)')
    parser.add_argument('-r', '--regions', dest='regions_file', help='Genomic regions to analyse')
    parser.add_argument('-t', '--signature', dest='signature_file', help='Trinucleatide signature file')
    parser.add_argument('-s', '--score', dest='score_file', help='Tabix score file')

    # Optional
    parser.add_argument('-o', '--output', dest='output_folder', default='output', help='Output folder')
    parser.add_argument('-n', '--name', dest='project_name', default=None, help='Project name')
    parser.add_argument('--cores', dest='cores', default=os.cpu_count(), help="Maximum CPU cores to use (default all available)")

    args = parser.parse_args()
    logging.debug(args)

    # Initialize OncodriveFM2
    ofm2 = OncodriveFM2(
        args.input_file,
        args.regions_file,
        args.signature_file,
        args.score_file,
        args.output_folder,
        project_name=args.project_name,
        cores=args.cores
    )

    # Run
    ofm2.run()


if __name__ == "__main__":
    cmdline()