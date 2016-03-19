import argparse

import logging
import os

from os.path import join, exists
from oncodrivefml.config import load_configuration, file_exists_or_die, file_name
from oncodrivefml.executors.bymutation import GroupByMutationExecutor
from oncodrivefml.executors.bysample import GroupBySampleExecutor
from oncodrivefml.load import load_and_map_variants
from oncodrivefml.mtc import multiple_test_correction
from oncodrivefml.store import store_tsv, store_png, store_html
from oncodrivefml.signature import load_signature
from multiprocessing.pool import Pool
from oncodrivefml.utils import executor_run, loop_logging


class OncodriveFML(object):

    def __init__(self, input_file, elements_file, output_folder, config_file, blacklist):
        """
        Initialize OncodriveFML analysis

        :param input_file: Mutations input file
        :param elements_file: Genomic element input file
        :param output_folder: Folder where the results will be store
        :param config_file: Configuration file
        :param blacklist: File with sample ids (one per line) to remove when loading the input file
        """

        # Required parameters
        self.input_file = file_exists_or_die(input_file)
        self.elements_file= file_exists_or_die(elements_file)
        self.config = load_configuration(config_file)
        self.blacklist = blacklist
        self.cores = self.config['settings']['cores']
        if self.cores is None:
            self.cores = os.cpu_count()
        self.statistic_method = self.config['statistic']['method']

        # Optional parameters
        self.output_folder = file_name(self.elements_file) if output_folder is None else output_folder
        self.output_file_prefix = join(self.output_folder, file_name(self.input_file) + '-oncodrivefml')

        # Output parameters
        self.variants = None
        self.elements = None
        self.signature = None

    def create_element_executor(self, element, muts):

        if self.statistic_method == 'maxmean':
            return GroupBySampleExecutor(element, muts, self.elements[element], self.signature, self.config)

        return GroupByMutationExecutor(element, muts, self.elements[element], self.signature, self.config)

    def run(self):
        """
        Run the OncodriveFML analysis
        """

        # Skip if done
        if exists(self.output_file_prefix + '.tsv'):
            logging.warning("Already calculated at '{}'".format(self.output_file_prefix + '.tsv'))
            return

        # Load variants mapping
        self.variants, self.elements = load_and_map_variants(self.input_file, self.elements_file, blacklist=self.blacklist)

        # Load signature
        self.signature = load_signature(self.input_file, self.config['signature'], blacklist=self.blacklist)

        # Create one executor per element
        element_executors = [self.create_element_executor(element, muts) for element, muts in self.variants.items()]

        # Remove elements that don't have mutations to compute
        element_executors = [e for e in element_executors if len(e.muts) > 0]

        # Sort executors to compute first the ones that have more mutations
        element_executors = sorted(element_executors, key=lambda e: -len(e.muts))

        # Run the executors
        with Pool(self.cores) as pool:
            results = {}
            logging.info("Computing OncodriveFML")
            for executor in loop_logging(pool.imap(executor_run, element_executors), size=len(element_executors), step=6*self.cores):
                if len(executor.result['mutations']) > 0:
                    results[executor.name] = executor.result

        # Run multiple test correction
        logging.info("Computing multiple test correction")
        results_mtc = multiple_test_correction(results, num_significant_samples=2)

        # Sort and store results
        logging.info("Storing results")
        if not exists(self.output_folder):
            os.makedirs(self.output_folder, exist_ok=True)
        result_file = self.output_file_prefix + '.tsv'
        store_tsv(results_mtc, result_file)

        logging.info("Creating figures")
        store_png(result_file, self.output_file_prefix + ".png")
        store_html(result_file, self.output_file_prefix + ".html")

        logging.info("Done")


def cmdline():

    # Command line arguments parser
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument('-i', '--input', dest='input_file', required=True, help='Variants file')
    parser.add_argument('-e', '--elements', dest='elements_file', required=True, help='Genomic elements to analyse')

    # Optional arguments
    parser.add_argument('-o', '--output', dest='output_folder', default=None, help="Output folder. Default to regions file name without extensions.")
    parser.add_argument('-c', '--config', dest='config_file', default=None, help="Configuration file. Default to 'oncodrivefml.conf' in the current folder if exists or to ~/.bbglab/oncodrivefml.conf if not.")
    parser.add_argument('--samples-blacklist', dest='samples_blacklist', default=None, help="Remove this samples when loading the input file")
    parser.add_argument('--debug', dest='debug', default=False, action='store_true', help="Show more progress details")

    # Parse arguments
    args = parser.parse_args()

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')
    logging.getLogger().setLevel(logging.DEBUG if args.debug else logging.INFO)
    logging.debug(args)

    # Load configuration file and prepare analysis
    logging.info("Loading configuration")
    analysis = OncodriveFML(args.input_file, args.elements_file, args.output_folder, args.config_file, args.samples_blacklist)

    # Run the analysis
    analysis.run()

if __name__ == "__main__":
    cmdline()