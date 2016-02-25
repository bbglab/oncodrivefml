import argparse
import gzip
import logging
import os
import pickle

from os.path import join, exists

from oncodrivefml.config import load_configuration, file_exists_or_die, file_name
from oncodrivefml.executors import ElementExecutor
from oncodrivefml.load import load_and_map_variants
from oncodrivefml.mtc import multiple_test_correction
from oncodrivefml.store import store_tsv, store_png, store_html
from oncodrivefml.signature import load_signature

from multiprocessing.pool import Pool

from oncodrivefml.utils import run_executor


class OncodriveFML(object):

    def __init__(self, input_file, elements_file, output_folder, config_file, blacklist):

        # Required parameters
        self.input_file = file_exists_or_die(input_file)
        self.elements_file= file_exists_or_die(elements_file)
        self.config = load_configuration(config_file)
        self.blacklist = blacklist

        # Optional parameters
        self.output_folder = file_name(self.elements_file) if output_folder is None else output_folder
        self.output_file_prefix = join(self.output_folder, file_name(self.input_file) + '-oncodrivefml')

    def run(self, debug=False):

        # Skip if done
        if exists(self.output_file_prefix + '.tsv'):
            logging.warning("Already calculated at '{}'".format(self.output_file_prefix + '.tsv'))
            return

        # Load variants mapping
        variants, elements = load_and_map_variants(self.input_file, self.elements_file, blacklist=self.blacklist)

        # Load signature
        signature = load_signature(self.input_file, self.config['signature'], blacklist=self.blacklist)

        # Run in parallel
        cores = self.config['settings'].get('cores')
        if cores is None:
            cores = os.cpu_count()
        info_step = 6*cores

        if cores > 1:
            pool = Pool(cores)
            map_func = pool.imap
        else:
            map_func = map

        # Run the executors
        element_executors_computed = {}
        element_executors = [ElementExecutor(element, muts, elements[element], signature, self.config) for element, muts in variants.items()]
        logging.info("Computing OncodriveFML")
        i = 0
        for i, executor in enumerate(map_func(run_executor, element_executors)):
            if i % info_step == 0:
                logging.info("[{} of {}]".format(i+1, len(element_executors)))
            element_executors_computed[executor.name] = executor

        logging.info("[{} of {}]".format(i+1, len(element_executors)))

        if cores > 1:
            pool.close()

        # Results
        results = {e.name: e.result for e in element_executors_computed.values()}

        # Run multiple test correction
        logging.info("Computing multiple test correction")
        results_mtc = multiple_test_correction(results, num_significant_samples=2)

        # Sort and store results
        logging.info("Storing results")
        if not exists(self.output_folder):
            os.makedirs(self.output_folder)
        result_file = self.output_file_prefix + '.tsv'
        store_tsv(results_mtc, result_file)

        logging.info("Creating figures")
        store_png(result_file, self.output_file_prefix + ".png")
        store_html(result_file, self.output_file_prefix + ".html")

        if debug:
            element_executors_file = self.output_file_prefix + "_executors.pickle.gz"
            with gzip.open(element_executors_file, 'wb') as fd:
                pickle.dump(element_executors_computed, fd)

        logging.info("Done")
        return 0


def cmdline():

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
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG if args.debug else logging.INFO)
    logging.debug(args)

    # Parse configuration file
    logging.info("Loading configuration")
    analysis = OncodriveFML(args.input_file, args.elements_file, args.output_folder, args.config_file, args.samples_blacklist)

    # Run the analysis
    analysis.run(debug=args.debug)

if __name__ == "__main__":
    cmdline()