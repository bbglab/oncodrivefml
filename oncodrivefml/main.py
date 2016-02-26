import argparse

import gzip
import logging
import os
import pickle

from os.path import join, exists

from ago import human
from datetime import datetime

from oncodrivefml.config import load_configuration, file_exists_or_die, file_name
from oncodrivefml.executors import ElementExecutor
from oncodrivefml.load import load_and_map_variants
from oncodrivefml.mtc import multiple_test_correction
from oncodrivefml.store import store_tsv, store_png, store_html
from oncodrivefml.signature import load_signature

from multiprocessing.pool import Pool

from oncodrivefml.utils import executor_run


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
        results = {}
        element_executors = [ElementExecutor(element, muts, elements[element], signature, self.config) for element, muts in variants.items()]
        element_executors = [e for e in element_executors if len(e.muts) > 0]
        element_executors = sorted(element_executors, key=lambda e: -len(e.muts))

        i = 0
        start_time = datetime.now()
        logging.info("Computing OncodriveFML")

        logging.info("Run executors")
        for i, executor in enumerate(map_func(executor_run, element_executors)):
            if i % info_step == 0:
                logging.info("[{} of {}]".format(i+1, len(element_executors)))
            if len(executor.result['mutations']) > 0:
                results[executor.name] = executor.result
        logging.info("[{} of {}]".format(i+1, len(element_executors)))

        logging.debug("Computation time: {}".format(human(start_time)))

        if cores > 1:
            pool.close()

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
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')
    logging.getLogger().setLevel(logging.DEBUG if args.debug else logging.INFO)
    logging.debug(args)

    # Parse configuration file
    logging.info("Loading configuration")
    analysis = OncodriveFML(args.input_file, args.elements_file, args.output_folder, args.config_file, args.samples_blacklist)

    # Run the analysis
    analysis.run(debug=args.debug)

if __name__ == "__main__":
    cmdline()