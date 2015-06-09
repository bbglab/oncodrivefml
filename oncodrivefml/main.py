import argparse
import logging
import functools
import pandas as pd
import os

from multiprocessing.pool import Pool
from os.path import expanduser
from oncodrivefml import utils
from oncodrivefml.utils import _file_name, _silent_mkdir, _multiple_test_correction, _sampling, _load_variants_dict, _compute_element


class OncodriveFM2(object):

    def __init__(self, variants_file, regions_file, signature_file, score_file, output_folder,
                 project_name=None, cores=os.cpu_count(), cache=None, min_samplings=10000, max_samplings=1000000):

        # Configuration
        self.cores = cores
        self.min_samplings = min_samplings
        self.max_samplings = max_samplings
        self.variants_file = expanduser(variants_file)
        self.regions_file = expanduser(regions_file)

        if signature_file in ['compute', 'none']:
            self.signature_type = signature_file
            self.signature_field = None
            self.signature_file = None
        else:
            self.signature_type = 'file'
            signature_conf = signature_file.split(":")
            self.signature_field = signature_conf[0]
            self.signature_file = expanduser(signature_conf[1])

        self.score_file = expanduser(score_file)
        self.output_folder = expanduser(output_folder)
        self.project_name = project_name if project_name is not None else _file_name(variants_file)
        self.signature_name = _file_name(signature_file)
        self.results_file = os.path.join(output_folder, self.project_name + '-oncodrivefm2.tsv')
        self.cache = cache
        if cache is None:
            self.cache = os.path.join(self.output_folder, "cache")

        # Some initializations
        _silent_mkdir(self.cache)
        _silent_mkdir(output_folder)

    def run(self):

        # Skip if done
        if os.path.exists(self.results_file):
            logging.info("Already calculated at '{}'".format(self.results_file))
            return

        # Load variants
        variants_dict = _load_variants_dict(self.variants_file, self.regions_file, signature_name=self.signature_name)

        # Signature
        signature_dict = None
        if self.signature_type == "none":
            # We don't use signature
            logging.warning("We are not using any signature")
        elif self.signature_type == "compute":
            #TODO compute the signature
            logging.warning("TODO: compute signature. Running without signature.")
        else:
            if not os.path.exists(self.signature_file):
                logging.error("Signature file {} not found.".format(self.signature_file))
                return -1
            else:
                logging.info("Loading signature")
                signature_probabilities = pd.read_csv(self.signature_file, sep='\t')
                signature_probabilities.set_index(['Signature_reference', 'Signature_alternate'], inplace=True)
                signature_dict = signature_probabilities.to_dict()[self.signature_field]

        # Compute elements statistics
        logging.info("Computing statistics")
        elements = [(e, muts) for e, muts in variants_dict.items() if len(muts) > 0]
        if len(elements) == 0:
            logging.error("There is no mutation at any element")
            return -1

        results = {}
        info_step = 12*self.cores
        pool = Pool(self.cores)
        compute_element_partial = functools.partial(_compute_element, self.regions_file, self.score_file, self.cache, signature_dict, self.min_samplings, self.max_samplings)
        for i, (element, item) in enumerate(pool.imap(compute_element_partial, elements)):
            if i % info_step == 0:
                logging.info("[{} of {}]".format(i+1, len(elements)))

            if type(item) == dict:
                results[element] = item
            else:
                logging.debug(item)

        logging.info("[{} of {}]".format(i+1, len(elements)))

        # Run multiple test correction
        logging.info("Computing multiple test correction")
        results_concat = _multiple_test_correction(results, num_significant_samples=2)

        # Sort and store results
        results_concat.sort('pvalue', 0, inplace=True)
        fields = ['muts', 'muts_recurrence', 'samples_mut', 'pvalue', 'qvalue']
        with open(self.results_file, 'wt') as fd:
            results_concat[fields].to_csv(fd, sep="\t", header=True, index=True)
        logging.info("Done")

        return 1


def cmdline():

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)

    # Parse the arguments
    parser = argparse.ArgumentParser()

    # Mandatory
    parser.add_argument('-i', '--input', dest='input_file', help='Variants file (maf, vcf or tab formated)')
    parser.add_argument('-r', '--regions', dest='regions_file', help='Genomic regions to analyse')
    parser.add_argument('-t', '--signature', dest='signature_file', default="none", help='Trinucleotide signature file')
    parser.add_argument('-s', '--score', dest='score_file', help='Tabix score file')

    # Optional
    parser.add_argument('-o', '--output', dest='output_folder', default='output', help='Output folder')
    parser.add_argument('-n', '--name', dest='project_name', default=None, help='Project name')
    parser.add_argument('-mins', '--min_samplings', dest='min_samplings', type=int, default=10000, help="Minimum number of randomizations")
    parser.add_argument('-maxs', '--max_samplings', dest='max_samplings', type=int, default=100000, help="Maximum number of randomizations")
    parser.add_argument('--cores', dest='cores', type=int, default=os.cpu_count(), help="Maximum CPU cores to use (default all available)")
    parser.add_argument('--cache', dest='cache', default=None, help="Folder to store some intermediate data to speed up further executions.")

    args = parser.parse_args()
    logging.debug(args)

    # Check global configuration
    if utils.TABIX is None:
        logging.error("Cannot find 'tabix' executable. Please define the TABIX_PATH global variable.")
        exit(-1)

    if utils.HG19 is None:
        logging.error("Cannot find full genome files. Please define the FULL_GENOME_PATH global variable.")
        exit(-1)

    # Initialize OncodriveFM2
    ofm2 = OncodriveFM2(
        args.input_file,
        args.regions_file,
        args.signature_file,
        args.score_file,
        args.output_folder,
        project_name=args.project_name,
        cores=args.cores,
        cache=args.cache,
        min_samplings=args.min_samplings,
        max_samplings=args.max_samplings
    )

    #TODO allow only one score format
    utils.SCORE_CONF = utils.SCORES[os.path.basename(args.score_file)]

    # Run
    return_code = ofm2.run()

    if return_code != 1:
        exit(return_code)


if __name__ == "__main__":
    cmdline()