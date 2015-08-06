import argparse
import gzip
import logging
import functools
import pickle
import os

from multiprocessing.pool import Pool
from os.path import expanduser
from configobj import ConfigObj
from validate import Validator
from oncodrivefml import signature
from oncodrivefml.drmaa import drmaa_run
from oncodrivefml.qqplot import qqplot_png, qqplot_html
from oncodrivefml.qqplot import add_symbol
from oncodrivefml.compute import file_name, silent_mkdir, multiple_test_correction, compute_element
from oncodrivefml.load import load_variants_dict, load_regions
from oncodrivefml.signature import load_signature


SCORES = {
    'whole_genome_SNVs.tsv.gz': {
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5, 'element': None
    },
    'hg19_wg_score.tsv.gz': {
        'chr': 0, 'chr_prefix': 'chr', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'hg19_rnasnp_scores.txt.gz': {
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 5, 'element': 6
    },
    'tfbs_creation.tsv.gz': {
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'tfbs_disruption.tsv.gz': {
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'disruption_v2.txt.bgz': {
        'chr': 0, 'chr_prefix': 'chr', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None, 'extra': 5
    }
}


class OncodriveFML(object):

    def __init__(self, variants_file, regions_file, signature_file, score_file, output_folder,
                 project_name=None, cores=os.cpu_count(), min_samplings=10000, max_samplings=1000000,
                 max_jobs=100, debug=False, trace=None, geometric=False, score_conf=None):

        # Configuration
        self.cores = cores
        self.debug = debug
        self.trace = [] if trace is None else trace
        self.max_jobs = max_jobs
        self.min_samplings = min_samplings
        self.max_samplings = max_samplings
        self.variants_file = expanduser(variants_file)
        self.regions_file = expanduser(regions_file)
        self.geometric = geometric

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
        self.score_conf = score_conf
        self.score_conf['file'] = expanduser(score_conf['file'])
        self.output_folder = expanduser(output_folder)
        self.project_name = project_name if project_name is not None else file_name(variants_file)
        self.signature_name = file_name(signature_file)
        self.results_file = os.path.join(output_folder, self.project_name + '-oncodrivefml.tsv')
        self.qqplot_file = os.path.join(output_folder, self.project_name + '-oncodrivefml')

        # Some initializations
        silent_mkdir(output_folder)

    def run(self, drmaa=None, figures=True):

        # Skip if done
        if os.path.exists(self.results_file):
            logging.info("Already calculated at '{}'".format(self.results_file))
            return

        # Load regions
        logging.info("Loading regions")
        regions = load_regions(self.regions_file)

        # Load variants
        logging.info("Loading and mapping mutations")
        variants_dict = load_variants_dict(self.variants_file, regions, signature_name=self.signature_name)

        # Signature
        signature_dict = load_signature(self.variants_file, self.signature_file, self.signature_field, self.signature_type, self.signature_name)

        # Run in a DRMAA cluster
        if drmaa is not None:
            return drmaa_run(variants_dict, signature_dict, self, drmaa, figures=figures)

        # Compute elements statistics
        logging.info("Computing statistics")
        elements = []
        for e, muts in variants_dict.items():
            if len(muts) > 0:
                if e in self.trace:
                    trace = os.path.join(self.output_folder, "{}.pickle.gz".format(e))
                else:
                    trace = None
                elements.append((e, muts, regions[e], trace))

        if len(elements) == 0:
            logging.error("There is no mutation at any element")
            return -1

        results = {}
        info_step = 6*self.cores
        pool = Pool(self.cores)

        compute_element_partial = functools.partial(compute_element, signature_dict, self.min_samplings, self.max_samplings, self.geometric, self.score_conf)

        all_missing_signatures = {}
        for i, (element, item, missing_signatures) in enumerate(pool.imap(compute_element_partial, elements)):
            if i % info_step == 0:
                logging.info("[{} of {}]".format(i+1, len(elements)))

            if type(item) == dict:
                results[element] = item
            else:
                logging.debug(item)

            for ref, alt in missing_signatures:
                all_missing_signatures[ref] = alt

        logging.info("[{} of {}]".format(i+1, len(elements)))

        if len(all_missing_signatures) > 0:
            logging.warning("There are positions without signature probability. We are using a score of zero at these positions.")
            logging.warning("If you are computing the signature from the input file, most probable this means that you don't have enough mutations.")
            logging.warning("Try using a precomputed signature of a similar cancer type, or don't use signature.")
            logging.warning("The missing signatures are:")
            for ref_triplet, alt_triplet in all_missing_signatures.items():
                logging.warning("\tref: '%s' alt: '%s'", ref_triplet, alt_triplet)

        # Store partial result
        if self.variants_file.endswith(".pickle.gz"):
            logging.info("Store partial result")
            with gzip.open(os.path.join(self.output_folder, self.project_name + '.pickle.gz'), 'wb') as fd:
                pickle.dump(results, fd)
            logging.info("Done")
            return 0

        # Run multiple test correction
        logging.info("Computing multiple test correction")
        results_concat = multiple_test_correction(results, num_significant_samples=2)

        # Sort and store results
        results_concat.sort('pvalue', 0, inplace=True)
        fields = ['muts', 'muts_recurrence', 'samples_mut', 'pvalue', 'qvalue']
        df = results_concat[fields].copy()
        df.reset_index(inplace=True)
        df = add_symbol(df)
        with open(self.results_file, 'wt') as fd:
            df.to_csv(fd, sep="\t", header=True, index=False)

        if figures:
            logging.info("Creating figures")
            qqplot_png(self.results_file, self.qqplot_file + ".png")
            qqplot_html(self.results_file, self.qqplot_file + ".html")

        logging.info("Done")
        return 0


def cmdline():

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
    parser.add_argument('--debug', dest='debug', default=False, action='store_true')
    parser.add_argument('--no-figures', dest='no_figures', default=False, action='store_true')
    parser.add_argument('--drmaa', dest='drmaa', type=int, default=None, help="Run in a DRMAA cluster using this value as the number of elements to compute per job.")
    parser.add_argument('--drmaa-max-jobs', dest='drmaa_max_jobs', type=int, default=100, help="Maximum parallell concurrent jobs")
    parser.add_argument('--trace', dest='trace', nargs='+', type=str, default=None, help="Elements IDs to store files to trace and reproduce the execution")
    parser.add_argument('--geometric', dest='geometric', default=False, action='store_true')
    args = parser.parse_args()

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG if args.debug else logging.INFO)
    logging.debug(args)

    if signature.HG19 is None:
        logging.error("Cannot find full genome files. Please define the FULL_GENOME_PATH global variable.")
        exit(-1)

    # Allow scores with different formats
    if args.score_file.endswith(".conf"):
        score_conf = ConfigObj(args.score_file, configspec={
            'file': 'string', 'chr': 'integer', 'chr_prefix': 'string', 'pos': 'integer', 'ref': 'integer',
            'alt': 'integer', 'score': 'integer', 'element': 'integer', 'extra': 'integer'
        })
        score_conf.validate(Validator(), preserve_errors=True)
    else:
        # TODO allow only one score format or move this to external configuration
        score_conf = SCORES.get(os.path.basename(args.score_file), SCORES['whole_genome_SNVs.tsv.gz'])
        score_conf['file'] = args.score_file

    # Initialize OncodriveFM2
    ofm2 = OncodriveFML(
        args.input_file,
        args.regions_file,
        args.signature_file,
        args.score_file,
        args.output_folder,
        project_name=args.project_name,
        cores=args.cores,
        min_samplings=args.min_samplings,
        max_samplings=args.max_samplings,
        max_jobs=args.drmaa_max_jobs,
        debug=args.debug,
        trace=args.trace,
        geometric=args.geometric,
        score_conf=score_conf
    )

    # Run
    return_code = ofm2.run(drmaa=args.drmaa, figures=not args.no_figures)

    if return_code != 0:
        exit(return_code)


if __name__ == "__main__":
    cmdline()