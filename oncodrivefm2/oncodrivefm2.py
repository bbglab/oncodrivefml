import argparse
import logging
import os
from configuration import Configuration
from oncodrivefm2_step1 import OncodriveFM2Scoring
from oncodrivefm2_step2 import OncodriveFM2Test
from signaturesampling import SignatureSampling


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class OncodriveFM2(object):
    def __init__(self, fmconfig):
        self._config = fmconfig
        self._sampling = SignatureSampling(fmconfig)

    def get_scoring(self):
        return OncodriveFM2Scoring(self._config, self._sampling)

    def get_test(self):
        return OncodriveFM2Test(self._config, self._sampling)

    def get_sampling(self):
        return self._sampling

    def get_config(self):
        return self._config


def run(config: Configuration, step0: bool, step1: bool, step2: bool):

    fm2 = OncodriveFM2(config)

    logger.info("SIGNATURE is: {}".format(config.signature))

    mapping = None
    scoring = None
    test = None

    if step0:
        # annotatation of variants to gene and segments
        raise Exception("step 0 unimplemented as of yet")

    # annotation of the score of the variants
    if step1:

        # Check if ensembl file is specified
        if config.ensembl_file is None:
            raise Exception("In order to run OncodriveFM2-Scoring (step 1) please supply a correct ensembl file")

        # Check input configuration
        if not os.path.isdir(config.input):
            raise Exception("The input for step 1 must be a directory containing the " + config.MUTS_FILE_NAME +
                            " and " + config.INDELS_FILE_NAME + "files")

        # Check if indel file is supplied for user-input indel-strategy
        if config.indel_strategy == config.INDEL_STRATEGY_USER_INPUT and config.indel_file == None:
            raise Exception("Indel file (file with scores for indels) must be supplied for USER INPUT indel-strategy")

        # Run Scoring (Step 1)
        scoring = fm2.get_scoring()
        scoring.run()

        discarded = scoring.get_discarded()
        logger.info("{} discarded elements: {}".format(len(discarded), discarded))


    if step2:
        # sampling and test

        test = fm2.get_test()

        if (scoring != None):
            config.input = scoring.get_output_file()

        if not os.path.isfile(config.input):
            raise Exception("The input for step 2 must be a file containing the sampling results.")

        test.run()

def cmdline():
    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')

    # Define the steps of the method
    step0 = 'variant-mapping'
    step1 = 'score-calc'
    step2 = 'test'
    test_steps = [step0, step1, step2]

    # Parse the arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("--cores", type=int, help="Number of CPUs to be used")
    parser.add_argument("--project", type=str, help="Project identifier/name")
    parser.add_argument("--tissue", type=str, help="Cancer type code")
    parser.add_argument("--feature", type=str, help="Feature: CDS, UTR, etc.")
    parser.add_argument("--signature", default="", type=str, help="Signature of mutation probabilities (by default 'project_tissue'), "
                                                      "'none' if not desired.")
    parser.add_argument('--indel-strategy', type=str, choices=Configuration.INDEL_STRATEGIES,
                        help="Choose a strategy for calculating indel score")
    parser.add_argument('--step', type=str, choices=test_steps, action='append', help="Three steps in the following order"
                                                                                      "have to be performed: " + ", ".join(test_steps))
    parser.add_argument("--score", type=str, help="Define desired score(s)")
    parser.add_argument("--input", dest="input_dir_or_file", type=str, help="Input file or directory, depending on which step is called.")
    parser.add_argument("--output-dir", type=str, help="Folder of intermediate and final results.")
    parser.add_argument("--cache-dir", type=str, help="Folder where variant-score mapping are cached.")
    parser.add_argument("--ensembl-file", type=str, help="Needed in step 1 (score-calc).", required=False)
    parser.add_argument("--indel-file", type=str, help="Needed in step 1 for user-input indel strategy.", required=False)
    parser.add_argument("--expression-file", type=str, help="Needed in step 2 (test).", required=False)
    parser.add_argument("--testing", help="Run a test : input files will be cut.", action="store_true", required=False)

    args = parser.parse_args()
    print(args)

    config = Configuration(project=args.project, tissue=args.tissue, score=args.score, feature=args.feature.lower(),
                           output_dir=args.output_dir, cache_dir=args.cache_dir, input_dir_or_file=args.input_dir_or_file,
                           ensembl_file=args.ensembl_file, cores=args.cores, indel_strategy=args.indel_strategy, testing=args.testing,
                           indel_file=args.indel_file, signature=args.signature)

    run(config, step0 in args.step, step1 in args.step, step2 in args.step)


if __name__ == "__main__":
    cmdline()