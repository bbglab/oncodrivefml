"""
Contains the command line parsing and the main class of the method
"""

import io
import os
import sys
import csv
import click
import logging
from os.path import join, exists
from collections import defaultdict
from multiprocessing.pool import Pool
from logging.config import dictConfig

from oncodrivefml import __version__, __logger_name__
from oncodrivefml.config import load_configuration, file_exists_or_die, file_name
from oncodrivefml.executors.bymutation import GroupByMutationExecutor
from oncodrivefml.executors.bysample import GroupBySampleExecutor
from oncodrivefml.load import load_and_map_variants, load_mutations
from oncodrivefml.mtc import multiple_test_correction
from oncodrivefml.scores import init_scores_module
from oncodrivefml.store import store_tsv, store_png, store_html
from oncodrivefml.signature import load_signature, yield_mutations, change_ref_build, load_trinucleotides_counts, \
    compute_regions_signature
from oncodrivefml.utils import executor_run, loop_logging
from oncodrivefml.indels import init_indels_module
from oncodrivefml.walker import flatten_partitions, compute_sampling, partitions_list

logger = logging.getLogger(__logger_name__)

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


class OncodriveFML(object):
    """

    Args:
       mutations_file: Mutations input file (see :mod:`~oncodrivefml.load` for details)
       elements_file: Genomic element input file (see :mod:`~oncodrivefml.load` for details)
       output_folder: Folder where the results will be stored
       config: configuration (see :ref:`configuration <project configuration>`)
       blacklist: File with sample ids (one per line) to remove when loading the input file
       generate_pickle (bool): flag to only generate pickle files and exit

    """

    def __init__(self, mutations_file, elements_file, output_folder, config, blacklist, generate_pickle):
        logger.debug('Using OncodriveFML version %s', __version__)

        # Required parameters
        self.mutations_file = file_exists_or_die(mutations_file)
        logger.debug('Mutations file: %s', self.mutations_file)
        self.elements_file = file_exists_or_die(elements_file)
        logger.debug('Elements file: %s', self.elements_file)
        self.configuration = config
        self.blacklist = blacklist
        if self.blacklist is not None:
            logger.debug('Blacklist file: %s', self.blacklist)
            logger.debug('Using a blacklist causes some pickle files not to be saved/loaded')
        self.generate_pickle = generate_pickle

        genome_reference_build = self.configuration['genome']['build']
        change_ref_build(genome_reference_build)

        self.cores = self.configuration['settings']['cores']
        if self.cores is None:
            self.cores = os.cpu_count()
        logger.debug('Using %s cores', self.cores)

        if self.configuration['signature']['method'] == 'bysample':
            self.configuration['signature']['method'] = 'complement'
            self.configuration['signature']['classifier'] = 'SAMPLE'

        self.samples_statistic_method = self.configuration['statistic']['per_sample_analysis']

        # Optional parameters
        self.output_folder = output_folder
        self.output_file_prefix = join(self.output_folder, file_name(self.mutations_file) + '-oncodrivefml')
        logger.debug('Output: %s', self.output_folder)

        # Output parameters
        self.mutations = None
        self.elements = None
        self.signatures = None

        self.avoid_parallel = False

        s = io.BytesIO()
        self.configuration.write(s)
        logger.debug('Configuration used:\n' + s.getvalue().decode())
        s.close()

    def __create_element_executor(self, element_id, element_mutations):
        """
        To enable parallelization, for each element ID,
        one :class:`~oncodrivefml.executors.element.ElementExecutor` is created.

        Args:
            element_id (str): ID of the element
            element_mutations (list): list with all the mutations observed in the element

        Returns:
            :class:`~oncodrivefml.executors.bymutation.GroupByMutationExecutor` or
            :class:`~oncodrivefml.executors.bysample.GroupBySampleExecutor`.

        """
        if self.samples_statistic_method is None:
            return GroupByMutationExecutor(element_id, element_mutations, self.elements[element_id], self.signatures,
                                           self.configuration)
        else:
            return GroupBySampleExecutor(element_id, element_mutations, self.elements[element_id], self.signatures,
                                         self.configuration)

    def __compute_signature(self):
        conf = self.configuration['signature']

        precomputed_signature = self.mutations_file + '_signature_' + conf['method'] + '_' + conf['classifier'] + ".pickle.gz"
        load_signature_pickle = precomputed_signature if self.blacklist is None else None
        save_signature_pickle = precomputed_signature if self.generate_pickle else None

        trinucleotides_counts = None
        if conf['only_mapped_mutations']:
            signature_function = lambda: yield_mutations(self.mutations)
            if conf['normalize_by_sites'] is not None:
                trinucleotides_counts = compute_regions_signature(self.elements.values(), self.cores)
        else:
            signature_function = lambda: load_mutations(self.mutations_file, blacklist=self.blacklist)
            if conf['normalize_by_sites'] is not None:
                trinucleotides_counts = load_trinucleotides_counts(conf['normalize_by_sites'])

        # Load signatures
        self.signatures = load_signature(conf, signature_function, trinucleotides_counts,
                                         load_pickle=load_signature_pickle,
                                         save_pickle=save_signature_pickle)

    def __compute_simulation_probs(self, counts):
        """
        In the simulation mutations are simulated as either substitutions or stops.

        Args:
            counts (dict): counts for each mutation type

        """
        # compute cohort percentages of subs and indels
        if self.configuration['statistic']['indels']['gene_exomic_frameshift_ratio']:
            self.configuration['p_indels'] = None
            self.configuration['p_subs'] = None
        else:
            if self.configuration['statistic']['indels']['include']:
                # subs_counter = counts['snp']
                subs_counter = counts['snp_mapped']

                if not self.configuration['statistic']['discard_mnp']:
                    # subs_counter += counts['mnp']
                    subs_counter += counts['mnp_mapped']

                # indels_counter = counts['indel']
                indels_counter = counts['indels_mapped'] - counts['indels_mapped_multiple_of_3']
                subs_counter += counts['indels_mapped_multiple_of_3']

                try:
                    p_indels = indels_counter / (indels_counter + subs_counter)
                    p_subs = 1 - p_indels
                    self.configuration['p_indels'] = p_indels
                    self.configuration['p_subs'] = p_subs
                except ZeroDivisionError:
                    logger.warning('Impossible to compute relative counts of indels and substitutions')
                    self.configuration['p_indels'] = None
                    self.configuration['p_subs'] = None

            else:
                self.configuration['p_indels'] = 0
                self.configuration['p_subs'] = 1

    def run(self):
        """
        Run the OncodriveFML analysis.
        """

        # TODO add explanation

        # Load mutations mapping
        mutations_data, self.elements = load_and_map_variants(self.mutations_file,
                                                              self.elements_file,
                                                              blacklist=self.blacklist,
                                                              save_pickle=self.generate_pickle)

        self.mutations = mutations_data['data']

        self.__compute_simulation_probs(mutations_data['metadata'])

        # Load signatures
        self.__compute_signature()

        if self.generate_pickle:
            logger.info('Pickles generated. Exiting')
            return

        init_scores_module(self.configuration['score'])

        # initialize the indels module
        init_indels_module(self.configuration['statistic']['indels'])

        # Create one executor per element
        element_executors = [self.__create_element_executor(element_id, muts) for
                             element_id, muts in self.mutations.items()]

        # Remove elements that don't have mutations to compute
        element_executors = [e for e in element_executors if len(e.muts) > 0]

        # Sort executors to compute first the ones that have more mutations
        element_executors = sorted(element_executors, key=lambda e: -len(e.muts))

        # Run the executors
        with Pool(self.cores) as pool:
            results = {}
            logger.info("Computing OncodriveFML")
            map_func = map if self.avoid_parallel else pool.imap
            for executor in loop_logging(map_func(executor_run, element_executors), size=len(element_executors), step=6*self.cores):
                if len(executor.result['mutations']) > 0:
                    results[executor.name] = executor.result

            # Flatten partitions
            partitions = list(flatten_partitions(results))

            i = 0
            while len(partitions) > 0 or i == 0:

                i += 1
                logger.info("Parallel sampling. Iteration %d, genes %d, partitions %d", i, len(set([n for n,p,r in partitions])), len(partitions))

                # Pending sampling execution
                for name, obs, neg_obs in loop_logging(map_func(compute_sampling, partitions), size=len(partitions), step=1):
                    result = results[name]
                    result['obs'] += obs
                    result['neg_obs'] += neg_obs

                # Increase sampling_size
                partitions = []
                for name, result in results.items():
                    sampling_size = float(result['sampling_size'])

                    if result['obs'] is not None and result['obs'] < self.configuration['statistic']['sampling_min_obs'] and sampling_size < self.configuration['statistic']['sampling_max']:
                        next_sampling_size = min(sampling_size*10, self.configuration['statistic']['sampling_max'])
                        pending_sampling_size = int(next_sampling_size - sampling_size)
                        chunk_count = (pending_sampling_size * result['muts_count']) // (self.configuration['statistic']['sampling_chunk'] * 10**6)
                        chunk_size = pending_sampling_size if chunk_count == 0 else pending_sampling_size // chunk_count

                        result['sampling_size'] = next_sampling_size
                        for partition in partitions_list(pending_sampling_size, chunk_size):
                            partitions.append((name, partition, result))

            # Compute p-values
            logger.info("Compute p-values")
            for result in results.values():
                sampling_size = float(result['sampling_size'])
                result['pvalue'] = max(1, result['obs']) / sampling_size if result['obs'] is not None else None
                result['pvalue_neg'] = max(1, result['neg_obs']) / sampling_size if result['neg_obs'] is not None else None

        if results == {}:
            logger.warning("Empty resutls, possible reason: no mutation from the dataset can be mapped to the provided regions.")
            sys.exit(0)

        # Run multiple test correction
        logger.info("Computing multiple test correction")
        results_mtc = multiple_test_correction(results, num_significant_samples=2)

        # Sort and store results
        logger.info("Storing results")
        if not exists(self.output_folder):
            os.makedirs(self.output_folder, exist_ok=True)
        result_file = self.output_file_prefix + '.tsv'
        store_tsv(results_mtc, result_file)

        lines = 0
        gene_ids = {None, ''}
        with open(result_file) as csvfile:
            fd = csv.DictReader(csvfile, delimiter='\t')
            for line in fd:
                lines += 1
                gene_ids.add(line['GENE_ID'])
        if lines+2 != len(gene_ids):
            logger.error('Number of genes does not match number of lines in the output file. Please check the logs of the execution to find more information.')

        logger.info("Creating figures")
        store_png(result_file, self.output_file_prefix + ".png")
        store_html(result_file, self.output_file_prefix + ".html")

        logger.info("Done")


def main(mutations_file, elements_file, output_folder, config_file, samples_blacklist, generate_pickle, config_override_dict=None):
    """
    Run OncodriveFML analysis

    Args:
        mutations_file (str): path to the mutations file
        elements_file (str): path to the elements file
        output_folder (str): path to the output folder. Set to :obj:`None` to create
           a folder with the elements file name in the current directory
        config_file (str): path to configuration file
        samples_blacklist (str): path to samples blacklist file. Set to :obj:`None`
           if you are not using any blacklist file
        generate_pickle (bool): whether run OncodriveFML to generate pickle files or full analysis
        config_override_dict (dict, optional): override configuration from file

    """

    configuration = load_configuration(config_file, override=config_override_dict)

    output_folder = file_name(elements_file) if output_folder is None else output_folder
    output_file = join(output_folder, file_name(mutations_file) + '-oncodrivefml.tsv')
    # Skip if done
    if exists(output_file):
        logging.warning("Already calculated at '{}'".format(output_file))
        return
    else:
        if not exists(output_folder):
            os.makedirs(output_folder, exist_ok=True)
        if 'logging' in configuration:
            if configuration['logging'].get('handlers', {}).get('file', None) is not None:
                logging_file = join(output_folder, file_name(mutations_file) + '__log.txt')
                configuration['logging']['handlers']['file']['filename'] = logging_file
            dictConfig(configuration['logging'])

    analysis = OncodriveFML(mutations_file, elements_file, output_folder, configuration,
                            samples_blacklist, generate_pickle)

    logger.info('Running analysis')
    # Run the analysis
    analysis.run()




@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input', 'mutations_file', type=click.Path(exists=True), help='Variants file', metavar='MUTATIONS_FILE',required=True)
@click.option('-e', '--elements', 'elements_file', type=click.Path(exists=True), metavar='ELEMENTS_FILE', help='Genomic elements to analyse', required=True)
@click.option('-t', '--type', type=click.Choice(['coding', 'noncoding']), help='Type of genomic elements file', required=True)
@click.option('-s', '--sequencing', type=click.Choice(['wgs', 'wes', 'targeted']), help='Type of sequencing: whole genome, whole exome or targeted.', required=True)
@click.option('-o', '--output', 'output_folder', type=click.Path(), metavar='OUTPUT_FOLDER', help="Output folder. Default to regions file name without extensions.", default=None)
@click.option('-c', '--configuration', 'config_file', default=None, type=click.Path(exists=True), metavar='CONFIG_FILE', help="Configuration file. Default to 'oncodrivefml_v2.conf' in the current folder if exists or to ~/.bbglab/oncodrivefml.conf if not.")
@click.option('--samples-blacklist', default=None, type=click.Path(exists=True), metavar='SAMPLES_BLACKLIST', help="Remove these samples when loading the input file.")
@click.option('--no-indels', help="Discard indels in your analysis", is_flag=True)
@click.option('--generate-pickle', help="Run OncodriveFML to generate pickle files that could speed up future executions and exit.", is_flag=True)
@click.option('--debug', help="Show more progress details", is_flag=True)
@click.version_option(version=__version__)
def cmdline(mutations_file, elements_file, type, sequencing, output_folder, config_file, samples_blacklist, no_indels, generate_pickle, debug):
    """
    Run OncodriveFML on the genomic regions in ELEMENTS FILE
    using the mutations in MUTATIONS FILE.


    """

    dd = lambda: defaultdict(dd)
    override_config = dd()

    # Overrride the configuration
    if no_indels:
        override_config['statistic']['indels']['include'] = False
    else:
        override_config['statistic']['indels']['include'] = True
    if type == 'coding':
        override_config['statistic']['indels']['method'] = 'stop'
    elif type == 'noncoding':
        override_config['statistic']['indels']['method'] = 'max'

    if sequencing == 'wex':
        override_config['signature']['normalize_by_sites'] = 'whole_exome'
    elif sequencing == 'wgs':
        override_config['signature']['normalize_by_sites'] = 'whole_genome'
    else:
        override_config['signature']['normalize_by_sites'] = None

    if debug:
        override_config['logging']['handlers']['console']['level'] = 'DEBUG'
    else:
        override_config['logging']['handlers']['console']['level'] = 'INFO'

    main(mutations_file, elements_file, output_folder, config_file, samples_blacklist, generate_pickle, override_config)


if __name__ == "__main__":
    cmdline()
