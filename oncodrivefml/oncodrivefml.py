"""
Contains the main class of the method
"""

import io
import os
import sys
import csv
import logging
import warnings
from os.path import join, exists
from multiprocessing.pool import Pool

import numpy as np

from oncodrivefml import __version__, __logger_name__, signature, load, reference
from oncodrivefml.config import file_exists_or_die, file_name
from oncodrivefml.executors.bymutation import GroupByMutationExecutor
from oncodrivefml.executors.bysample import GroupBySampleExecutor
from oncodrivefml.mtc import multiple_test_correction
from oncodrivefml.scores import init_scores_module
from oncodrivefml.store import store_tsv, store_png, store_html
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

    """

    def __init__(self, mutations_file, elements_file, output_folder, config, blacklist, cores, seed):
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

        genome_reference_build = self.configuration['genome']['build']
        reference.change_build(genome_reference_build)

        self.cores = cores if cores is not None else self.configuration['settings']['cores']
        if self.cores is None:
            self.cores = os.cpu_count()
        logger.debug('Using %s cores', self.cores)

        if self.configuration['signature']['method'] == 'bysample':
            self.configuration['signature']['method'] = 'complement'
            self.configuration['signature']['classifier'] = 'SAMPLE'

        if self.configuration['statistic']['indels']['include'] and \
                'simulate_with_signature' in self.configuration['statistic']['indels']:
            warnings.warn('"simulate_with_signature" option deprecated. '
                          'Indels simulated as substitutions follow the signature',
                          DeprecationWarning)

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

        np.random.seed(seed)

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
        seed = np.random.randint(0, 2**32-1)
        if self.samples_statistic_method is None:
            return GroupByMutationExecutor(element_id, element_mutations, self.elements[element_id], self.signatures,
                                           self.configuration, seed)
        else:
            return GroupBySampleExecutor(element_id, element_mutations, self.elements[element_id], self.signatures,
                                         self.configuration, seed)

    def __compute_signature(self):
        conf = self.configuration['signature']

        if conf.get('only_mapped_mutations', False):
            warnings.warn('only_mapped_mutations for signature configuration is no longer used. Build your signature as a separate file.', DeprecationWarning)

        if conf.get('include_mnp', False):
            warnings.warn(
                'include_mnp for signature configuration is no longer used',
                DeprecationWarning)

        if conf['method'] == 'none':
            logger.warning('Not using any signature')
            self.signatures = None

        elif conf['method'] == 'file':
            if conf['normalize_by_sites'] is not None:
                logger.warning('Cannot normalize precomputed signatures. Normalization of signatures omitted.')

            self.signatures = signature.load(conf['path'])

        else:
            self.signatures = signature.compute(load.snp(self.mutations_file, blacklist=self.blacklist),
                                               conf['method'], conf['classifier'], conf['normalize_by_sites'])

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

        # Load mutations mapping
        mutations_data, self.elements = load.mutations_and_elements(self.mutations_file,
                                                               self.elements_file,
                                                               blacklist=self.blacklist,
                                                                indels_max_size=self.configuration['statistic']['indels'].get('max_size', None))

        self.mutations = mutations_data['data']

        self.__compute_simulation_probs(mutations_data['metadata'])

        # Load signatures
        self.__compute_signature()

        if self.configuration['statistic']['indels']['include'] and \
            self.configuration['statistic']['indels']['method'] == 'stop':
            stops_file_required = True
        else:
            stops_file_required = False
        init_scores_module(self.configuration['score'], stops_required=stops_file_required)

        # initialize the indels module
        init_indels_module(self.configuration['statistic']['indels'])

        # Create one executor per element
        element_executors = [self.__create_element_executor(element_id, muts) for
                             element_id, muts in sorted(self.mutations.items())]

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
                logger.info("Parallel sampling. Iteration %d, genes %d, partitions %d", i, len(set([n for n, p, r, s in partitions])), len(partitions))

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
                            partitions.append((name, partition, result, np.random.randint(0, 2**32-1)))

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
