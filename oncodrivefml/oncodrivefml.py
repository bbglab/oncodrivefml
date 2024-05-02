"""
Contains the main class of the method
"""
import gzip
import io
import os
import sys
import csv
import logging
import warnings
import json
from os.path import join, exists
from multiprocessing.pool import Pool

import numpy as np

from oncodrivefml import __version__, __logger_name__, signature, load, reference
from oncodrivefml.config import file_exists_or_die, file_name
from oncodrivefml.executors.bymutation import GroupByMutationExecutor
from oncodrivefml.executors.bysample import GroupBySampleExecutor
from oncodrivefml.mtc import multiple_test_correction
from oncodrivefml.groups import add_groups, add_groups_from_values
from oncodrivefml.scores import init_scores_module
from oncodrivefml.mutability import init_mutabilities_module
from oncodrivefml.depths import init_depths_module
from oncodrivefml.store import store_tsv, store_png, store_html, store_scores_tsv, store_groups_tsv
from oncodrivefml.utils import executor_run, loop_logging # , load_depths, load_mutability
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

    def __init__(self, mutations_file, elements_file, output_folder, config, blacklist):
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

        self.cores = self.configuration['settings']['cores']
        if self.cores is None:
            self.cores = os.cpu_count()
        logger.debug('Using %s cores', self.cores)

        # TODO decide whether this is provided in the config file or in the commmand line
        # so far I am providing it in the command line
        self.store_bckg = self.configuration['settings']['store_bckg']
        if self.store_bckg:
            logger.debug('The background means are going to be stored.')

        if self.configuration['signature']['method'] == 'bysample':
            self.configuration['signature']['method'] = 'complement'
            self.configuration['signature']['classifier'] = 'SAMPLE'

        if self.configuration['statistic']['indels']['include'] and \
                'simulate_with_signature' in self.configuration['statistic']['indels']:
            warnings.warn('"simulate_with_signature" option deprecated. '
                          'Indels simulated as substitutions follow the signature',
                          DeprecationWarning)

        self.samples_statistic_method = self.configuration['statistic']['per_sample_analysis']

        # grouping
        if self.configuration['grouping']['group_genes']:
            file_exists_or_die(self.configuration['grouping']['json_file'])
            logger.info("After processing the genes individually, OncodriveFML will compute the results per groups of genes.")

        # mutability
        if self.configuration['mutability']['adjusting']:
            file_exists_or_die(self.configuration['mutability']['file'])
            logger.info("Mutability per site will be used for adjusting the background probabilities of SNVs.")
            logger.info("Signature and depth values will be ignored for SNVs.")
            self.configuration['signature']['method'] = 'none'

            self.configuration['mutability_info'] = True
        else:
            self.configuration['mutability_info'] = False

        # depths
        if self.configuration['depth']['adjusting']:
            file_exists_or_die(self.configuration['depth']['depth_file'])
            logger.info("Depth per site will be used for adjusting the background probabilities.")
            # self.configuration['depths_loaded'] = load_depths(self.configuration['depth']['depth_file'],
            #                                                     self.configuration['depth']['chr_prefix'])
            # logger.info("Depths file loaded")

            self.configuration['depth_info'] = True
        else:
            self.configuration['depth_info'] = False


        ## Report additional information usage status
        # both depths and mutabilities
        if self.configuration['depth_info'] and not self.configuration['mutability_info'] :
            logger.info("Indels background probability will be adjusted by the sequencing depth.")
            logger.info("SNVs background probability will be adjusted by the mutabilities.")

        # mutabilities only
        elif not self.configuration['depth_info'] and self.configuration['mutability_info'] :
            logger.info("Indels background probability will be uniform.")
            logger.info("SNVs background probability will be adjusted by the mutabilities.")

        # depths only
        elif self.configuration['depth_info'] and not self.configuration['mutability_info'] :
            logger.info("Indels background probability will be adjusted by the sequencing depth.")
            logger.info("SNVs background probability will be adjusted by the sequencing depth and the mutational profile.")

        # none of them available
        else:
            logger.info("Neither depth nor mutabilities provided.")
            self.configuration['mutability_info'] = None
            self.configuration['depth_info'] = None

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

        np.random.seed(self.configuration['settings']['seed'])

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

        if self.configuration['score'].get('minimum_number_of_stops', None) is not None:
            warnings.warn('minimum_number_of_stops should be specified in the indels section',
                          DeprecationWarning)
        else:
            self.configuration['score']['minimum_number_of_stops'] = self.configuration['statistic']['indels']['minimum_number_of_stops']

        init_scores_module(self.configuration['score'], stops_required=stops_file_required)
        if self.configuration['mutability_info']:
            init_mutabilities_module(self.configuration['mutability'])
        
        if self.configuration['depth_info']:
            init_depths_module(self.configuration['depth'])

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
                for name, obs, neg_obs, back_means, internal_values in loop_logging(map_func(compute_sampling, partitions), size=len(partitions), step=1):
                    result = results[name]
                    result['obs'] += obs
                    result['neg_obs'] += neg_obs
                    
                    # print(obs, neg_obs, back_mean)
                    if type(result['back_means']) == list:
                        result['back_means'] = result['back_means'] + back_means
                    else:
                        result['back_means'] = back_means

                    if type(result['internal_values']) == list:
                        result['internal_values'] = result['internal_values'] + internal_values
                    else:
                        result['internal_values'] = internal_values

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
            logger.warning("Empty results, possible reason: no mutation from the dataset can be mapped to the provided regions.")
            sys.exit(0)


        # Compute z-score
        logger.info("Compute z-score")
        for gene in results.keys():
            result = results[gene]
            mean_of_observed = np.mean(result['scores'])
            std_distribution_of_means = result['population_std'] / np.sqrt( len(result['scores']) )

            result['std_distribution_of_means'] = std_distribution_of_means

            result['z-score'] = round( ( mean_of_observed - result['population_mean'] ) / std_distribution_of_means, 5)


        if self.configuration['grouping']['group_genes']:
            # Compute groups scores
            logger.info("Compute groups")

            with open(self.configuration['grouping']['json_file'], 'rt') as f:
                groups_dict = json.load(f)
            results_groups, results_groups_n_indv = add_groups_from_values(results, groups_dict)


            # TODO decide if the mtc has to be done for the groups alone
            # or with all the genes
            results_all_mtc = multiple_test_correction(results_groups_n_indv, num_significant_samples=2)

        # Run multiple test correction
        logger.info("Computing multiple test correction")
        
        # TODO revise the concept of num_significant_samples for normal tissues where we assume a same sample can have multiple clones
        # revise it also above
        results_mtc = multiple_test_correction(results, num_significant_samples=2)

        
        # TODO see if it would be worth doing a .store() method 
        # that runs this block of code that has nothing to do with the computations
        # and we could include the part that generates the plots

        # Sort and store results
        logger.info("Storing results")
        if not exists(self.output_folder):
            os.makedirs(self.output_folder, exist_ok=True)
        result_file = self.output_file_prefix + '.tsv.gz'
        store_tsv(results_mtc, result_file)

        if self.configuration['grouping']['group_genes']:
            result_groups_file = self.output_file_prefix + '.groups.tsv.gz'
            store_groups_tsv(results_all_mtc, result_groups_file, include_indv = False)
            

        if self.store_bckg:
            logger.info("Storing background means")
            scoresresult_file = self.output_file_prefix + '.bckg.tsv.gz'
            store_scores_tsv(results_mtc, scoresresult_file)

        lines = 0
        gene_ids = {None, ''}
        with gzip.open(result_file, 'rt') as csvfile:
            fd = csv.DictReader(csvfile, delimiter='\t')
            for line in fd:
                lines += 1
                gene_ids.add(line['GENE_ID'])
        if lines+2 != len(gene_ids):
            logger.error('Number of genes does not match number of lines in the output file. Please check the logs of the execution to find more information.')

        if self.output_folder is not None:
            logger.info("Creating figures")
            store_png(result_file, self.output_file_prefix + ".png")
            store_html(result_file, self.output_file_prefix + ".html")

        logger.info("Done")