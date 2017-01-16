"""
Contains the command line parsing and the main class of the method
"""


import os
import sys
import argparse
import logging
from multiprocessing.pool import Pool
from os.path import join, exists

from oncodrivefml.config import load_configuration, file_exists_or_die, file_name
from oncodrivefml.executors.bymutation import GroupByMutationExecutor
from oncodrivefml.executors.bysample import GroupBySampleExecutor
from oncodrivefml.load import load_and_map_variants, load_mutations, count_mutations
from oncodrivefml.mtc import multiple_test_correction
from oncodrivefml.store import store_tsv, store_png, store_html
from oncodrivefml.signature import load_signature, yield_mutations, change_ref_build
from oncodrivefml.utils import executor_run, loop_logging
from oncodrivefml.indels import _init_indels
from oncodrivefml.walker import flatten_partitions, compute_sampling, partitions_list


class OncodriveFML(object):
    """

    Args:
       mutations_file: Mutations input file (see :mod:`~oncodrivefml.load` for details)
       elements_file: Genomic element input file (see :mod:`~oncodrivefml.load` for details)
       output_folder: Folder where the results will be stored
       configuration_file: Configuration file (see :ref:`configuration <project configuration>`)
       blacklist: File with sample ids (one per line) to remove when loading the input file
       pickle_save (bool): save pickle files

    """

    def __init__(self, mutations_file, elements_file, output_folder, configuration_file, blacklist, pickle_save):

        # Required parameters
        self.mutations_file = file_exists_or_die(mutations_file)
        self.elements_file = file_exists_or_die(elements_file)
        self.configuration = load_configuration(configuration_file)
        self.blacklist = blacklist
        self.save_pickle = pickle_save

        genome_reference_build = self.configuration['genome']['build']
        change_ref_build(genome_reference_build)

        self.cores = self.configuration['settings']['cores']
        if self.cores is None:
            self.cores = os.cpu_count()
        self.samples_statistic_method = self.configuration['statistic']['per_sample_analysis']

        if self.configuration['signature']['method'] == 'bysample':
            self.configuration['signature']['method'] = 'complement'
            self.configuration['signature']['classifier'] = 'SAMPLE'

        # Optional parameters
        self.output_folder = file_name(self.elements_file) if output_folder is None else output_folder
        self.output_file_prefix = join(self.output_folder, file_name(self.mutations_file) + '-oncodrivefml')

        # Output parameters
        self.mutations = None
        self.elements = None
        self.signatures = None

        self.debug = False



    def create_element_executor(self, element_id, muts_for_an_element):
        """
        To enable paralelization, for each element ID,
        one :class:`~oncodrivefml.executors.bymutation.ElementExecutor` is created.

        Args:
            element_id (str): ID of the element
            muts_for_an_element (list): list with all the mutations observed in the element

        Returns:
            :class:`~oncodrivefml.executors.bymutation.ElementExecutor`:
            returns :class:`~oncodrivefml.executors.bymutation.GroupByMutationExecutor` if
            the there is no method to apply to the samples;
            :class:`~oncodrivefml.executors.bysamplen.GroupBySampleExecutor` otherwise.


        """
        if self.samples_statistic_method is not None:
            return GroupBySampleExecutor(element_id, muts_for_an_element, self.elements[element_id], self.signatures,
                                         self.configuration)

        return GroupByMutationExecutor(element_id, muts_for_an_element, self.elements[element_id], self.signatures,
                                       self.configuration)

    def run(self):
        """
        Run the OncodriveFML analysis.
        """

        # TODO add explanation

        # Skip if done
        if exists(self.output_file_prefix + '.tsv'):
            logging.warning("Already calculated at '{}'".format(self.output_file_prefix + '.tsv'))
            return

        # Load mutations mapping
        self.mutations, self.elements = load_and_map_variants(self.mutations_file,
                                                              self.elements_file,
                                                              blacklist=self.blacklist,
                                                              save_pickle=self.save_pickle)

        if self.configuration['statistic']['use_gene_mutations']:
            self.configuration['p_indels'] = None
            self.configuration['p_subs'] = None
        else:
            if self.configuration['statistic']['subs'] and self.configuration['statistic']['indels']['enabled']:
                # In case we are using indels and subs. Ohterwise it is pointless to get the counts of each
                subs_counter, indels_counter = count_mutations(self.mutations_file, blacklist=self.blacklist)

                if self.configuration['statistic']['indels']['method'] == 'stop' and self.configuration['statistic']['indels']['enable_frame']:
                    logging.info('Indels identified as in frame, are going to be simulated as subs')
                    # count how many indels are frameshift
                    discarded = 0
                    for element, mutations in self.mutations.items():
                        for mut in mutations:
                            if mut['TYPE'] == 'indel':
                                indel_size = max(len(mut['REF']), len(mut['ALT']))
                                if indel_size % 3 == 0:
                                    discarded +=1
                    indels_counter -= discarded
                    subs_counter += discarded

                p_indels = indels_counter/(indels_counter + subs_counter)
                p_subs = subs_counter/(subs_counter + indels_counter)

                self.configuration['p_indels'] = p_indels
                self.configuration['p_subs'] = p_subs
            else:
                self.configuration['p_indels'] = 1
                self.configuration['p_subs'] = 1


        if self.configuration['signature']['use_only_mapped_mutations']:
            signature_function = lambda: yield_mutations(self.mutations)

        else:
            signature_function = lambda: load_mutations(self.mutations_file, show_warnings=False, blacklist=self.blacklist)


        save_signature_pickle = self.save_pickle
        if save_signature_pickle:
            if self.blacklist is not None or self.configuration['signature']['use_only_mapped_mutations']:
                save_signature_pickle = False
                logging.warning('Signature pickle not saved because a blacklist for the mutations was provided or the use_only_mapped_mutations flag was set to true')

        load_signature_pickle = True if self.blacklist is None and not self.configuration['signature']['use_only_mapped_mutations']  else False

        # Load signatures
        self.signatures = load_signature(self.mutations_file, signature_function, self.configuration['signature'],
                                         save_pickle=save_signature_pickle,
                                         load_pickle=load_signature_pickle)

        # Create one executor per element
        element_executors = [self.create_element_executor(element_id, muts) for
                             element_id, muts in self.mutations.items()]

        # Remove elements that don't have mutations to compute
        element_executors = [e for e in element_executors if len(e.muts) > 0]

        # Sort executors to compute first the ones that have more mutations
        element_executors = sorted(element_executors, key=lambda e: -len(e.muts))

        # initialize the indels module
        indels_config = self.configuration['statistic']['indels']
        if indels_config['enabled']:
            _init_indels(indels_config)

        # Run the executors
        with Pool(self.cores) as pool:
            results = {}
            logging.info("Computing OncodriveFML")
            map_func = pool.imap if not self.debug else map
            for executor in loop_logging(map_func(executor_run, element_executors), size=len(element_executors), step=6*self.cores):
                if len(executor.result['mutations']) > 0:
                    results[executor.name] = executor.result

            # Flatten partitions
            partitions = list(flatten_partitions(results))

            i = 0
            while len(partitions) > 0 or i == 0:

                i += 1
                logging.info("Parallel sampling. Iteration {}, genes {}, partitions {}".format(i, len(set([n for n,p,r in partitions])), len(partitions)))

                # Pending sampling execution
                for name, obs, neg_obs in loop_logging(map_func(compute_sampling, partitions), size=len(partitions), step=1):
                    result = results[name]
                    result['obs'] += obs
                    result['neg_obs'] += neg_obs

                # Increase sampling_size
                partitions = []
                for name, result in results.items():
                    sampling_size = float(result['sampling_size'])

                    if result['obs'] < self.configuration['statistic']['sampling_min_obs'] and sampling_size < self.configuration['statistic']['sampling_max']:
                        next_sampling_size = min(sampling_size*10, self.configuration['statistic']['sampling_max'])
                        pending_sampling_size = int(next_sampling_size - sampling_size)
                        chunk_size = (pending_sampling_size * result['muts_count']) // self.configuration['statistic']['sampling_chunk']
                        chunk_size = pending_sampling_size if chunk_size == 0 else pending_sampling_size // chunk_size

                        result['sampling_size'] = next_sampling_size
                        for partition in partitions_list(pending_sampling_size, chunk_size):
                            partitions.append((name, partition, result))

            # Compute p-values
            logging.info("Compute p-values")
            for result in results.values():
                sampling_size = float(result['sampling_size'])
                result['pvalue'] = max(1, result['obs']) / sampling_size
                result['pvalue_neg'] = max(1, result['neg_obs']) / sampling_size

        if results == {}:
            logging.warning("Empty resutls, possible reason: no mutation from the dataset can be mapped to the provided regions.")
            sys.exit(0)

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
    """
    Parses the command and runs the analysis. See :meth:`~OncodriveFML.run`.

    Required arguments:

    -i, --input file        mutations file
    -e, --elements file     elements file

    Optional arguments:

    -o, --output folder     output folder
    -c, --configuration file
                            configuration file
    --samples-blacklist file
                            file with blacklisted samples
    --save-pickle           store intermediate information in pickle files
    --debug                 show more progress

    """

    # Command line arguments parser
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument('-i', '--input', dest='mutations_file', required=True, help='Variants file')
    parser.add_argument('-e', '--elements', dest='elements_file', required=True, help='Genomic elements to analyse')

    # Optional arguments
    parser.add_argument('-o', '--output', dest='output_folder', default=None, help="Output folder. Default to regions file name without extensions.")
    parser.add_argument('-c', '--configuration', dest='config_file', default=None, help="Configuration file. Default to 'oncodrivefml.conf' in the current folder if exists or to ~/.bbglab/oncodrivefml.conf if not.")
    parser.add_argument('--samples-blacklist', dest='samples_blacklist', default=None, help="Remove this samples when loading the input file")
    parser.add_argument('--save-pickle', dest='save_pickle', action='store_true', default=False, help="Save intermediate information as pickle files.")
    parser.add_argument('--debug', dest='debug', default=False, action='store_true', help="Show more progress details")

    # Parse arguments
    args = parser.parse_args()

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')
    logging.getLogger().setLevel(logging.DEBUG if args.debug else logging.INFO)
    logging.debug(args)

    # Load configuration file and prepare analysis
    logging.info("Loading configuration")
    analysis = OncodriveFML(args.mutations_file, args.elements_file, args.output_folder, args.config_file,
                            args.samples_blacklist, args.save_pickle)


    # Run the analysis
    analysis.run()

if __name__ == "__main__":
    cmdline()
