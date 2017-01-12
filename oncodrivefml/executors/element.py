import math
import logging
import numpy as np
from collections import Counter

from oncodrivefml.indels import Indel
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS


def flatten_partitions(results):

    for name, result in results.items():
        for partition in result['partitions']:
            yield (name, partition, result)


def compute_sampling(value):
    name, samples, result = value

    try:
        scores = result['simulation_scores']
        muts_count = result['muts_count']
        probs = result['simulation_probs']
        observed = result['observed']
        statistic_test = STATISTIC_TESTS.get(result['statistic_name'])
        background = np.random.choice(scores, size=(samples, muts_count), p=probs, replace=True)
        obs, neg_obs = statistic_test.calc_observed(background, np.array(observed))
    except Exception as e:
        logging.error("At {} error {}, result = {}".format(name, str(e), result.keys()))

    return name, obs, neg_obs


def partitions_list(total_size, chunk_size):
    """
    Create a list of values less or equal to chunk_size that sum total_size

    :param total_size: Total size
    :param chunk_size: Chunk size
    :return: list of integers
    """
    partitions = [chunk_size for _ in range(total_size // chunk_size)]

    res = total_size % chunk_size
    if res != 0:
        partitions += [res]

    return partitions


class ElementExecutor(object):
    """
    Base class for executors that do the analysis per genomic element.
    The mutations are simulated to compute a p_value
    by comparing the functional impact of the observed mutations
    and the expected.
    The simulation parameters are taken from the configuration file.

    Args:
        element_id (str): element ID
        muts (list): list of mutations belonging to that element (see :ref:`mutations <mutations dict>` inner items)
        segments (list): list of segments belonging to that element (see :ref:`elements <elements dict>` inner items)
        signature (dict): probabilities of each mutation (see :ref:`signatures <signature dict>`)
        config (dict): configurations

    """

    def __init__(self, element_id, muts, segments, signature, config):
        # Input attributes
        self.name = element_id
        self.indels_conf = config['statistic']['indels']
        self.use_indels = self.indels_conf['enabled']
        self.use_subs = config['statistic']['subs']
        self.use_mnp = config['statistic']['mnp'] and self.use_subs
        # MNP mutations are only used if subs are enabled

        if self.use_subs and self.use_indels and self.use_mnp:
            self.muts = muts
        else:
            self.muts = []
            if self.use_subs:
                self.muts += [m for m in muts if m['TYPE'] == 'subs']
            if self.use_mnp:
                self.muts += [m for m in muts if m['TYPE'] == 'mnp']
            if self.use_indels:
                self.muts += [m for m in muts if m['TYPE'] == 'indel']

        self.signature = signature
        self.segments = segments
        self.symbol = self.segments[0].get('SYMBOL', None)
        self.is_positive_strand = False if segments[0].get('STRAND', '+') == '-' else True
        # When the strand is unknown is considered the same as positive
        #TODO fix

        # Configuration parameters
        self.score_config = config['score']
        self.sampling_size = config['statistic']['sampling']
        self.min_obs = config['statistic']['sampling_min_obs']
        self.sampling_chunk = config['statistic']['sampling_chunk']
        self.statistic_name = config['statistic']['method']
        self.signature_column = config['signature']['classifier']
        self.samples_method = config['statistic']['per_sample_analysis']

        # Output attributes
        self.obs = 0
        self.neg_obs = 0
        self.result = None
        self.scores = None

        self.p_subs = config['p_subs']
        self.p_indels = config['p_indels']

    def compute_muts_statistics(self, indels=False):
        """
        Assigns scores to the mutations and retreive the ones that are going
        to be simulated

        Args:
            indels (:obj:`~oncodrivefml.indels.Indel`): Indels class if indels are considered. False otherwise

        Returns:
            dict: several information about the mutations and a list of them with the scores

        """
        raise RuntimeError("The classes that extend ElementExecutor must override, at least the compute_muts_statistics() method")

    def compute_mutation_score(self, mutation, indels):

        # Get substitutions scores
        if mutation['TYPE'] == "subs":
            mutation['POSITION'] = int(mutation['POSITION'])
            values = self.scores.get_score_by_position(mutation['POSITION'])
            for v in values:
                if v.ref == mutation['REF'] and v.alt == mutation['ALT']:
                    mutation['SCORE'] = v.value
                    break
            else:
                logging.warning('Discrepancy in SUBS at position {} of chr {}'.format(mutation['POSITION'], mutation['CHROMOSOME']))

        if mutation['TYPE'] == "mnp":
            pos = int(mutation['POSITION'])
            mnp_scores = []
            for index, nucleotides in enumerate(zip(mutation['REF'], mutation['ALT'])):
                ref_nucleotide, alt_nucleotide = nucleotides
                values = self.scores.get_score_by_position(pos + index)
                for v in values:
                    if v.ref == ref_nucleotide and v.alt == alt_nucleotide:
                        mnp_scores.append(v.value)
                        break
                else:
                    logging.warning('Discrepancy in MNP at position {} of chr {}'.format(pos, mutation['CHROMOSOME']))
            if not mnp_scores:
                return
            mutation['SCORE'] = max(mnp_scores)
            mutation['POSITION'] = pos + mnp_scores.index(mutation['SCORE'])

        if indels is not False and mutation['TYPE'] == "indel":

            # very long indels are discarded
            if max(len(mutation['REF']), len(mutation['ALT'])) > 20:
                return

            score = indels.get_indel_score(mutation)

            mutation['SCORE'] = score if not math.isnan(score) else None

    def run(self):
        """
        Loads the scores and compute the statistics for the observed mutations.
        For all positions around the mutation position, gets the scores and
        probabilities of mutations in those positions.
        Generates a random set of mutations (for each mutations, randomizes
        certain amount).
        Combining those random values (make the transpose to get a matrix
        number_of_real_amount_of_mutation x randomization_amount) computes
        how many are above and how many below the value of the current
        mutation (statistical value of all mutations).
        Computes the p-values.

        :return: p values
        """

        # Load element scores
        self.scores = Scores(self.name, self.segments, self.score_config)

        if self.use_indels:
            indels = Indel(self.scores, self.signature, self.signature_column,
                                self.indels_conf['method'], self.is_positive_strand)
        else:
            indels = False

        # Compute observed mutations statistics and scores
        self.result = self.compute_muts_statistics(indels=indels)

        if len(self.result['mutations']) > 0:
            statistic_test = STATISTIC_TESTS.get(self.statistic_name)

            positions = self.scores.get_all_positions()

            observed = []
            subs_scores = []
            subs_probs = []
            indels_scores = []
            indels_probs = []

            subs_probs_by_signature = {}
            signature_ids = []

            indels_simulated_as_subs = 0

            for mut in self.result['mutations']:
                observed.append(mut['SCORE'])  # Observed mutations
                if mut['TYPE'] == 'subs' or mut['TYPE'] == 'mnp':
                    if self.signature is not None:
                        # Count how many signature ids are and prepare a vector for each
                        # IMPORTANT: this implies that only the signature of the observed mutations is taken into account
                        signature_ids.append(mut.get(self.signature_column, self.signature_column))

                # Indels treated as subs also count for the subs probs
                elif mut['TYPE'] == 'indel' and indels.simulated_as_subs:
                    self.p_subs = 1
                    self.p_indels = 0
                    if self.signature is not None:
                        if self.indels_conf['indels_simulated_with_signature']:
                            signature_ids.append(mut.get(self.signature_column, self.signature_column))
                        else:
                            signature_ids.append('indels_having_no_signature')
                # When only in frame indels are simulated as subs
                elif mut['TYPE'] == 'indel' and indels.in_frame_simulated_as_subs:
                    if max(len(mut['REF']), len(mut['ALT'])) % 3 == 0:
                        indels_simulated_as_subs += 1
                    if self.signature is not None:
                        if self.indels_conf['indels_simulated_with_signature']:
                            signature_ids.append(mut.get(self.signature_column, self.signature_column))
                        else:
                            signature_ids.append('indels_having_no_signature')

            for signature_id in set(signature_ids):
                subs_probs_by_signature.update({signature_id: []})

            # When the probabilities of subs and indels are None, they are taken from
            # mutations seen in the gene
            if self.p_subs is None or self.p_indels is None: # use the probabilities based on observed mutations
                self.p_subs = (self.result['subs'] + indels_simulated_as_subs) / len(self.result['mutations'])
                self.p_indels = 1 - self.p_subs

            # Compute the values for the substitutions
            if self.use_subs and self.p_subs > 0:
                for pos in positions:
                    for s in self.scores.get_score_by_position(pos):
                        subs_scores.append(s.value)
                        for k, v in subs_probs_by_signature.items():
                            if k in self.signature.keys():
                                v.append(self.signature[k].get((s.ref_triplet, s.alt_triplet), 0.0))
                            else:
                                v.append(1/192)

                if len(subs_probs_by_signature) > 0:
                    signature_ids_counter = Counter(signature_ids)
                    total_ids = len(signature_ids)
                    subs_probs = np.array([0.0] * len(subs_scores))
                    for k, v in subs_probs_by_signature.items():
                        subs_probs += (np.array(v) * signature_ids_counter[k] / total_ids)
                    tot = sum(subs_probs)
                    subs_probs = subs_probs * self.p_subs / tot
                    subs_probs = list(subs_probs)
                else: # Same prob for all values (according to the prob of substitutions
                    subs_probs = [self.p_subs / len(subs_scores)] * len(subs_scores) if len(subs_scores) > 0 else []

            if self.use_indels and self.p_indels > 0:
                indels_scores = indels.get_background_indel_scores(mutations=[m for m in self.result['mutations'] if m['TYPE'] == 'indel'])

                # All indels have the same probability
                indels_probs = [self.p_indels/len(indels_scores)] * len(indels_scores) if len(indels_scores) > 0 else []

            simulation_scores = subs_scores + indels_scores
            simulation_probs = subs_probs + indels_probs

            muts_count = len(self.result['mutations'])
            chunk_size = (self.sampling_size * muts_count) // self.sampling_chunk
            chunk_size = self.sampling_size if chunk_size == 0 else self.sampling_size // chunk_size

            # Calculate sampling parallelization partitions
            self.result['partitions'] = partitions_list(self.sampling_size, chunk_size)
            self.result['sampling'] = {}

            # Run first partition
            first_partition = self.result['partitions'].pop(0)
            background = np.random.choice(simulation_scores, size=(first_partition, muts_count), p=simulation_probs, replace=True)
            obs, neg_obs = statistic_test.calc_observed(background, np.array(observed))
            self.result['obs'] = obs
            self.result['neg_obs'] = neg_obs

            # Sampling parallelization (if more than one partition)
            if len(self.result['partitions']) > 1 or obs < self.min_obs:
                self.result['muts_count'] = muts_count
                self.result['simulation_scores'] = simulation_scores
                self.result['simulation_probs'] = simulation_probs
                self.result['observed'] = observed
                self.result['statistic_name'] = self.statistic_name

        self.result['sampling_size'] = self.sampling_size
        self.result['symbol'] = self.symbol

        return self
