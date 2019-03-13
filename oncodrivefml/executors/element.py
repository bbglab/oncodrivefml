import math
import logging
import numpy as np

from oncodrivefml import __logger_name__
from oncodrivefml.executors import sig2probs
from oncodrivefml.indels import Indel
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS
from oncodrivefml.walker import partitions_list

logger = logging.getLogger(__logger_name__)


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

    def __init__(self, element_id, muts, segments, signature, config, seed):
        # Input attributes
        self.name = element_id
        self.seed = seed

        self.use_mnp = not config['statistic']['discard_mnp']
        self.indels_conf = config['statistic']['indels']
        self.use_indels = self.indels_conf['include']
        self.indels = None

        if self.use_indels and self.use_mnp:
            self.muts = muts
        else:
            self.muts = [m for m in muts if m['ALT_TYPE'] == 'snp']
            if self.use_mnp:
                self.muts += [m for m in muts if m['ALT_TYPE'] == 'mnp']
            if self.use_indels:
                self.muts += [m for m in muts if m['ALT_TYPE'] == 'indel']

        self.signature = signature
        self.segments = segments
        self.symbol = self.segments[0].get('SYMBOL', None)

        # Configuration parameters
        self.score_config = config['score']
        self.sampling_size = config['statistic']['sampling']
        self.min_obs = config['statistic']['sampling_min_obs']
        self.sampling_chunk = config['statistic']['sampling_chunk'] * 10**6
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

    def compute_muts_statistics(self):
        """
        Assigns scores to the mutations and retrieve the ones that are going
        to be simulated

        Returns:
            dict: several information about the mutations and a list of them with the scores

        """
        raise RuntimeError("The classes that extend ElementExecutor must override, at least the compute_muts_statistics() method")

    def compute_mutation_score(self, mutation):

        # Get substitutions scores
        if mutation['ALT_TYPE'] == "snp":
            mutation['POSITION'] = int(mutation['POSITION'])
            values = self.scores.get_score_by_position(mutation['POSITION'])
            for v in values:
                if v.ref == mutation['REF'] and v.alt == mutation['ALT']:
                    mutation['SCORE'] = v.value
                    break
            else:
                logger.warning('Reference mismatch or missing score in SNP at position {} of chr {}'.format(mutation['POSITION'], mutation['CHROMOSOME']))

        elif mutation['ALT_TYPE'] == "mnp":
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
                    logger.warning('Reference mismatch or missing score in MNP at position {} of chr {}'.format(pos, mutation['CHROMOSOME']))
            if not mnp_scores:
                return
            mutation['SCORE'] = max(mnp_scores)
            mutation['POSITION'] = pos + mnp_scores.index(mutation['SCORE'])

        elif mutation['ALT_TYPE'] == "indel":

            score = self.indels.get_indel_score(mutation)

            mutation['SCORE'] = score if not math.isnan(score) else None

    def run(self):
        """
        Loads the scores and compute the statistics for the observed mutations.
        Perform simulations to compute a set of randomized scores.
        By comparing those scores with the observed through different
        statistical tests, a p-value can be obtained.

        """
        # as this method in run as a separate process, we need to fix a seed here too
        np.random.seed(self.seed)

        # Load element scores
        self.scores = Scores(self.name, self.segments, self.score_config)

        if self.use_indels:
            self.indels = Indel(self.scores)

        # Compute observed mutations statistics and scores
        self.result = self.compute_muts_statistics()

        if len(self.result['mutations']) > 0:
            statistic_test = STATISTIC_TESTS.get(self.statistic_name)

            positions = self.scores.get_all_positions()

            observed = []
            subs_scores = []
            subs_probs = []

            indels_scores = []
            indels_probs = []

            probs_of_subs = sig2probs.build(self.signature, self.signature_column)

            indels_simulated_as_subs = 0

            for mut in self.result['mutations']:
                observed.append(mut['SCORE'])  # Observed mutations

                # Indels treated as subs also count for the subs probs
                if mut['ALT_TYPE'] == 'indel' and self.indels.simulated_as_subs:
                    self.p_subs = 1
                    self.p_indels = 0
                    probs_of_subs.add_observed(mut)
                # When only in frame indels are simulated as subs
                elif mut['ALT_TYPE'] == 'indel' and self.indels.in_frame_simulated_as_subs:
                    if max(len(mut['REF']), len(mut['ALT'])) % 3 == 0:
                        indels_simulated_as_subs += 1
                        probs_of_subs.add_observed(mut)
                elif mut['ALT_TYPE'] == 'indel':
                    pass
                else:  # SNP or MNP
                    # Count how many signature ids are and prepare a vector for each
                    # IMPORTANT: this implies that only the signature of the observed mutations is taken into account
                    probs_of_subs.add_observed(mut)

            # When the probabilities of subs and indels are None, they are taken from
            # mutations seen in the gene
            if self.p_subs is None or self.p_indels is None: # use the probabilities based on observed mutations
                self.p_subs = (self.result['snps'] + self.result['mnps'] + indels_simulated_as_subs) / len(self.result['mutations'])
                self.p_indels = 1 - self.p_subs

            # Compute the values for the substitutions
            if self.p_subs > 0:
                for pos in positions:
                    for s in self.scores.get_score_by_position(pos):
                        subs_scores.append(s.value)
                        probs_of_subs.add_background(s.change)

                if probs_of_subs.size > 0:
                    subs_probs = probs_of_subs.probs
                    if sum(subs_probs) == 0.0:
                        logger.warning('Probability of substitutions equal to 0 in {}'.format(self.name))
                        self.result['partitions'] = []
                        self.result['sampling_size'] = self.sampling_size
                        self.result['obs'] = None
                        self.result['neg_obs'] = None
                        return self
                    subs_probs = subs_probs * self.p_subs
                    subs_probs = list(subs_probs)
                else:  # Same prob for all values (according to the prob of substitutions)
                    subs_probs = [self.p_subs / len(subs_scores)] * len(subs_scores) if len(subs_scores) > 0 else []

            if self.use_indels and self.p_indels > 0:
                indels_scores = self.indels.get_background_indel_scores()

                # All indels have the same probability
                indels_probs = [self.p_indels/len(indels_scores)] * len(indels_scores) if len(indels_scores) > 0 else []

            simulation_scores = subs_scores + indels_scores
            simulation_probs = subs_probs + indels_probs

            muts_count = len(self.result['mutations'])
            chunk_count = (self.sampling_size * muts_count) // self.sampling_chunk
            chunk_size = self.sampling_size if chunk_count == 0 else self.sampling_size // chunk_count

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
            if len(self.result['partitions']) > 0 or obs < self.min_obs:
                self.result['muts_count'] = muts_count
                self.result['simulation_scores'] = simulation_scores
                self.result['simulation_probs'] = simulation_probs
                self.result['observed'] = observed
                self.result['statistic_name'] = self.statistic_name

        self.result['sampling_size'] = self.sampling_size
        self.result['symbol'] = self.symbol

        return self
