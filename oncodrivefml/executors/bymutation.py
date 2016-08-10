import logging

import numpy as np
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS
from oncodrivefml.signature import get_ref

from oncodrivefml.executors.indels import Indel, complements_dict
import math

class ElementExecutor(object):

    @staticmethod
    def compute_muts_statistics(muts, scores, indels=False, positive_strand=True):
        """
        Gets the score of each mutation

        Args:
            muts (list): list of mutations
            scores (dict): scores for all possible substitutions
            indels (bool): are indels taken into account or not. Default to False.
            positive_strand (bool): the element where the mutations occur has positive strand or not. Defaults to True.

        Returns:
            dict: several information about the mutations and a list of them with the scores

        """

        # Add scores to the element mutations
        scores_by_sample = {}
        scores_list = []
        scores_subs_list = []
        scores_indels_list = []
        subs_counter = 0 #counts how many are type substition
        total_subs_score = 0 #counts how many substitutons have a score value
        positions = []
        mutations = []
        for m in muts:

            # Get substitutions scores
            if m['TYPE'] == "subs":
                subs_counter += 1
                m['POSITION'] = int(m['POSITION'])
                values = scores.get_score_by_position(m['POSITION'])
                for v in values:
                    if v.ref == m['REF'] and v.alt == m['ALT']:
                        m['SCORE'] = v.value
                        total_subs_score += 1
                        break

            if indels and m['TYPE'] == "indel":

                score = Indel.get_indel_score(m, scores, int(m['POSITION']), positive_strand)

                m['SCORE'] = score if not math.isnan(score) else None

            # Update scores
            if m.get('SCORE', None) is not None:

                sample = m['SAMPLE']
                if sample not in scores_by_sample:
                    scores_by_sample[sample] = []

                scores_by_sample[sample].append(m['SCORE'])
                scores_list.append(m['SCORE'])

                if m['TYPE'] == "subs":
                    scores_subs_list.append(m['SCORE'])
                elif m['TYPE'] == "indel":
                    scores_indels_list.append(m['SCORE'])

                positions.append(m['POSITION'])
                mutations.append(m)

        # Aggregate scores
        num_samples = len(scores_by_sample)

        item = {
            'samples_mut': num_samples,
            'muts': len(scores_list),
            'muts_recurrence': len(set(positions)),
            'subs': subs_counter,
            'subs_score': total_subs_score,
            'scores': scores_list,
            'scores_subs': scores_subs_list,
            'scores_indels': scores_indels_list,
            'positions': positions,
            'mutations': mutations
        }

        return item

    def run(self):
        """

        Raises:
            RuntimeError: method must be implemented

        """
        raise RuntimeError("The classes that extend ElementExecutor must override the run() method")


def detect_repeatitive_seq(chrom, seq, pos):
    """
    How many times a sequences is repeated (in the reference genome)
    starting in a certain position

    Args:
        chrom (str): chromosome ID
        seq (str): sequence of bases to be detected
        pos (int): position  where to start looking

    Returns:
        int: how many consequitive times the sequence appears

    """
    size = len(seq)
    repeats = 0
    while seq == get_ref(chrom, pos, size=size):
        pos += size
        repeats += 1
    return repeats


class GroupByMutationExecutor(ElementExecutor):
    """
    Excutor to simulate each mutation in a set. The simulation parameters are taken from the
    configuration file.

    Args:
        element_id (str): element ID
        muts (list): list of mutations belonging to that element (see :ref:`mutations <mutations dict>` inner items)
        segments (list): list of segments belonging to that element (see :ref:`elements <elements dict>` inner items)
        signature (dict): probabilites of each mutation (see :ref:`signatures <signature dict>`)
        config (dict): configurations

    Indels where the sequence is repeated a predefined number of times are
    discarded.
    """
    def __init__(self, element_id, muts, segments, signature, config):
        # Input attributes
        self.name = element_id
        self.indels = config['statistic']['indels'].get('enabled', False)
        subs = config['statistic'].get('subs', False)

        if subs and self.indels:
            self.muts = muts
        else:
            if subs:
                self.muts = [m for m in muts if m['TYPE'] == 'subs']
            else:
                self.muts =[m for m in muts if m['TYPE'] == 'indel']

        self.signature = signature
        self.segments = segments
        self.is_positive_strand = True if segments[0].get('strand', '+') == '+' else False

        # Configuration parameters
        self.score_config = config['score']
        self.sampling_size = config['background'].get('sampling', 100000)
        self.statistic_name = config['statistic'].get('method', 'amean')
        self.simulation_range = config['background'].get('range', None)
        self.signature_column = config['signature']['classifier']

        # Output attributes
        self.obs = 0
        self.neg_obs = 0
        self.result = None
        self.scores = None

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
        self.scores = Scores(self.name, self.segments, self.signature, self.score_config)

        # Compute observed mutations statistics and scores
        self.result = self.compute_muts_statistics(self.muts, self.scores, indels=self.indels, positive_strand=self.is_positive_strand)

        min_background_size = 0

        if len(self.result['mutations']) > 0:
            statistic_test = STATISTIC_TESTS.get(self.statistic_name)
            observed = []
            background = []

            for mut in self.result['mutations']:

                simulation_scores = []
                simulation_signature = []

                if self.simulation_range is not None:
                    positions = range(mut['POSITION'] - self.simulation_range, mut['POSITION'] + self.simulation_range)
                else:
                    positions = self.scores.get_all_positions()

                signature = self.signature

                if mut['TYPE'] == 'subs':
                    for pos in positions:
                        for s in self.scores.get_score_by_position(pos):
                            simulation_scores.append(s.value)
                            #TODO KeyError
                            if signature is not None:
                                simulation_signature.append(self.signature[mut.get(self.signature_column, self.signature_column)].get((s.ref_triplet, s.alt_triplet), 0.0))

                else: #indels
                    #TODO change the method to compute first the position and then the scores only for those
                    indel_size = max(len(mut['REF']), len(mut['ALT']))
                    length = Indel.compute_window_size(indel_size)

                    mutation_pattern = Indel.get_pattern(mut, self.is_positive_strand, length)
                    signature = None

                    sampling_positions = positions

                    for pos in sampling_positions:
                        score = Indel.get_indel_score_for_background(self.scores, pos, length, mut['CHROMOSOME'],
                                                                     mutation_pattern, indel_size, self.is_positive_strand)
                        if not math.isnan(score):
                            simulation_scores.append(score)

                    if len(simulation_scores) < 100:
                        logging.warning("Element {} and mutation {} has only {} valid background scores".format(self.name, mut, len(simulation_scores)))

                simulation_scores = np.array(simulation_scores)

                if signature is not None:
                    simulation_signature = np.array(simulation_signature)
                    simulation_signature = simulation_signature / simulation_signature.sum()
                else:
                    simulation_signature = None

                observed.append(mut['SCORE'])
                background.append(np.random.choice(simulation_scores, size=self.sampling_size, p=simulation_signature, replace=True))

                if min_background_size == 0 or min_background_size > len(simulation_scores):
                    min_background_size = len(simulation_scores)

            self.obs, self.neg_obs = statistic_test.calc_observed(np.array(background).transpose(), np.array(observed))

        # Calculate p-values
        self.result['min_background_size'] = min_background_size
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
