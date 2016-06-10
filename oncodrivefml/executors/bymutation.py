import logging

import numpy as np
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS
from oncodrivefml.signature import get_ref

from oncodrivefml.executors.indels import Indel
import math

class ElementExecutor(object):

    @staticmethod
    def compute_muts_statistics(muts, scores, indels=False):
        """
        For each mutation in muts, get the score for that position and that
        mutation

        :param muts: list of mutations corresponding to the element in the
        variants_dict
        :param scores: { pos : [ ( ref, alt, scoreValue,  {signatures:prob} ) ] }
        :param indels:
        :return:
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

                score = Indel.get_indel_score(m, scores, int(m['POSITION']))

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
        Computes the element p-value
        """
        raise RuntimeError("The classes that extend ElementExecutor must override the run() method")


def detect_repeatitive_seq(chrom, seq, pos):
    size = len(seq)
    repeats = 0
    while seq == get_ref(chrom, pos, size=size):
        pos += size
        repeats += 1
    return repeats


class GroupByMutationExecutor(ElementExecutor):
    """
    This executor simulates each mutation independently following the signatures probability
    and within a range if it's provided.
    """

    def __init__(self, element_id, muts, segments, signature, config):
        """

        :param element_id: element id
        :param muts: list of mutations corresponding to the element in the
        variants_dict
        :param segments: list of values from the elements dict
        :param signature: dict with all signatures
        :param config: dict with the configuration
        """

        # Input attributes
        self.name = element_id
        self.indels = config['statistic']['indels'] != 'none'
        self.muts = [m for m in muts if m['TYPE'] == 'subs']

        # Add only indels that are not in a repeatitive sequence
        if self.indels:
            indels_max_repeats = config['statistic']['indels_max_repeats']
            indels_set = set()
            for m in [m for m in muts if m['TYPE'] == 'indel']:
                chrom = m['CHROMOSOME']
                ref = m['REF']
                alt = m['ALT']
                pos = m['POSITION']

                # Skip recurrent indels
                if (chrom, pos, ref, alt) not in indels_set:
                    indels_set.add((chrom, pos, ref, alt))

                    # Check if it's repeated
                    seq = alt if '-' in ref else ref
                    repeats = detect_repeatitive_seq(chrom, seq, pos)
                    if repeats <= indels_max_repeats:
                        self.muts.append(m)

        self.signature = signature
        self.segments = segments

        # Configuration parameters
        self.score_config = config['score']
        self.sampling_size = config['background'].get('sampling', 100000)
        self.statistic_name = config['statistic'].get('method', 'amean')
        self.simulation_range = config['background'].get('range', None)
        self.signature_column = 'SAMPLE' if config['signature'].get('method', 'full') == 'bysample' else 'SIGNATURE'

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
        self.result = self.compute_muts_statistics(self.muts, self.scores, indels=self.indels)

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
                            simulation_signature.append(s.signature.get(mut[self.signature_column]))
                else: #indels
                    is_insertion = True if '-' in mut['REF'] else False
                    signature = None
                    if is_insertion:
                        for pos in positions:
                            if get_ref(pos) == get_ref(mut['POSITION']) and get_ref(pos + 1) == get_ref(mut['POSITION'] + 1):
                                score = Indel.get_indel_score(mut, self.scores, pos)
                                if not math.isnan(score):
                                    simulation_scores.append(score)

                    else:
                        del_size = max(len(mut['REF']), len(mut['ALT']))
                        for pos in positions:
                            if get_ref(pos) == get_ref(mut['POSITION']) and get_ref(pos + del_size) == get_ref(mut['POSITION'] + del_size):
                                score = Indel.get_indel_score(mut, self.scores, pos)
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

            self.obs, self.neg_obs = statistic_test.calc_observed(np.array(background).transpose(), np.array(observed))

        # Calculate p-values
        self.result['background_size'] = len(simulation_scores)
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
