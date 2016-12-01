import logging
from collections import Counter

import numpy as np
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS
from oncodrivefml.signature import get_ref

from oncodrivefml.indels import Indel
import math


class ElementExecutor(object):

    @staticmethod
    def compute_muts_statistics(muts, scores, indels=False):
        """
        Gets the score of each mutation

        Args:
            muts (list): list of mutations
            scores (dict): scores for all possible substitutions
            indels (:obj:`~oncodrivefml.indels.Indel`): Indels class if indels are considered. False otherwise.
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

            if m['TYPE'] == "MNP":
                subs_counter += 1  #It is counted as a single substitution
                pos = int(m['POSITION'])
                mnp_scores = []
                for index, nucleotides in enumerate(zip(m['REF'], m['ALT'])):
                    ref_nucleotide, alt_nucleotide = nucleotides
                    values = scores.get_score_by_position(pos+index)
                    for v in values:
                        if v.ref == ref_nucleotide and v.alt == alt_nucleotide:
                            mnp_scores.append(v.value)
                            break
                    else:
                        logging.warning('Discrepancy in MNP at position {} of chr {}'.format(pos, m['CHROMOSOME']))
                if not mnp_scores:
                    continue
                m['SCORE'] = max(mnp_scores)
                m['POSITION'] = pos + mnp_scores.index(m['SCORE'])

            if indels is not False and m['TYPE'] == "indel":

                # very long indels are discarded
                if max(len(m['REF']), len(m['ALT'])) > 20:
                    continue

                score = indels.get_indel_score(m)

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
        self.use_indels = config['statistic']['indels'].get('enabled', False)
        self.indels_conf = config['statistic']['indels']
        self.use_subs = config['statistic'].get('subs', True)
        self.use_mnp = config['statistic'].get('mnp', True)

        if self.use_subs and self.use_indels and self.use_mnp:
            self.muts = muts
        else:
            self.muts = []
            if self.use_subs:
                self.muts += [m for m in muts if m['TYPE'] == 'subs']
            if self.use_mnp:
                self.muts += [m for m in muts if m['TYPE'] == 'MNP']
            if self.use_indels:
                self.muts += [m for m in muts if m['TYPE'] == 'indel']

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

        self.p_subs = config['p_subs']
        self.p_indels = config['p_indels']

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
        self.result = self.compute_muts_statistics(self.muts, self.scores, indels=indels)

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


            for mut in self.result['mutations']:
                observed.append(mut['SCORE']) # Observed mutations
                if mut['TYPE'] == 'subs' or mut['TYPE']=='MNP':
                    if self.signature is not None:
                        # Count how many signature ids are and prepare a vector for each
                        # IMPORTANT: this implies that only the signature of the observed mutations is taken into account
                        signature_id = mut.get(self.signature_column, self.signature_column)
                        if signature_id not in signature_ids:
                            subs_probs_by_signature.update({signature_id: []})
                        signature_ids.append(signature_id)

            # When the probabilities of subs and indels are None, they are taken from
            # mutations seen in the gene
            if self.p_subs is None or self.p_indels is None: # use the probabilities based on observed mutations
                self.p_subs = self.result['subs'] / len(self.result['mutations'])
                self.p_indels = 1 - self.p_subs


            # Compute the values for the substitutions
            if self.use_subs and self.p_subs > 0:
                for pos in positions:
                    for s in self.scores.get_score_by_position(pos):
                        subs_scores.append(s.value)
                        for k, v in subs_probs_by_signature.items():
                            v.append(self.signature[k].get((s.ref_triplet, s.alt_triplet), 0.0))

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


            background = np.random.choice(simulation_scores, size=(self.sampling_size, len(self.result['mutations'])), p=simulation_probs, replace=True)


            self.obs, self.neg_obs = statistic_test.calc_observed(background, np.array(observed))

        # Calculate p-values
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
