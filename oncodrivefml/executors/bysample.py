from collections import defaultdict, Counter

import logging
import numpy as np

from oncodrivefml.executors.bymutation import ElementExecutor
from oncodrivefml.indels import Indel
from oncodrivefml.scores import Scores
from oncodrivefml.signature import get_ref
from oncodrivefml.stats import STATISTIC_TESTS



class GroupBySampleExecutor(ElementExecutor):
    """
    Executor that simulates one mutation per sample (if the sample exists in the observed mutations).
    The simulation parameters are taken from the configuration file.

    Args:
        element_id (str): element ID
        muts (list): list of mutations belonging to that element (see :ref:`mutations <mutations dict>` inner items)
        segments (list): list of segments belonging to that element (see :ref:`elements <elements dict>` inner items)
        signature (dict): probabilities of each mutation (see :ref:`signatures <signature dict>`)
        config (dict): configurations

    """

    def __init__(self, element_id, muts, segments, signature, config):
        super(GroupBySampleExecutor, self).__init__(element_id, muts, segments, signature, config)

    def compute_muts_statistics(self, indels=False):
        """
        Gets the score of each mutation.
        The mutation (per sample) with the highest score is used as a reference for the simulation.
        It means that this mutation is used to get the signature, position...

        Args:
            indels (:obj:`~oncodrivefml.indels.Indel`): Indels class if indels are considered. False otherwise

        Returns:
            dict: several information about the mutations and a list of them with the scores

        """

        samples_statistic_test = STATISTIC_TESTS.get(self.samples_method)

        # Add scores to the element mutations
        scores_by_sample = defaultdict(list)
        scores_list = []
        positions = []
        mutations = []
        amount_of_subs = 0
        amount_of_indels = 0

        mut_per_sample = {}

        for m in self.muts:

            self.compute_mutation_score(m, indels)

            # Update scores
            if m.get('SCORE', None) is not None:

                sample = m['SAMPLE']

                if sample not in mut_per_sample.keys() or m['SCORE'] > mut_per_sample[sample]['SCORE']:
                    mut_per_sample[sample] = m

                scores_by_sample[sample].append(m['SCORE'])

        for sample, m in mut_per_sample.items():

            m['SCORE'] = samples_statistic_test.calc(scores_by_sample[sample])
            # Create an imaginary mutation with the same parameters as the one with the maximum score
            # and the score result from the statistical test on all the scores for all mutation in that sample

            scores_list.append(m['SCORE'])

            if m['TYPE'] == "subs" or m['TYPE'] == "mnp":
                amount_of_subs += 1
            elif m['TYPE'] == "indel":
                amount_of_indels += 1

            positions.append(m['POSITION'])
            mutations.append(m)

        # Aggregate scores
        num_samples = len(scores_by_sample)

        item = {
            'samples_mut': num_samples,
            'muts': len(scores_list),
            'muts_recurrence': len(set(positions)),
            'subs': amount_of_subs,
            'indels': amount_of_indels,
            'scores': scores_list,
            'positions': positions,
            'mutations': mutations,
            'scores_by_sample': scores_by_sample
        }

        return item
