import logging

import numpy as np

from oncodrivefml.executors.bymutation import GroupByMutationExecutor
from oncodrivefml.scores import Scores

from oncodrivefml.executors.indels import Indel
import math


class GroupByMutationWeigthedExecutor(GroupByMutationExecutor):

    def __init__(self, element_id, muts, segments, signature, muts_by_sample, config):
        super().__init__(element_id, muts, segments, signature, config)
        self.muts_by_sample = muts_by_sample

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

            observed = []
            background = []

            weights = []
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
                weights.append(1 / self.muts_by_sample[mut['SAMPLE']])
                background.append(np.random.choice(simulation_scores, size=self.sampling_size, p=simulation_signature, replace=True))

                if min_background_size == 0 or min_background_size > len(simulation_scores):
                    min_background_size = len(simulation_scores)

            # Calculate statistic
            values = np.array(background).transpose()
            observed = np.array(observed)

            # Weighted
            observed_value = np.average(observed, weights=weights)
            values = np.averages(values, axis=1, weights=weights)

            self.obs = len(values[values >= observed_value])
            self.neg_obs = len(values[values <= observed_value])


        # Calculate p-values
        self.result['min_background_size'] = min_background_size
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
