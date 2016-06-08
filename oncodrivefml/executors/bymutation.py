import numpy as np
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS
from oncodrivefml.signature import get_ref

import math

class Weight:
    '''
    In 0 it will always be 1
    In length + 1 is 0
    '''
    def __init__(self, length, type):
        self.length = length
        if type == 'cte':
            self.funct = self.cte
        if type=='linear':
            self.intercept = 1
            self.slope = -1.0/length
            self.funct = self.linear

    def function(self, x):
        return self.funct(x)

    def cte(self, x):
        return 1

    def linear(self, x):
        return self.slope * x + self.intercept

def weigth_indel_scores(indel_scores, indel_size, indel_window_size, funct):
    # if we have and indel multiple of 3 all have weight 1
    if not(indel_size % 3):
        return indel_scores
    #During the indel size, weight is one, so we ignore those
    for i in range(indel_size, indel_window_size):#if indel_size is bigger than the window it is an empty range
        # We can pass to the weight function the value i or the value i-indel_size+1. The first means using the value of the function as it starts in 0, the second as it start when the indel ends
        indel_scores[i] = indel_scores[i] * funct(i-indel_size+1) # first element has x = 1 TODO check



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
                m['POSITION'] = int(m['POSITION'])
                indel_size = max(len(m['REF']), len(m['ALT']))
                is_insertion = True if '-' in m['REF'] else False

                indel_window_size = 10
                indel_scores = [math.nan]* indel_window_size

                weight = Weight(indel_window_size, 'cte')
                # TODO move this outside
                if not (indel_size % 3): # frame shift
                    indel_window_size = (indel_size + 2) if (indel_size+2) <= indel_window_size else indel_window_size
                    #TODO adapt window size for triplets

                index = 0
                for position in range(m['POSITION'], m['POSITION'] + indel_window_size):
                    scores_in_position = scores.get_score_by_position(position) # if position is outside the element, it has not score so the indel score remains nan
                    for score in scores_in_position:
                        if is_insertion:
                            alteration = m['ALT'][index] if index < indel_size else get_ref(position-indel_size)
                        else:  # deletion
                            try:
                                alteration = get_ref(position+indel_size)
                            except ValueError:  # generated when we are outside the element
                                break # ignore it (score = nan)
                        #alteration = (m['ALT'] if is_insertion else m['REF'])[index] if index < indel_size else get_ref(position-indel_size)
                        if score.ref == get_ref(position) and score.alt == alteration:
                            indel_scores[index] = score.value
                            break
                    index += 1

                #indel_window_size = index #adjust the window size to the looped size. Only useful if we are at the end of the element

                indel_scores = weigth_indel_scores(indel_scores, indel_size, indel_window_size, weight.function)

                m['SCORE'] = max(indel_scores)

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

                for pos in positions:
                    for s in self.scores.get_score_by_position(pos):
                        simulation_scores.append(s.value)
                        simulation_signature.append(s.signature.get(mut[self.signature_column]))

                simulation_scores = np.array(simulation_scores)

                if self.signature is not None:
                    simulation_signature = np.array(simulation_signature)
                    simulation_signature = simulation_signature / simulation_signature.sum()
                else:
                    simulation_signature = None

                observed.append(mut['SCORE'])
                background.append(np.random.choice(simulation_scores, size=self.sampling_size, p=simulation_signature, replace=True))

            self.obs, self.neg_obs = statistic_test.calc_observed(zip(*background), observed)

        # Calculate p-values
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
