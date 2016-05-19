import numpy as np
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS
from oncodrivefml.signature import get_ref


class ElementExecutor(object):

    @staticmethod
    def compute_muts_statistics(muts, scores, indels=False):
        """
        For each mutation in muts, get the score for that position and that
        mutation

        :param muts: list of mutations corresponding to the element in the
        variants_dict
        :param scores: { pos : [ ( ref, alt, scoreValue,  {signature:prob} ) ] }
        :param indels:
        :return:
        """

        # Add scores to the element mutations
        scores_by_sample = {}
        scores_list = []
        scores_subs_list = []
        scores_indels_list = []
        total_subs = 0 #counts how many are type substition
        total_subs_score = 0 #counts how many substitutons have a score value
        positions = []
        mutations = []
        for m in muts:

            # Get substitutions scores
            if m['TYPE'] == "subs":
                total_subs += 1
                m['POSITION'] = int(m['POSITION'])
                values = scores.get_score_by_position(m['POSITION'])
                for v in values:
                    if v.ref == m['REF'] and v.alt == m['ALT']:
                        #ASK why signature independent
                        m['SCORE'] = v.value
                        total_subs_score += 1
                        break

            if indels and m['TYPE'] == "indel":
                m['POSITION'] = int(m['POSITION'])
                indel_size = max(len(m['REF']), len(m['ALT']))

                indel_scores = []
                for pos in range(m['POSITION'], m['POSITION'] + indel_size):
                    indel_scores += [s.value for s in scores.get_score_by_position(pos)]

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
            'subs': total_subs,
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
    This executor simulates each mutation independently following the signature probability
    and within a range if it's provided.
    """

    def __init__(self, name, muts, segments, signature, config):
        """

        :param name: element name
        :param muts: list of mutations corresponding to the element in the
        variants_dict
        :param segments: list of values from the elements dict
        :param signature: dict with all signatures
        :param config: dict with the configuration
        """

        # Input attributes
        self.name = name
        self.indels = config['statistic']['indels'] != 'none'
        self.muts = [m for m in muts if m['TYPE'] == 'subs']
        #ASK won't this avoid getting indels below??

        # Add only indels if there is at least one substitution
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
        For all positions around the mutation position


        :return: GroupByMutationExecutor
        """

        # Load element scores
        self.scores = Scores(self.name, self.segments, self.signature, self.score_config)

        # Compute observed mutations statistics and scores
        self.result = self.compute_muts_statistics(self.muts, self.scores, indels=self.indels)

        if len(self.result['mutations']) > 0:
            statistic_test = STATISTIC_TESTS.get(self.statistic_name)
            observed = []
            background = []

            for mut in self.result['mutations']:#ASK acts only as a counter
                # if no simulation rate is set?

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
                        #ASK what if the signatrue of the mut is different
                        # than the signature of in ts

                simulation_scores = np.array(simulation_scores)
                simulation_signature = np.array(simulation_signature)
                simulation_signature = simulation_signature / simulation_signature.sum()

                observed.append(mut['SCORE'])
                background.append(np.random.choice(simulation_scores, size=self.sampling_size, p=simulation_signature, replace=True))

            self.obs, self.neg_obs = statistic_test.calc_observed(zip(*background), observed)

        # Calculate p-values
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
