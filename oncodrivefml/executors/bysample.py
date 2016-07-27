from collections import defaultdict

import numpy as np

from oncodrivefml.executors.bymutation import ElementExecutor, detect_repeatitive_seq
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS


class GroupBySampleExecutor(ElementExecutor):
    """
    The scores within an element are grouped by the sample_id.
    Within all the mutations with the same sample_id, the max is taken
    and used for the scores list.

    Args:
        name:
        muts:
        segments:
        signature:
        config:
    """

    def __init__(self, name, muts, segments, signature, config):

        # Input attributes
        self.name = name
        self.signature = signature
        self.segments = segments
        subs = config['statistic']['subs'] != 'none'
        if subs:
            self.muts = [m for m in muts if m['TYPE'] == 'subs']
        else:
            self.muts = []
        self.indels = config['statistic']['indels'] != 'none'
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

        self.is_positive_strand = True if segments[0].get('strand', '+') == '+' else False

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

        # Load element scores
        self.scores = Scores(self.name, self.segments, self.signature, self.score_config)

        # Compute observed mutations statistics and scores
        self.result = self.compute_muts_statistics(self.muts, self.scores, indels=self.indels, positive_strand=self.is_positive_strand)

        min_background_size = 0

        if len(self.result['mutations']) > 0:
            statistic_test = STATISTIC_TESTS.get(self.statistic_name)
            observed = []
            background = []

            scores_by_signature = defaultdict(list)
            for mut in self.result['mutations']:
                scores_by_signature[mut[self.signature_column]].append(mut['SCORE'])

            for sample, scores in scores_by_signature.items():

                simulation_scores = []
                simulation_signature = []

                positions = self.scores.get_all_positions()

                signature = self.signature

                for pos in positions:
                    for s in self.scores.get_score_by_position(pos):
                        simulation_scores.append(s.value)
                        simulation_signature.append(self.signature[sample].get((s.ref_triplet, s.alt_triplet), 0.0))

                simulation_scores = np.array(simulation_scores)
                simulation_signature = np.array(simulation_signature)
                simulation_signature = simulation_signature / simulation_signature.sum()

                observed.append(statistic_test.calc(scores))
                random_scores = np.random.choice(simulation_scores, size=(self.sampling_size, len(scores)), p=simulation_signature, replace=True)
                # matrix sampling_size, len(scores)
                background.append(np.max(random_scores, axis=1))  # gets the max of each row

                if min_background_size == 0 or min_background_size > len(simulation_scores):
                    min_background_size = len(simulation_scores)

            self.obs, self.neg_obs = statistic_test.calc_observed(np.array(background).transpose(), np.array(observed))

        # Calculate p-values
        self.result['min_background_size'] = min_background_size
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
