
import numpy as np

from oncodrivefml.executors.bymutation import ElementExecutor
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS


class GroupBySampleExecutor(ElementExecutor):
    """
    This executor simulates the observed mutations aggregated by sample.
    and within a range if it's provided.
    """

    def __init__(self, name, muts, segments, signature, config):

        # Input attributes
        self.name = name
        self.muts = [m for m in muts if m['TYPE'] == 'subs']
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

        # Load element scores
        self.scores = Scores(self.name, self.segments, self.signature, self.score_config)

        # Compute observed mutations statistics and scores
        self.result = self.compute_muts_statistics(self.muts, self.scores)

        if len(self.result['mutations']) > 0:
            statistic_test = STATISTIC_TESTS.get(self.statistic_name)
            observed = []
            background = []

            for sample, scores in self.result['scores_by_sample'].items():

                simulation_scores = []
                simulation_signature = []

                positions = self.scores.get_all_positions()

                for pos in positions:
                    for s in self.scores.get_score_by_position(pos):
                        simulation_scores.append(s.value)
                        simulation_signature.append(s.signature.get(sample))

                simulation_scores = np.array(simulation_scores)
                simulation_signature = np.array(simulation_signature)
                simulation_signature = simulation_signature / simulation_signature.sum()

                observed.append(statistic_test.calc(scores))
                random_scores = np.random.choice(simulation_scores, size=(self.sampling_size, len(scores)), p=simulation_signature, replace=True)
                background.append(np.max(random_scores, axis=1))

            self.obs, self.neg_obs = statistic_test.calc_observed(zip(*background), observed)

        # Calculate p-values
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
