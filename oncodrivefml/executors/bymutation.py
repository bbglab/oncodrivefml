import numpy as np
from oncodrivefml.scores import Scores
from oncodrivefml.stats import STATISTIC_TESTS


class ElementExecutor(object):

    @staticmethod
    def compute_muts_statistics(muts, scores):

        # Add scores to the element mutations
        scores_by_sample = {}
        scores_list = []
        scores_subs_list = []
        scores_indels_list = []
        total_subs = 0
        total_subs_score = 0
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
                        m['SCORE'] = v.value
                        total_subs_score += 1
                        break

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
            'scores_by_sample': scores_by_sample,
            'positions': positions,
            'mutations': mutations
        }

        return item

    def run(self):
        """
        Computes the element p-value
        """
        raise RuntimeError("The classes that extend ElementExecutor must override the run() method")


class GroupByMutationExecutor(ElementExecutor):
    """
    This executor simulates each mutation independently following the signature probability
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
                simulation_signature = np.array(simulation_signature)
                simulation_signature = simulation_signature / simulation_signature.sum()

                observed.append(mut['SCORE'])
                background.append(np.random.choice(simulation_scores, size=self.sampling_size, p=simulation_signature, replace=True))

            self.obs, self.neg_obs = statistic_test.calc_observed(zip(*background), observed)

        # Calculate p-values
        self.result['pvalue'] = max(1, self.obs) / float(self.sampling_size)
        self.result['pvalue_neg'] = max(1, self.neg_obs) / float(self.sampling_size)

        return self
