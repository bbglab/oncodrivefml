import numpy as np
from oncodrivefml.scores import Scores
from oncodrivefml.statistics import STATISTIC_TESTS


class ElementExecutor(object):

    def __init__(self, name, muts, segments, signature, config):

        # Input attributes
        self.name = name
        self.muts = muts
        self.signature = signature
        self.segments = segments
        self.config = config

        # Output attributes
        self.result = None
        self.scores = None

    def run(self):

        self.scores = Scores(self.name, self.segments, self.signature, self.config['score'])

        # Compute elements statistics
        self.result = self.compute_muts_statistics()

        # Run only if there are valid mutations
        if self.result['muts'] > 0:

            sampling_size = self.config['background'].get('sampling', 100000)
            obs, neg_obs = self.sampling(sampling_size)

            self.result['pvalue'] = max(1, obs) / float(sampling_size) if obs is not None else None
            self.result['pvalue_neg'] = max(1, neg_obs) / float(sampling_size) if neg_obs is not None else None

        return self

    def sampling(self, sampling_size):

        statistic_test = STATISTIC_TESTS.get(self.config['statistic'].get('method', 'amean'), STATISTIC_TESTS['amean'])
        simulation_range = self.config['background'].get('range', 1000)

        observed = []
        background = []
        for mut in self.result['mutations']:
            if mut['TYPE'] != 'subs':
                continue

            simulation_scores = []
            simulation_signature = []
            for pos in range(mut['POSITION'] - simulation_range, mut['POSITION'] + simulation_range):
                for s in self.scores.get_score_by_position(pos):
                    simulation_scores.append(s.value)
                    simulation_signature.append(s.signature.get(mut['SIGNATURE']))

            simulation_scores = np.array(simulation_scores)
            simulation_signature = np.array(simulation_signature)
            simulation_signature = simulation_signature / simulation_signature.sum()

            observed.append(mut['SCORE'])
            background.append(np.random.choice(simulation_scores, size=sampling_size, p=simulation_signature, replace=True))

        observed_test = statistic_test.calc(observed)
        obs, neg_obs = 0.0, 0.0
        for values in zip(*background):
            background_test = statistic_test.calc(values)
            if background_test >= observed_test:
                obs += 1
            if background_test <= observed_test:
                neg_obs += 1

        return obs, neg_obs

    def compute_muts_statistics(self):

        # Skip features without mutations
        if len(self.muts) == 0:
            return None

        # Add scores to the element mutations
        scores_by_sample = {}
        scores_list = []
        scores_subs_list = []
        scores_indels_list = []
        total_subs = 0
        total_subs_score = 0
        positions = []
        mutations = []
        for m in self.muts:

            # Get substitutions scores
            if m['TYPE'] == "subs":
                total_subs += 1
                m['POSITION'] = int(m['POSITION'])
                values = self.scores.get_score_by_position(m['POSITION'])
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
            'positions': positions,
            'mutations': mutations
        }

        return item
