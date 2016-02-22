from oncodrivefml.compute import load_scores, compute_muts_statistics, sampling


def run_executor(executor):
    return executor.run()


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

    def run(self):

        statistic_name = self.config['statistic'].get('method', 'amean')
        recurrence = self.config['background'].get('recurrence', True)

        # Load all element scores
        scores_by_position, scores_by_segment, signature_by_segment, missing_signatures = load_scores(
            self.name,
            self.segments,
            self.signature,
            self.config['score']
        )

        # Compute elements statistics
        self.result = compute_muts_statistics(self.muts, scores_by_position, statistic_name, recurrence)

        # Run only if there are valid mutations
        if self.result['muts'] > 0:

            sampling_size = self.config['background'].get('sampling', 100000)
            obs, neg_obs = sampling(sampling_size, scores_by_segment, signature_by_segment, self.name, self.result, statistic_name)

            self.result['pvalue'] = max(1, obs) / float(sampling_size) if obs is not None else None
            self.result['pvalue_neg'] = max(1, neg_obs) / float(sampling_size) if neg_obs is not None else None

        # Remove this key to simplify serialization
        del self.result['muts_by_tissue']

        return self
