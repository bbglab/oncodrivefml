from collections import defaultdict

from oncodrivefml.executors.bymutation import ElementExecutor
from oncodrivefml.stats import STATISTIC_TESTS


class GroupBySampleExecutor(ElementExecutor):
    """
    Executor that simulates one mutation per sample (if the sample exists in the observed mutations).
    The simulation parameters are taken from the configuration file.

    """

    def compute_muts_statistics(self):
        """
        Gets the score of each mutation.
        The mutation (per sample) with the highest score is used as a reference for the simulation.
        It means that this mutation is used to get the signature, position...

        Returns:
            dict: several information about the mutations and a list of them with the scores

        """

        samples_statistic_test = STATISTIC_TESTS.get(self.samples_method)

        # Add scores to the element mutations
        scores_by_sample = defaultdict(list)
        scores_list = []
        positions = []
        mutations = []
        amount_of_snps = 0
        amount_of_mnps = 0
        amount_of_indels = 0

        mut_per_sample = {}

        for m in self.muts:

            self.compute_mutation_score(m)

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

            if m['ALT_TYPE'] == "snp" :
                amount_of_snps += 1
            elif m['ALT_TYPE'] == "mnp":
                amount_of_mnps += 1
            elif m['ALT_TYPE'] == "indel":
                amount_of_indels += 1

            positions.append(m['POSITION'])
            mutations.append(m)

        # Aggregate scores
        num_samples = len(scores_by_sample)

        item = {
            'samples_mut': num_samples,
            'muts': len(scores_list),
            'muts_recurrence': len(set(positions)),
            'snps': amount_of_snps,
            'mnps': amount_of_mnps,
            'indels': amount_of_indels,
            'scores': scores_list,
            'positions': positions,
            'mutations': mutations,
            'scores_by_sample': scores_by_sample
        }

        return item
