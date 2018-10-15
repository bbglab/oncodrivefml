from collections import defaultdict

from oncodrivefml.executors.element import ElementExecutor


class GroupByMutationExecutor(ElementExecutor):
    """
    Executor that simulates the same number of mutations as the ones
    observed in the element.
    The simulation parameters are taken from the configuration file.

    """

    def compute_muts_statistics(self):
        """
        Gets the score of each mutation

        Returns:
            dict: several information about the mutations and a list of them with the scores

        """

        # Add scores to the element mutations
        scores_by_sample = defaultdict(list)
        scores_list = []
        positions = []
        mutations = []
        amount_of_snps = 0
        amount_of_mnps = 0
        amount_of_indels = 0

        for m in self.muts:

            self.compute_mutation_score(m)

            # Update scores
            if m.get('SCORE', None) is not None:

                sample = m['SAMPLE']

                scores_by_sample[sample].append(m['SCORE'])

                scores_list.append(m['SCORE'])

                if m['ALT_TYPE'] == "snp":
                    amount_of_snps += 1
                elif m['ALT_TYPE'] == "mnp":
                    amount_of_mnps +=1
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
            'mutations': mutations
        }

        return item
