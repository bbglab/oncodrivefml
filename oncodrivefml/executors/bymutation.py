from collections import defaultdict

from oncodrivefml.executors.element import ElementExecutor


class GroupByMutationExecutor(ElementExecutor):
    """
    Executor that simulates the same number of mutations as the ones
    observed in the element.
    The simulation parameters are taken from the configuration file.

    Args:
        element_id (str): element ID
        muts (list): list of mutations belonging to that element (see :ref:`mutations <mutations dict>` inner items)
        segments (list): list of segments belonging to that element (see :ref:`elements <elements dict>` inner items)
        signature (dict): probabilities of each mutation (see :ref:`signatures <signature dict>`)
        config (dict): configurations

    """
    def __init__(self, element_id, muts, segments, signature, config):
        super(GroupByMutationExecutor, self).__init__(element_id, muts, segments, signature, config)

    def compute_muts_statistics(self, indels=False):
        """
        Gets the score of each mutation

        Args:
            indels (:obj:`~oncodrivefml.indels.Indel`): Indels class if indels are considered. False otherwise.

        Returns:
            dict: several information about the mutations and a list of them with the scores

        """

        # Add scores to the element mutations
        scores_by_sample = defaultdict(list)
        scores_list = []
        positions = []
        mutations = []
        amount_of_subs = 0
        amount_of_indels = 0

        for m in self.muts:

            self.compute_mutation_score(m, indels)

            # Update scores
            if m.get('SCORE', None) is not None:

                sample = m['SAMPLE']

                scores_by_sample[sample].append(m['SCORE'])

                scores_list.append(m['SCORE'])

                if m['TYPE'] == "subs" or m['TYPE'] == "mnp":
                    amount_of_subs += 1
                elif m['TYPE'] == "indel":
                    amount_of_indels += 1

                positions.append(m['POSITION'])
                mutations.append(m)


        # Aggregate scores
        num_samples = len(scores_by_sample)

        item = {
            'samples_mut': num_samples,
            'muts': len(scores_list),
            'muts_recurrence': len(set(positions)),
            'subs': amount_of_subs,
            'indels': amount_of_indels,
            'scores': scores_list,
            'positions': positions,
            'mutations': mutations
        }

        return item
