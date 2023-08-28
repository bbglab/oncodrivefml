import math
import logging
import numpy as np

from oncodrivefml import __logger_name__
from oncodrivefml.executors import sig2probs
from oncodrivefml.indels import Indel
from oncodrivefml.scores import Scores
from oncodrivefml.mutability import Mutabilities
from oncodrivefml.stats import STATISTIC_TESTS
from oncodrivefml.walker import partitions_list

logger = logging.getLogger(__logger_name__)


class ElementExecutor(object):
    """
    Base class for executors that do the analysis per genomic element.
    The mutations are simulated to compute a p_value
    by comparing the functional impact of the observed mutations
    and the expected.
    The simulation parameters are taken from the configuration file.

    Args:
        element_id (str): element ID
        muts (list): list of mutations belonging to that element (see :ref:`mutations <mutations dict>` inner items)
        segments (list): list of segments belonging to that element (see :ref:`elements <elements dict>` inner items)
        signature (dict): probabilities of each mutation (see :ref:`signatures <signature dict>`)
        config (dict): configurations

    """

    def __init__(self, element_id, muts, segments, signature, config, seed):
        # Input attributes
        self.name = element_id
        self.seed = seed

        self.use_mnp = not config['statistic']['discard_mnp']
        self.indels_conf = config['statistic']['indels']
        self.use_indels = self.indels_conf['include']
        self.indels = None

        if self.use_indels and self.use_mnp:
            self.muts = muts
        else:
            self.muts = [m for m in muts if m['ALT_TYPE'] == 'snp']
            if self.use_mnp:
                self.muts += [m for m in muts if m['ALT_TYPE'] == 'mnp']
            if self.use_indels:
                self.muts += [m for m in muts if m['ALT_TYPE'] == 'indel']

        self.signature = signature
        self.segments = segments
        self.symbol = self.segments[0].get('SYMBOL', None)
        self.chromosome = self.segments[0].get('CHROMOSOME', None)
        # print(type(self.chromosome))

        # Configuration parameters
        self.score_config = config['score']
        self.sampling_size = config['statistic']['sampling']
        self.min_obs = config['statistic']['sampling_min_obs']
        self.sampling_chunk = config['statistic']['sampling_chunk'] * 10**6
        self.statistic_name = config['statistic']['method']
        self.signature_column = config['signature']['classifier']
        self.samples_method = config['statistic']['per_sample_analysis']
        
        self.mutability_info = config['mutability_info']
        self.depths_info = config['depths_info']


        # Output attributes
        self.obs = 0
        self.neg_obs = 0
        self.result = None
        self.scores = None

        self.p_subs = config['p_subs']
        self.p_indels = config['p_indels']

        if self.mutability_info:
            self.mutability_config = config['mutability']
        elif self.depths_info:
            self.depths = config['depths_loaded'][self.chromosome]

    def compute_muts_statistics(self):
        """
        Assigns scores to the mutations and retrieve the ones that are going
        to be simulated

        Returns:
            dict: several information about the mutations and a list of them with the scores

        """
        raise RuntimeError("The classes that extend ElementExecutor must override, at least the compute_muts_statistics() method")

    def compute_mutation_score(self, mutation):

        # Get substitutions scores
        if mutation['ALT_TYPE'] == "snp":
            mutation['POSITION'] = int(mutation['POSITION'])
            values = self.scores.get_score_by_position(mutation['POSITION'])
            for v in values:
                if v.ref == mutation['REF'] and v.alt == mutation['ALT']:
                    mutation['SCORE'] = v.value
                    break
            else:
                logger.warning('Reference mismatch or missing score in SNP at position {} of chr {}'.format(mutation['POSITION'], mutation['CHROMOSOME']))

        elif mutation['ALT_TYPE'] == "mnp":
            pos = int(mutation['POSITION'])
            mnp_scores = []
            for index, nucleotides in enumerate(zip(mutation['REF'], mutation['ALT'])):
                ref_nucleotide, alt_nucleotide = nucleotides
                values = self.scores.get_score_by_position(pos + index)
                for v in values:
                    if v.ref == ref_nucleotide and v.alt == alt_nucleotide:
                        mnp_scores.append(v.value)
                        break
                else:
                    # you should not worry about this warning,
                    # it is not a problem that some internal positions in the MNV do not have a score
                    # the nucleotide of the reference and the alternative might coincide in a given position
                    # and then there is no score
                    
                    # logger.warning('Reference mismatch or missing score in MNP at position {} of chr {}'.format(pos, mutation['CHROMOSOME']))
                    pass

                    # logger.warning('Reference mismatch MNP at {}:{} {} {}'.format(mutation['CHROMOSOME'], pos,
                    #                                                                 [mutation['REF'], mutation['ALT']] , nucleotides))
                    # for v in values:
                    #     logger.warning('Inner issue report: {}:{} {} {} {}'.format(mutation['CHROMOSOME'], pos,
                    #                                                                 [mutation['REF'], mutation['ALT']], nucleotides, [v.ref, v.alt] ))

            if not mnp_scores:
                logger.warning('Missing scores for MNP at position {} of chr {}'.format(pos, mutation['CHROMOSOME']))
                return
            
            mutation['SCORE'] = max(mnp_scores)
            mutation['POSITION'] = pos + mnp_scores.index(mutation['SCORE'])

        elif mutation['ALT_TYPE'] == "indel":

            score = self.indels.get_indel_score(mutation)

            mutation['SCORE'] = score if not math.isnan(score) else None

    def run(self):
        """
        Loads the scores and compute the statistics for the observed mutations.
        Perform simulations to compute a set of randomized scores.
        By comparing those scores with the observed through different
        statistical tests, a p-value can be obtained.

        """
        # as this method in run as a separate process, we need to fix a seed here too
        np.random.seed(self.seed)

        # Load element scores
        self.scores = Scores(self.name, self.segments, self.score_config)

        if self.mutability_info:
            self.mutability = Mutabilities(self.name, self.segments, self.mutability_config)


        if self.use_indels:
            self.indels = Indel(self.scores)

        # Compute observed mutations statistics and scores
        self.result = self.compute_muts_statistics()

        if len(self.result['mutations']) > 0:
            statistic_test = STATISTIC_TESTS.get(self.statistic_name)

            positions = self.scores.get_all_positions()

            # if the mutability is provided,
            # load it and check that the positions are the same ones as in the scores
            if self.mutability_info:
                mutability_positions = self.mutability.get_all_positions()
                if set(mutability_positions) != set(positions):
                    logger.warning("Mutability positions differ from the reference positions.")
                    logger.warning("Mutability information will not be used in this element {}.".format(self.name))
                    self.mutability_info = False

            observed = []
            subs_scores = []
            subs_probs = []
            subs_depths = []
            subs_mutability = []

            indels_scores = []
            indels_probs = []

            probs_of_subs = sig2probs.build(self.signature, self.signature_column)

            indels_simulated_as_subs = 0

            for mut in self.result['mutations']:
                observed.append(mut['SCORE'])  # Observed mutations

                # Indels treated as subs also count for the subs probs
                if mut['ALT_TYPE'] == 'indel' and self.indels.simulated_as_subs:
                    self.p_subs = 1
                    self.p_indels = 0
                    probs_of_subs.add_observed(mut)
                # When only in frame indels are simulated as subs
                elif mut['ALT_TYPE'] == 'indel' and self.indels.in_frame_simulated_as_subs:
                    if max(len(mut['REF']), len(mut['ALT'])) % 3 == 0:
                        indels_simulated_as_subs += 1
                        probs_of_subs.add_observed(mut)
                elif mut['ALT_TYPE'] == 'indel':
                    pass
                else:  # SNP or MNP
                    # Count how many signature ids are and prepare a vector for each
                    # IMPORTANT: this implies that only the signature of the observed mutations is taken into account
                    probs_of_subs.add_observed(mut)

            # When the probabilities of subs and indels are None, they are taken from
            # mutations seen in the gene
            if self.p_subs is None or self.p_indels is None: # use the probabilities based on observed mutations
                self.p_subs = (self.result['snps'] + self.result['mnps'] + indels_simulated_as_subs) / len(self.result['mutations'])
                self.p_indels = 1 - self.p_subs

            # Compute the values for the substitutions
            if self.p_subs > 0:
                if self.mutability_info:
                    for pos in positions:
                        for s, m in zip(self.scores.get_score_by_position(pos),
                                            self.mutability.get_mutability_by_position(pos)):
                            subs_scores.append(s.value)
                            probs_of_subs.add_background(s.change)

                            subs_mutability.append(m.value)

                            # if that position is not in the list of mutability, use value 0
                            # mutab = self.mutability.get(str(pos), 0)
                            # subs_mutability.append( mutab[num] if type(mutab) == list else 0 )

                elif self.depths_info:
                    for pos in positions:
                        for s in self.scores.get_score_by_position(pos):
                            subs_scores.append(s.value)
                            probs_of_subs.add_background(s.change)
                            # if that position is not in the list of depths, use value 0
                            subs_depths.append(self.depths.get(pos, 0)) 
                else:
                    for pos in positions:
                        for s in self.scores.get_score_by_position(pos):
                            subs_scores.append(s.value)
                            probs_of_subs.add_background(s.change)

                if probs_of_subs.size > 0:
                    subs_probs = probs_of_subs.probs
                    if sum(subs_probs) == 0.0:
                        logger.warning('Probability of substitutions equal to 0 in {}'.format(self.name))
                        self.result['partitions'] = []
                        self.result['sampling_size'] = self.sampling_size
                        self.result['obs'] = None
                        self.result['neg_obs'] = None
                        return self
                    subs_probs = subs_probs * self.p_subs
                    subs_probs = list(subs_probs)
                else:  # Same prob for all values (according to the prob of substitutions)
                    subs_probs = [self.p_subs / len(subs_scores)] * len(subs_scores) if len(subs_scores) > 0 else []

            if self.use_indels and self.p_indels > 0:
                indels_scores = self.indels.get_background_indel_scores()

                # All indels have the same probability
                indels_probs = [self.p_indels/len(indels_scores)] * len(indels_scores) if len(indels_scores) > 0 else []

            simulation_scores = subs_scores + indels_scores
            simulation_probs = subs_probs + indels_probs

            # if mutability is provided use it for
            # adjusting the probabilities of sampling positions
            if self.mutability_info:
                
                # only SNVs provided as input
                if len(simulation_probs) == len(subs_mutability):
                    # if np.sum(subs_mutability) > 0:
                    simulation_probs = np.nan_to_num(np.array(simulation_probs, dtype = np.float64) * np.array(subs_mutability, dtype = np.float64))
                    simulation_probs = list(np.nan_to_num(simulation_probs / simulation_probs.sum()))
                    # logger.info("SNVs probabilities corrected by mutability")
                    # else:
                    #     logger.info("ERROR")
                        # print(subs_mutability)
                        # print(self.segments)

                # SNVs and indels provided as input
                elif len(simulation_probs) == 2 * len(subs_mutability):
                    simulation_probs_1st_half = simulation_probs[:len(simulation_probs)//2]
                    simulation_probs_2nd_half = simulation_probs[:len(simulation_probs)//2]
                    simulation_probs_1st_half_total = sum(simulation_probs_1st_half)
                    
                    simulation_probs_1st_half = np.nan_to_num(np.array(simulation_probs_1st_half, dtype = np.float64) * np.array(subs_mutability, dtype = np.float64))
                    simulation_probs_1st_half = np.nan_to_num(simulation_probs_1st_half / simulation_probs_1st_half.sum())
                    simulation_probs_1st_half_norm = simulation_probs_1st_half * simulation_probs_1st_half_total
                    
                    simulation_probs = list(simulation_probs_1st_half_norm) + list(simulation_probs_2nd_half)
                    
                    simulation_probs = list(np.nan_to_num(simulation_probs / np.sum(simulation_probs)))
                    # logger.info("SNVs probabilities corrected by mutability")


            # if depths are provided use them for
            # adjusting the probabilities of sampling positions
            elif self.depths_info:
                if len(subs_probs) == len(subs_depths):
                    subs_probs = np.nan_to_num(np.array(subs_probs) * np.array(subs_depths))
                    subs_probs = list(np.nan_to_num(subs_probs / np.sum(subs_probs)))
                    logger.info("SNVs probabilities corrected by depth")
                
                elif len(simulation_probs) == 2 * len(subs_depths):
                    simulation_probs_1st_half = simulation_probs[:len(simulation_probs)//2]
                    simulation_probs_2nd_half = simulation_probs[:len(simulation_probs)//2]
                    simulation_probs_1st_half_total = sum(simulation_probs_1st_half)
                    
                    simulation_probs_1st_half = np.nan_to_num(np.array(simulation_probs_1st_half, dtype = np.float64) * np.array(subs_depths, dtype = np.float64))
                    simulation_probs_1st_half = np.nan_to_num(simulation_probs_1st_half / simulation_probs_1st_half.sum())
                    simulation_probs_1st_half_norm = simulation_probs_1st_half * simulation_probs_1st_half_total
                    
                    simulation_probs = list(simulation_probs_1st_half_norm) + list(simulation_probs_2nd_half)
                    
                    simulation_probs = list(np.nan_to_num(simulation_probs / np.sum(simulation_probs)))

                    logger.info("SNVs probabilities corrected by mutability")



            muts_count = len(self.result['mutations'])
            chunk_count = (self.sampling_size * muts_count) // self.sampling_chunk
            chunk_size = self.sampling_size if chunk_count == 0 else self.sampling_size // chunk_count

            # Calculate sampling parallelization partitions
            self.result['partitions'] = partitions_list(self.sampling_size, chunk_size)
            self.result['sampling'] = {}

            # Run first partition
            first_partition = self.result['partitions'].pop(0)
            background = np.random.choice(simulation_scores, size=(first_partition, muts_count), p=simulation_probs, replace=True)
            obs, neg_obs, back_means = statistic_test.calc_observed(background, np.array(observed))
            self.result['obs'] = obs
            self.result['neg_obs'] = neg_obs
            self.result['back_means'] = list(back_means)

            # Sampling parallelization (if more than one partition)
            if len(self.result['partitions']) > 0 or obs < self.min_obs:
                self.result['muts_count'] = muts_count
                self.result['simulation_scores'] = simulation_scores
                self.result['simulation_probs'] = simulation_probs
                self.result['observed'] = observed
                self.result['statistic_name'] = self.statistic_name

        self.result['sampling_size'] = self.sampling_size
        self.result['symbol'] = self.symbol

        return self
