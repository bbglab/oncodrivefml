"""
This module contains all utilities to process
insertions and deletions.

Currently 3 methods have been implemented to compute
the impact of the indels.


.. _indels methods:

1. As a set of substitutions ('max'):

   The indel is treated as set of substitutions.
   It is used for non-coding regions

   The functional impact of the observed mutation is the maximum
   of all the substitutions.
   The background is simulated as substitutions are.

#. As a stop ('stop'):

   The indel is expected to produce a stop in the genome,
   unless it is a frame-shift indel.
   It is used for coding regions.

   The functional impact is derived from the function impact
   of the stops of the gene.
   The background is simulated also as stops.

"""

import math

import numpy as np

from oncodrivefml.reference import get_ref


# Global variables
analysis = None
max_repeats = None
stop_function = None


def init_indels_module(indels_config):
    """
    Initialize the indels module

    Args:
        indels_config (dict): configuration of how to compute the impact of indels

    """
    global analysis, stop_function, max_repeats

    if indels_config['method'] == 'stop':
        analysis = 'coding'
    elif indels_config['method'] == 'max':
        analysis = 'noncoding'

    stop = StopsScore(indels_config['stops_function'])
    stop_function = stop.function

    max_repeats = indels_config['max_consecutive']


class StopsScore:
    def __init__(self, funct_type):
        """
        Contain the function that operates on the scores of the stops
        when the indel is treated as stop

        Args:
            funct_type (str): indentifier of the function

        """
        if funct_type == 'mean':
            self.funct = self.mean
        elif funct_type == 'median':
            self.funct = self.median
        elif funct_type == 'random':
            self.funct = self.random
        elif funct_type == 'random_choice':
            self.funct = self.choose

    def function(self, x):
        return self.funct(x)

    def mean(self, x):
        return np.mean(x)

    def median(self, x):
        return np.median(x)

    def random(self, x):
        p2 = max(x)
        p1 = min(x)
        return np.random.uniform(p1, p2)

    def choose(self, x):
        return np.random.choice(x)


class Indel:
    """
    Methods to compute the impact of indels
    for the observed and the background

    Args:
        scores (:class:`~oncodrivefml.scores.Scores`): functional impact per position
        signature (dict): see :ref:`signature <signature dict>`
        signature_id (str): classifier for the signatures
        method (str): identifies which method to use to compute the functional impact
            (see :ref:`methods <indels methods>`)
        strand (str): if the element being analysed has positive, negative or unknown strand (+,-,.)

    """

    def __init__(self, scores):
        self.scores = scores
        self.simulated_as_subs = False
        self.in_frame_simulated_as_subs = False

        if analysis == "coding":
            self.in_frame_simulated_as_subs = True
            self.get_indel_score = self.get_indel_score_from_stop
            self.get_background_indel_scores = self.get_background_indel_scores_as_stops
            self.scores.get_stop_scores()
        elif analysis == 'noncoding':
            self.simulated_as_subs = True
            self.get_indel_score = self.get_indel_score_max_of_subs
            self.get_background_indel_scores = self.get_background_indel_scores_as_substitutions_without_signature

    @staticmethod
    def is_frameshift(size):
        """

        Args:
            size (int): length of the indel

        Returns:
            bool. Whether the size is multiple of 3 (in the frames have been
            enabled in the configuration)

        """
        if (size % 3) == 0:
            return False
        return True

    def is_in_repetitive_region(self, mutation):
        """
        Check if  an indel falls in a repetitive region

        Looking in the window with the indel in the middle, check if the
        same sequence of the indel appears at least a certain number of times
        specified in the configuration.
        The window where to look has twice the size of the indel multiplied by
        the number of times already mentioned.

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`

        Returns:
            bool: Whether the indel falls in a repetitive region or not

        """

        if max_repeats == 0:
            # 0 means do not search for repetitive regions
            return False

        chrom = mutation['CHROMOSOME']
        ref = mutation['REF']
        alt = mutation['ALT']
        pos = mutation['POSITION']

        # Check if it's repeated
        seq = alt if '-' in ref else ref
        size = max_repeats * len(seq)
        ref = get_ref(chrom, pos, size)
        return ref.count(seq) >= max_repeats

    def get_mutation_sequences(self, mutation, size):
        """
        Get the reference and altered sequence of the indel
        along the window size

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`
            size (int): window length

        Returns:
            tuple: Reference and alternated sequences

        """
        position = mutation['POSITION']
        is_insertion = True if '-' in mutation['REF'] else False
        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        reference = get_ref(mutation['CHROMOSOME'], position, indel_size + size)  # ensure we have enough values
        # TODO check if we are looking for values outside the element
        if is_insertion:
            alteration = mutation['ALT'] + reference
        else:
            alteration = reference[indel_size:]
        return reference[:size], alteration[:size]

    def compute_scores(self, reference, alternation, initial_position, size):
        """
        Compute the scores of all substitution between the reference and altered sequences

        Args:
            reference (str): sequence
            alternation (str): sequence
            initial_position (int): position where the indel occurs
            size (int): number of position to look

        Returns:
            list: Scores of the substitution in the indel. :obj:`~math.nan` when it is not possible
            to compute a value.

        """
        computed_scores = [math.nan] * size
        for index, position in enumerate(range(initial_position, initial_position + size)):
            scores_in_position = self.scores.get_score_by_position(position)
            for score in scores_in_position:
                if score.ref == reference[index] and score.alt == alternation[index]:
                    computed_scores[index] = score.value
                    break
        return computed_scores

    def get_indel_score_max_of_subs(self, mutation):
        """
        Compute the score of an indel by treating each alteration as
        a substitution.

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`

        Returns:
            float: Maximum value of all substitutions

        """

        if self.is_in_repetitive_region(mutation):
            return math.nan

        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        position = int(mutation['POSITION'])

        ref, alt = self.get_mutation_sequences(mutation, indel_size)

        init_pos = position

        indel_scores = self.compute_scores(ref, alt, init_pos, indel_size)

        cleaned_scores = [score for score in indel_scores if not math.isnan(score)]
        return max(cleaned_scores) if cleaned_scores else math.nan

    def get_indel_score_from_stop(self, mutation):
        """
        Compute the indel score as a stop

        A function is applied to the values of the scores in the gene

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`

        Returns:
            float: Score value. :obj:`~math.nan` if is not possible to compute it

        """

        if self.is_in_repetitive_region(mutation):
            return math.nan

        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        if Indel.is_frameshift(indel_size):
            return stop_function(self.scores.stop_scores)
        else:
            position = int(mutation['POSITION'])
            ref, alt = self.get_mutation_sequences(mutation, indel_size)

            init_pos = position
            indel_scores = self.compute_scores(ref, alt, init_pos, indel_size)

            cleaned_scores = [score for score in indel_scores if not math.isnan(score)]
            return max(cleaned_scores) if cleaned_scores else math.nan

    def get_background_indel_scores_as_substitutions_without_signature(self):
        """
        Return the values of scores of all possible substitutions
        Returns:
            list.

        """
        indel_scores = []
        for pos in self.scores.get_all_positions():
            for s in self.scores.get_score_by_position(pos):
                indel_scores.append(s.value)
        return indel_scores

    def get_background_indel_scores_as_stops(self):
        """

        Returns:
            list: Values of the stop scores of the gene

        """
        return self.scores.stop_scores

    def not_found(self, mutation):
        return np.nan
