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

#. Using the pattern ('pattern'):

   The indel produces a set of changes with a limited impact (in space).
   To do so, a window size is specified, and the changes within that
   window are analysed.

   The functional impact of the observed indel is taken as if it was a
   set of substitutions.
   For the background, the same pattern is applied to all possible position
   and a functional impact derived from the changes it produces.


"""

import math
import random
import numpy as np
from math import exp
from enum import Enum

from oncodrivefml.signature import get_ref


# Global variables
window_size = 10
weight = None
weighting_function = lambda x: 1
frame_length = 1
indels_max_repeats = 0
stop_function = None


complements_dict = {"A": "T", "T": "A", "G": "C", "C": "G", 'N': 'N'}

transition_dict = {"A": "G", "G": "A", "T": "C", "C": "T", 'N': 'N'}

transversion_dict = {"A": "C", "C": "A", "T": "G", "G": "T", 'N': 'N'}



def _init_indels(indels_config):
    """
    Initialize the indels module

    Args:
        indels_config (dict): configuration of how to compute the impact of indels

    """
    global window_size, weight, weighting_function, frame_length, stop_function, indels_max_repeats
    if indels_config['method'] == 'pattern':
        window_size = indels_config['window_size']
        function = indels_config['weight_function']
        weight = Weight(window_size, function)
        weighting_function = weight.function

    elif indels_config['method'] == 'stop':
        stop = StopsScore(indels_config['mean'])
        stop_function = stop.function

    if indels_config['enable_frame']:
        frame_length = 3

    indels_max_repeats = indels_config['max_repeats']


class Weight:
    """
    Contain the function that weights the score of the indel
    when the pattern method is used

    Args:
        length (int): expected length of the window
        type (str): function identifier

    Functions are implemented to work on the interval between 0 and the window length.
    This means that the function is close to 1 at x=0 and close to 0 at x=length.

    If x = actual position of the element in the window you get the original function value
    If x = actual position of the element in the window - indel_size you get the function value starting at the first element outside the indel size
    If x = actual position of the element in the window - indel_size +1 you get the function value starting at the last element of the indel_size
    """
    def __init__(self, length, type):
        self.length = length
        if type == 'constant':
            self.funct = self.constant
        if type =='linear':
            self.intercept = 1
            self.slope = -1.0/length
            self.funct = self.linear
        if type == 'logistic':
            self.beta = 0.5
            self.funct = self.logistic

    def function(self, x):
        return self.funct(x)

    def constant(self, x):
        return 1

    def linear(self, x):
        return self.slope * x + self.intercept

    def logistic(self, x):
        return 1.0/(1+exp(self.beta*(x-self.length/2)))


class StopsScore:
    def __init__(self, type):
        """
        Contain the function that operates on the scores of the stops
        when the indel is treated as stop

        Args:
            type (str): indentifier of the function

        """
        if type == 'mean':
            self.funct = self.mean
        elif type =='median':
            self.funct = self.median
        elif type == 'random':
            self.funct = self.random
        elif type == 'random_choice':
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
        return random.uniform(p1, p2)

    def choose(self, x):
        return random.choice(x)





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
        has_positive_strand (bool): if the element being analysed has positive strand
    """

    def __init__(self, scores, signature, signature_id, method, has_positive_strand=True):
        self.scores = scores
        self.signature = signature
        self.signature_id = signature_id
        self.has_positive_strand = has_positive_strand

        if method == 'pattern':
            self.get_background_indel_scores = self.get_background_indel_scores_from_pattern
            self.get_indel_score = self.get_indel_score_from_pattern
        elif method == "stop":
            self.get_background_indel_scores = self.get_background_indel_scores_as_stops
            self.get_indel_score = self.get_indel_score_from_stop
            self.scores.get_stop_scores()
        elif method == 'max':
            self.get_background_indel_scores = self.get_background_indel_scores_as_substitutions_without_signature
            self.get_indel_score = self.get_indel_score_max_of_subs


    @staticmethod
    def is_frameshift(size):
        """

        Args:
            size (int): length of the indel

        Returns:
            bool. Whether the size is multiple of 3 (in the frames have been
            enabled in the configuration)

        """
        if frame_length != 1 and (size % frame_length) == 0:
            return False
        return True

    @staticmethod
    def compute_window_size(indel_size):
        """
        From an size of an indel compute size of the
        window frame to analyse

        Args:
            indel_size (int): length of the indel

        Returns:
            int.

        """
        size = window_size

        if not Indel.is_frameshift(indel_size):  # in-frame
            size = (indel_size + frame_length - 1)

        if indel_size > size:
            size = indel_size

        return size


    @staticmethod
    def is_in_repetitive_region(mutation):
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
            bool. Whether the indel falls in a repetitive region or to

        """

        if indels_max_repeats == 0:
            # 0 means do not search for repetitive regions
            return False

        chrom = mutation['CHROMOSOME']
        ref = mutation['REF']
        alt = mutation['ALT']
        pos = mutation['POSITION']

        # Check if it's repeated
        seq = alt if '-' in ref else ref
        size = indels_max_repeats * 2 * len(seq)
        ref = get_ref(chrom, pos-(size//2)-1, size+1)
        return ref.count(seq) >= indels_max_repeats


    def get_mutation_sequences(self, mutation, size):
        """
        Get the reference and altered sequence of the indel
        along the window size

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`
            size (int): window length

        Returns:
            tuple. Reference and alterned sequences

        """
        # TODO pass indel_size as parameter because it is already computed in the get_indel_score method
        position = mutation['POSITION']
        is_insertion = True if '-' in mutation['REF'] else False
        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        if self.has_positive_strand:
            reference = get_ref(mutation['CHROMOSOME'], position, indel_size + size)  # ensure we have enough values
            # TODO check if we are looking for values outside the element
            if is_insertion:
                alteration = mutation['ALT'] + reference
            else:
                alteration = reference[indel_size:]
            return reference[:size], alteration[:size]
        else:  # negative strand
            reference = get_ref(mutation['CHROMOSOME'], position - (indel_size + size), indel_size + size)  # ensure we have enough values
            if is_insertion:
                alteration = reference + mutation['ALT']
            else:
                alteration = reference[:-indel_size]  # remove the last indels_size elements
            return reference[-size:], alteration[-size:]  # return only last size elements

    def weight(self, score_values, indel_size, total_size):
        """
        Apply the weight function to the measured scores

        The function is only applied to values beyond the lenght of the indel
        if the indel is a frameshift

        Args:
            score_values (list): scores of the indel pattern
            indel_size (int): length of the indel
            total_size (int): length of the window

        Returns:
            list. Values with the function applied

        """
        if Indel.is_frameshift(indel_size):
            if self.has_positive_strand:
                for i in range(indel_size, total_size):  # if indel_size is bigger than the window it is an empty range
                    # We can pass to the weight function the value i or the value i-indel_size+1. The first means using the value of the function as it starts in 0, the second as it start when the indel ends
                    score_values[i] *= weighting_function(i - indel_size + 1)  # first element has x = 1
            else:
                for i in range(total_size - indel_size):
                    score_values[i] *= weighting_function((total_size - indel_size) - i)

        return score_values

    def compute_scores(self, reference, alternation, initial_position, size):
        """
        Compute the scores of all substitution between the reference and altered sequences

        Args:
            reference (str): sequence
            alternation (str): sequence
            initial_position (int): position where the indel occurs
            size (int): number of position to look

        Returns:
            list. Scores of the substitution in the indel. :obj:`~math.nan` when it is not possible
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

    def get_indel_score_from_pattern(self, mutation):
        """
        Compute the score of an indel from its pattern

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`

        Returns:
            float. Impact of the indel

        """

        if Indel.is_in_repetitive_region(mutation):
            return math.nan

        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        position = int(mutation['POSITION'])

        indel_window_size = Indel.compute_window_size(indel_size)

        ref, alt = self.get_mutation_sequences(mutation, indel_window_size)

        if self.has_positive_strand:
            init_pos = position
        else:
            init_pos = position - indel_window_size

        indel_scores = self.compute_scores(ref, alt, init_pos, indel_window_size)

        cleaned_scores = [score for score in indel_scores if not math.isnan(score)]
        return max(cleaned_scores) if cleaned_scores else math.nan

    def get_indel_score_max_of_subs(self, mutation):
        """
        Compute the score of an indel by treating each alteration as
        a substitution.

        Is similar to :meth:`get_indel_score_from_pattern` but forcing the window size
        to be the same as the indel length.

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`

        Returns:
            float. Maximum value of all substitution

        """

        if Indel.is_in_repetitive_region(mutation):
            return math.nan

        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        position = int(mutation['POSITION'])

        indel_window_size = indel_size

        ref, alt = self.get_mutation_sequences(mutation, indel_window_size)

        if self.has_positive_strand:
            init_pos = position
        else:
            init_pos = position - indel_window_size

        indel_scores = self.compute_scores(ref, alt, init_pos, indel_window_size)

        indel_scores = self.weight(indel_scores, indel_size, indel_window_size)

        cleaned_scores = [score for score in indel_scores if not math.isnan(score)]
        return max(cleaned_scores) if cleaned_scores else math.nan

    def get_indel_score_substitutions(self, mutation, mutation_position):
        indel_size = max(len(mutation['REF']), len(mutation['ALT']))

        indel_scores = []
        for pos in range(mutation_position, mutation_position + indel_size):
            indel_scores += [s.value for s in self.scores.get_score_by_position(pos)]

        return max(indel_scores)

    def get_indel_score_for_background_from_pattern(self, position, length, chromosome, pattern, indel_size):
        """
        Compute the background score applying the pattern of the mutation to a position

        Args:
            position (int): position where the mutation is simulated
            length (int): window size
            chromosome (str): chromosome
            pattern (dict): dict with the changes
            indel_size (int): size of the indel

        Returns:

        """
        # position is the starting position of the sequence where to apply the pattern
        # no matter if it is positve or negative strand
        reference = get_ref(chromosome, position, length)
        alteration = Indel.apply_pattern(reference, pattern)

        indel_scores = self.compute_scores(reference, alteration, position, length)

        indel_scores = self.weight(indel_scores, indel_size, length)

        cleaned_scores = [score for score in indel_scores if not math.isnan(score)]
        return max(cleaned_scores) if cleaned_scores else math.nan

    def get_background_indel_scores_from_pattern(self, mutations):
        """
        Compute the background score (:meth:`get_indel_score_for_background_from_pattern`)
        for all possible positions

        Args:
            mutations (dict): a mutation object as in :ref:`here <mutations dict>`

        Returns:
            list. Indel scores

        """
        indel_scores = []
        for mutation in mutations:
            mut_scores = []
            indel_size = max(len(mutation['REF']), len(mutation['ALT']))
            length = Indel.compute_window_size(indel_size)
            mutation_pattern = self.get_pattern(mutation, length)

            for pos in self.scores.get_all_positions():
                score = self.get_indel_score_for_background_from_pattern(pos, length, mutation['CHROMOSOME'], mutation_pattern,
                                                            indel_size)

                if not math.isnan(score):
                    mut_scores.append(score)

            indel_scores += mut_scores

        if indel_scores == []:#TODO decide the best way to hadle this
            # When no mutations are given, simulate the results as substitutions
            indel_scores = self.get_background_indel_scores_as_substitutions_without_signature()

        return indel_scores

    @staticmethod
    def compute_pattern(reference, alternate, length):
        """
        Obtain the pattenr of changes in a sequence

        Args:
            reference (str):
            alternate (str):
            length (int):

        Returns:
            list.

        """
        pattern = []

        for i in range(length):
            if reference[i] == alternate[i]:
                pattern.append(ChangePattern.same)
            elif reference[i] == complements_dict[alternate[i]]:
                pattern.append(ChangePattern.complementary)
            elif reference[i] == transition_dict[alternate[i]]:
                pattern.append(ChangePattern.transition)
            else:
                pattern.append(ChangePattern.transversion)

        return pattern

    def get_pattern(self, mutation, length):
        """
        Obtain the pattern that an indel produces

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`
            length (int): window

        Returns:
            list.

        """
        ref, alt = self.get_mutation_sequences(mutation, length)
        return Indel.compute_pattern(ref, alt, length)

    @staticmethod
    def apply_pattern(sequence, pattern):
        """
        Apply a pattern to a sequence

        Args:
            sequence (str): genome sequence
            pattern (list):

        Returns:
            str. Altered sequence

        """
        new_seq = ''
        for i in range(len(sequence)):
            new_seq += pattern_change[pattern[i]](sequence[i])
        return new_seq

    def get_indel_score_from_stop(self, mutation):
        """
        Compute the indel score as a stop

        A function is applied to the values of the scores in the gene

        Args:
            mutation (dict): a mutation object as in :ref:`here <mutations dict>`

        Returns:
            float. Score value. :obj:`~math.nan` if is not possible to compute it

        """

        if Indel.is_in_repetitive_region(mutation):
            return math.nan

        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        if Indel.is_frameshift(indel_size):
            return stop_function(self.scores.stop_scores)
        else:
            position = int(mutation['POSITION'])
            ref, alt = self.get_mutation_sequences(mutation, indel_size)

            if self.has_positive_strand:
                init_pos = position
            else:
                init_pos = position - indel_size
            indel_scores = self.compute_scores(ref, alt, init_pos, indel_size)

            cleaned_scores = [score for score in indel_scores if not math.isnan(score)]
            return max(cleaned_scores) if cleaned_scores else math.nan


    # def get_background_indel_scores_as_substitutions(self, mutation, positions):
    #     indel_scores = []
    #     signatures = []
    #     for pos in positions:
    #         for s in self.scores.get_score_by_position(pos):
    #             indel_scores.append(s.value)
    #             if self.signature is not None:
    #                 signatures.append(
    #                     self.signature[mutation.get(self.signature_id, self.signature_id)].get(
    #                         (s.ref_triplet, s.alt_triplet), 0.0))
    #     return indel_scores, signatures

    def get_background_indel_scores_as_substitutions_without_signature(self, **kwargs):
        """
        Return the values of scores of all possible substitutions

        Args:
            **kwargs:

        Returns:
            list.

        """
        # TODO if these values have already been computed for the subs, we can optimize the
        # code avoiding this second loop
        indel_scores = []
        for pos in self.scores.get_all_positions():
            for s in self.scores.get_score_by_position(pos):
                indel_scores.append(s.value)
        return indel_scores

    def get_background_indel_scores_as_stops(self, **kwargs):
        """
        Args:
            **kwargs:

        Returns:
            list. Values of the stop scores of the gene

        """
        return self.scores.stop_scores


    def not_found(self, mutation):
        return np.nan

class ChangePattern(Enum):
    same = 1,
    complementary = 2,
    transition = 3,
    transversion = 4

pattern_change = {ChangePattern.same: lambda x: x,
                  ChangePattern.complementary: lambda x: complements_dict[x],
                  ChangePattern.transition: lambda x: transition_dict[x],
                  ChangePattern.transversion: lambda x: transversion_dict[x]
                  }
