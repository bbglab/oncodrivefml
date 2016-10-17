import math
from enum import Enum

import logging

from oncodrivefml.signature import get_ref
from math import exp

window_size = 10
weight = None
weighting_function = lambda x: 1
in_frame_length = 1
#TODO rename as inframe


complements_dict = {"A": "T", "T": "A", "G": "C", "C": "G", 'N': 'N'}

transition_dict = {"A": "G", "G": "A", "T": "C", "C": "T", 'N': 'N'}

transversion_dict = {"A": "C", "C": "A", "T": "G", "G": "T", 'N': 'N'}



def _init_indels(indels_config):
    global window_size, weight, weighting_function, in_frame_length
    if indels_config['method'] == 'pattern':
        window_size = indels_config.get('window_size', 10)
        function = indels_config.get('weight_function', 'constant')
        weight = Weight(window_size, function)
        weighting_function = weight.function
        if indels_config.get('frame_shift', False):
            in_frame_length = 3


class Weight:
    """
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


class Indel:

    def __init__(self, scores, method, has_positive_strand=True):
        self.scores = scores
        self.has_positive_strand = has_positive_strand

        if method == 'pattern':
            self.get_background_indel_scores = self.get_background_indel_scores_from_pattern
            self.get_indel_score = self.get_indel_score_from_pattern

    @staticmethod
    def compute_window_size(indel_size):
        size = window_size

        if in_frame_length != 1 and (indel_size % in_frame_length) == 0:  # in-frame
            size = (indel_size + in_frame_length - 1)

        if indel_size > size:
            size = indel_size

        return size

    def get_mutation_sequences(self, mutation, size):
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

    def weight(self, scores, indel_size, total_size):
        if in_frame_length == 1 or (indel_size % in_frame_length) != 0:

            if self.has_positive_strand:
                for i in range(indel_size, total_size):  # if indel_size is bigger than the window it is an empty range
                    # We can pass to the weight function the value i or the value i-indel_size+1. The first means using the value of the function as it starts in 0, the second as it start when the indel ends
                    scores[i] *= weighting_function(i - indel_size + 1)  # first element has x = 1 TODO check
            else:
                for i in range(total_size - indel_size):
                    scores[i] *= weighting_function((total_size - indel_size) - i)

        return scores

    def compute_scores(self, reference, alternation, initial_position, size):
        computed_scores = [math.nan] * size
        for index, position in enumerate(range(initial_position, initial_position + size)):
            scores_in_position = self.scores.get_score_by_position(position)
            for score in scores_in_position:
                if score.ref == reference[index] and score.alt == alternation[index]:
                    computed_scores[index] = score.value
                    break
        return computed_scores

    def get_indel_score_from_pattern(self, mutation):
        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        position = int(mutation['POSITION'])

        indel_window_size = Indel.compute_window_size(indel_size)

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
        # position is the starting position of the sequence where to apply the pattern
        # no matter if it is positve or negative strand
        reference = get_ref(chromosome, position, length)
        alteration = Indel.apply_pattern(reference, pattern)

        indel_scores = self.compute_scores(reference, alteration, position, length)

        indel_scores = self.weight(indel_scores, indel_size, length)

        cleaned_scores = [score for score in indel_scores if not math.isnan(score)]
        return max(cleaned_scores) if cleaned_scores else math.nan

    def get_background_indel_scores_from_pattern(self, mutation, positions):
        scores = []
        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        length = Indel.compute_window_size(indel_size)
        mutation_pattern = self.get_pattern(mutation, length)

        for pos in positions:
            score = self.get_indel_score_for_background_from_pattern(pos, length, mutation['CHROMOSOME'], mutation_pattern,
                                                        indel_size)

            if not math.isnan(score):
                scores.append(score)

        return scores

    @staticmethod
    def compute_pattern(reference, alternate, length):
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
        ref, alt = self.get_mutation_sequences(mutation, length)
        return Indel.compute_pattern(ref, alt, length)

    @staticmethod
    def apply_pattern(sequence, pattern):
        new_seq = ''
        for i in range(len(sequence)):
            new_seq += pattern_change[pattern[i]](sequence[i])
        return new_seq


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
