import math
from oncodrivefml.signature import get_ref
from math import exp

window_size = 10
weight = None
weighting_function = lambda x: 1
frame_shift = 1


def _init_indels(lenght=10, method='cte', shift=1):
    global window_size, weight, weighting_function, frame_shift
    window_size = lenght
    weight = Weight(window_size, method)
    weighting_function = weight.function
    frame_shift = shift


class Weight:
    '''
    Functions are implemented to work on the interval between 0 and the window length.
    This means that the function is close to 1 at x=0 and close to 0 at x=length.

    If x = actual position of the element in the window you get the original function value
    If x = actual position of the element in the window - indel_size you get the function value starting at the first element outside the indel size
    If x = actual position of the element in the window - indel_size +1 you get the function value starting at the last element of the indel_size

    '''
    def __init__(self, length, type):
        self.length = length
        if type == 'cte':
            self.funct = self.cte
        if type =='linear':
            self.intercept = 1
            self.slope = -1.0/length
            self.funct = self.linear
        if type == 'logistic':
            self.beta = 10
            self.funct = self.logistic

    def function(self, x):
        return self.funct(x)

    def cte(self, x):
        return 1

    def linear(self, x):
        return self.slope * x + self.intercept

    def logistic(self, x):
        return 1.0/(1+exp(self.beta*(x-self.length/2)))


class Indel:

    @staticmethod
    def get_indel_score(mutation, scores, mutation_position, positive_strand):
        indel_size = max(len(mutation['REF']), len(mutation['ALT']))
        is_insertion = True if '-' in mutation['REF'] else False


        indel_window_size = window_size

        if frame_shift != 1 and (indel_size % frame_shift) == 0:  # frame shift
            indel_window_size = (indel_size + frame_shift-1) if (indel_size + frame_shift-1) <= window_size else window_size

        if indel_size > indel_window_size:
            indel_window_size = indel_size

        indel_scores = [math.nan] * indel_window_size

        if positive_strand:
            for index, position in enumerate(range(mutation_position, mutation_position + indel_window_size)):
                scores_in_position = scores.get_score_by_position(position)  # if position is outside the element, it has not score so the indel score remains nan
                for score in scores_in_position:
                    if is_insertion:
                        alteration = mutation['ALT'][index] if index < indel_size else get_ref(mutation['CHROMOSOME'], position - indel_size)
                    else:  # deletion
                        try:
                            alteration = get_ref(mutation['CHROMOSOME'], position + indel_size)
                        except ValueError:  # generated when we are outside the element
                            break  # ignore it (score = nan)
                    if score.ref == get_ref(mutation['CHROMOSOME'], position) and score.alt == alteration:
                        indel_scores[index] = score.value
                        break
        else:
            for index, position in enumerate(range(mutation_position-indel_window_size, mutation_position)):
                scores_in_position = scores.get_score_by_position(position)  # if position is outside the element, it has not score so the indel score remains nan
                for score in scores_in_position:
                    if is_insertion:
                        alteration = get_ref(mutation['CHROMOSOME'], position + indel_size) if index < (indel_window_size - indel_size) else mutation['ALT'][index - (indel_window_size - indel_size)]
                    else:  # deletion
                        try:
                            alteration = get_ref(mutation['CHROMOSOME'], position - indel_size)
                        except ValueError:  # generated when we are outside the element
                            break  # ignore it (score = nan)
                    if score.ref == get_ref(mutation['CHROMOSOME'], position) and score.alt == alteration:
                        indel_scores[index] = score.value
                        break


        if frame_shift == 1 or (indel_size % frame_shift) != 0:
            if positive_strand:
                for i in range(indel_size, indel_window_size):  # if indel_size is bigger than the window it is an empty range
                    # We can pass to the weight function the value i or the value i-indel_size+1. The first means using the value of the function as it starts in 0, the second as it start when the indel ends
                    indel_scores[i] *= weighting_function(i - indel_size + 1)  # first element has x = 1 TODO check
            else:
                for i in range(indel_window_size-indel_size):
                    indel_scores[i] *= weighting_function((indel_window_size - indel_size)-i)

        cleaned_scores = [score for score in indel_scores if not math.isnan(score)]
        return max(cleaned_scores) if cleaned_scores else math.nan
