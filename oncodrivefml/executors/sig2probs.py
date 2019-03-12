import collections

import numpy as np


class SubstitutionProbs:

    def __init__(self, signature):
        self._signature = signature
        self._probs = []

    def add_observed(self, mutation):
        pass

    def add_background(self, change):
        self._probs.append(self._signature.get(change, 0.0))

    @property
    def size(self):
        return len(self._probs)

    @property
    def probs(self):
        return np.array(self._probs) / sum(self._probs)


class NoSignature(SubstitutionProbs):

    def __init__(self):
        self._size = 0

    @property
    def size(self):
        return self._size

    def add_background(self, change):
        self._size += 1

    @property
    def probs(self):
        return np.array([1 / self._size] * self._size)


class GroupSignature(SubstitutionProbs):

    def __init__(self, signature, classifier):
        super().__init__(signature)
        self._size = 0
        self._classifier = classifier
        self._seen_signatures_counter = collections.defaultdict(int)
        self._probs_by_sig = collections.defaultdict(list)

    @property
    def size(self):
        return self._size

    def add_observed(self, mutation):
        signature_id = mutation[self._classifier]
        self._seen_signatures_counter[signature_id] += 1

    def add_background(self, change):
        for k in self._seen_signatures_counter.keys():
            value = self._signature[k].get(change, 0.0)
            self._probs_by_sig[k].append(value)
        self._size += 1

    @property
    def probs(self):
        total_ids = sum(v for v in self._seen_signatures_counter.values())
        probs = np.array([0.0] * self.size)
        for k, v in self._probs_by_sig.items():
            probs += (np.array(v) * self._seen_signatures_counter[k] / total_ids)
        total = sum(probs)
        return probs / total


def build(signature=None, classifier=None):
    if signature is None:
        return NoSignature()
    elif classifier is None:
        return SubstitutionProbs(signature)
    else:
        return GroupSignature(signature, classifier)
