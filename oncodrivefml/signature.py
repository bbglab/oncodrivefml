"""
This module contains information related with the signature.

The signature is a way of assigning probabilities to certain mutations that have some
relation amongst them (e.g. cancer type, sample...).
This relation is identified by the **signature_id**.

The ``classifier`` parameter in the :ref:`configuration <project configuration>` of the signature
specifies which column of the mutations file (:data:`~oncodrivefml.load.MUTATIONS_HEADER`) is used as
the identifier for the different signature groups.
If not provided, all mutations contribute to one global signature.

The probabilities are taken only from substitutions. For them, the two bases that
surround the mutated one are taken into account. This is called the triplet.
For a certain mutation in a position *x* the reference triplet is the base in the
reference genome in position *x-1*, the base in *x* and the base in the *x+1*. The altered triplet
of the same mutation is equal for the bases in *x-1* and *x+1* but the base in *x* is the one
observed in the mutation.


.. _signature dict:

signature (:obj:`dict`)

    .. code-block:: python

        { signature_id:
            {
                (ref_triplet, alt_triplet): prob
            }
        }

"""
import collections
import logging
from os import path

from bgsignature.file import load as load_signature
from bgsignature.count.mutation import count as signature_count

from oncodrivefml import __logger_name__, reference
from oncodrivefml.error import OncodriveFMLError

logger = logging.getLogger(__logger_name__)


def collapse(counts):
    d = collections.defaultdict(int)
    for k, v in counts.items():
        d[k] += v
        d[reference.reverse_complementary_sequence(k)] += v
    return d


def load(file):
    return load_signature(file)


def compute(mutations, method, classifier=None, normalize=None):
    genome = reference.get_build()
    logger.info("Computing signatures")
    counts = signature_count(mutations, genome=genome, size=3, group=classifier)
    if method != 'full':
        counts = counts.collapse()

    if normalize is None:
        return counts.sum1()
    elif not genome.startswith('hg'):
        logger.warning('Cannot normalize genome %s. Use bgsignature package to build your custom signatures', genome)
        return counts.sum1()
    else:
        if normalize in ('whole_genome', 'wgs'):
            file = path.join(path.dirname(__file__), "genome.json.gz")
        elif normalize in ('whole_exome', 'wes', 'wxs'):
            file = path.join(path.dirname(__file__), "exome.json.gz")
        else:
            raise OncodriveFMLError('Siganture normalization {} not valid'.format(normalize))
        normalization_counts = load(file)
        if method != 'full':
            normalization_counts = collapse(normalization_counts)
        return counts.normalize(normalization_counts)
