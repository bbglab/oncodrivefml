"""
This module contatins information related with the signature.

The signature is a way of assigning probabilities to certain mutations that have some
relation amongst them (e.g. cancer type, sample...).

This relation is identified by the **signature_id**.
It uses the *SIGNATURE* field in the `mutations dict`_.

.. warning::

    The value of the *SIGNATURE* in the `mutations dict`_ can be different than the value
     in the *SIGNATURE* column of the mutations file (see :class:`oncodrivefml.OncodriveFML`).

The ``classifier`` parameter in the :ref:`configuration <project configuration>` of the signature
specifies which column of the mutations file (:data:`oncodrivefml.load.MUTATIONS_HEADER`) is used to replace the
*SIGNATURE* column. If the column does not exist the ``classifier`` itself is used as value for the
*SIGNATURE*.

The probabilities are taken using only substitutions. For them, the two bases that
surround the mutated one are taken into account. This is called the triplet.
For a certain mutation in a position *x* the reference triplet is the base in the
refernce genome in position *x-1, the base in *x* and the base in the *x+1*. The altered triplet
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

import gzip
import logging
import os
import mmap
import pickle

import bgdata
import pandas as pd
from os.path import join, exists
from oncodrivefml.load import load_mutations
from collections import defaultdict

HG19 = None
"""
Path to the reference genome file.
See: :func:`get_hg19_dataset`.
"""

HG19_MMAP_FILES = {}
"""
Dictionary with chromosome as keys and memory maps with the reference genome as values.
See: :func:`get_hg19_mmap`.
"""

__CB = {"A": "T", "T": "A", "G": "C", "C": "G"}


def get_hg19_dataset():
    """
    Sets the path to the reference genome

    Returns:
        :attr:`HG19`

    """
    global HG19

    if HG19 is None:
        HG19 = bgdata.get_path('datasets', 'genomereference', 'hg19')

    return HG19


def get_hg19_mmap(chromosome):
    """
    Get a memory map with the reference genome of the chromosome.

    Args:
        chromosome (str): chromosome identifier (number or X or Y)

    Returns:
        mmap: memory map with the reference genome of the chromosome


    If the chromosome was not in :attr:`HG19_MMAP_FILES` before, it is added.
    """
    if chromosome not in HG19_MMAP_FILES:
        fd = open(join(get_hg19_dataset(), "chr{0}.txt".format(chromosome)), 'rb')
        HG19_MMAP_FILES[chromosome] = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    return HG19_MMAP_FILES[chromosome]


def get_ref_triplet(chromosome, start):
    """

    Args:
        chromosome (str): chromosome identifier
        start (int): starting position

    Returns:
        str: 3 bases from the reference genome

    """
    return get_ref(chromosome, start, size=3)


def get_ref(chromosome, start, size=1):
    """

    Args:
        chromosome (str): chromosome identifier
        start (int): starting position
        size (int): amount of bases. Default to 1.

    Returns:
        str: bases in the reference genome

    """
    mm_file = get_hg19_mmap(chromosome)
    mm_file.seek(start - 1)
    return mm_file.read(size).decode().upper()


def get_reference_signature(line):
    """

    Args:
        line (dict): contatins the chromosome and the position

    Returns:
        str: triplet around certain positions

    """
    return get_ref_triplet(line['CHROMOSOME'], line['POSITION'] - 1)


def get_alternate_signature(line):
    """

    Args:
        line (dict): contains the previous base, the alteration and the next base

    Returns:
        str: triplet with the central base replaced by the alteration indicated in the line

    """
    return line['Signature_reference'][0] + line['ALT'] + line['Signature_reference'][2]


def complementary_sequence(seq):
    """

    Args:
        seq (str): sequence of bases

    Returns:
        str: complementary sequence

    """
    return "".join([__CB[base] if base in __CB else base for base in seq.upper()])


def collapse_complementaries(signature):
    """
    Add to the amount of a certain pair (ref_triplet, alt_triplet) the amount of the complementary.

    Args:
        signature (dict): { (ref_triplet, alt_triplet): amount }

    Returns:
        dict: { (ref_triplet, alt_triplet): new_amount }. New_amount is the addition of the amount
        for (ref_triplet, alt_triplet) and the amount for (complementary_ref_triplet, complementary_alt_triplet)

    """
    comp_sig = defaultdict(int)
    for k, v in signature.items():
        comp_sig[k] += v
        comp_k = (complementary_sequence(k[0]), complementary_sequence(k[1]))
        comp_sig[comp_k] += v
    return comp_sig


def signature_probability(signature_counts):
    """
    Associates to each key (tuple(reference_tripet, altered_triplet)) the value divided by the total amount

    Args:
        signature_counts (dict): pair key-amount {(ref_triplet, alt_triplet): value}

    Returns:
        dict: pair key-(amount/total_amount)

    """
    total = sum([v for v in signature_counts.values()])
    return {k: v/total for k, v in signature_counts.items()}


def compute_signature(signature_function, classifier, blacklist, collapse=False):
    """
    Gets the probability of each substitution that occurs for a certain signature_id.

    Each substitution is identified by the pair (reference_triplet, altered_triplet).

    The signature_id is taken from the mutations ``SIGNATURE`` field.

    Args:
        signature_function: function that yields one mutation each time
        classifier (str): passed to :func:`oncodrivefml.load.load_mutations`
            as parameter ``signature_classifier``.
        blacklist: file with blacklisted samples (see :class:`oncodrivefml.main.OncodriveFML`).
            Used by :func:`oncodrivefml.load.load_mutations`
        collapse (bool): consider one substitutions and the complementary one as the same. Defaults to True.

    Returns:
        dict: probability of each substitution (measured by the triplets) grouped by the signature_classifier

        .. code-block:: python

            { signature_id:
                {
                    (ref_triplet, alt_triplet): prob
                }
            }

    .. warning::

        Only substitutions are taken into account

    """
    signature_count = defaultdict(lambda: defaultdict(int))
    for mut in signature_function():
        if mut['TYPE'] != 'subs':
            continue

        signature_ref = get_ref_triplet(mut['CHROMOSOME'], mut['POSITION'] - 1)
        signature_alt = signature_ref[0] + mut['ALT'] + signature_ref[2]

        signature_count[mut.get(classifier, classifier)][(signature_ref, signature_alt)] += 1

    signature = {}
    for k, v in signature_count.items():
        if collapse:
            signature[k] = signature_probability(collapse_complementaries(v))
        else:
            signature[k] = signature_probability(v)


    return signature


def load_signature(mutations_file, signature_function, signature_config, blacklist=None, save_pickle=False):
    """
    Computes the probability that certain mutation occurs.

    Args:
        mutations_file: mutations file
        signature_function: function that yields one mutation each time
        signature_config (dict): information of the signature (see :ref:`configuration <project configuration>`)
        blacklist (optional): file with blacklisted samples (see :class:`oncodrivefml.main.OncodriveFML`). Defaults to None.
            Used by :func:`oncodrivefml.load.load_mutations`
        save_pickle (:obj:`bool`, optional): save pickle files

    Returns:
        dict: probability of each substitution (measured by the triplets) grouped by the signature_id

        .. code-block:: python

            { signature_id:
                {
                    (ref_triplet, alt_triplet): prob
                }
            }

    Before computing the signature, it is checked whether a pickle file with the signature already exists or not.

    """
    method = signature_config['method']
    classifier = signature_config['classifier']
    path = signature_config.get('path', None)
    column_ref = signature_config.get('column_ref', None)
    column_alt = signature_config.get('column_alt', None)
    column_probability = signature_config.get('column_probability', None)

    if path is not None and path.endswith(".pickle.gz"):
        with gzip.open(path, 'rb') as fd:
            return pickle.load(fd)

    signature_dict = None
    if method == "none":
        # We don't use signature
        logging.warning("We are not using any signature")

    elif method == "full" or method == "complement":

        signature_dict_precomputed = mutations_file + "_signature_full_"+classifier+".pickle.gz"
        if exists(signature_dict_precomputed):
            logging.info("Using precomputed signatures")
            with gzip.open(signature_dict_precomputed, 'rb') as fd:
                signature_dict = pickle.load(fd)
        else:
            logging.info("Computing full global signatures")
            signature_dict = compute_signature(signature_function, classifier, blacklist)
            if save_pickle:
                try:
                    # Try to store as precomputed
                    with gzip.open(signature_dict_precomputed, 'wb') as fd:
                        pickle.dump(signature_dict, fd)
                except OSError:
                    logging.debug("Imposible to write precomputed full signature here: {}".format(signature_dict_precomputed))

        if method == "complement":
            signature_dict = collapse_complementaries(signature_dict)

    elif method == "bysample":
        signature_dict_precomputed = mutations_file + "_signature_bysample.pickle.gz"
        if exists(signature_dict_precomputed):
            logging.info("Using precomputed per sample signatures")
            with gzip.open(signature_dict_precomputed, 'rb') as fd:
                signature_dict = pickle.load(fd)
        else:
            logging.info("Computing signatures per sample")
            signature_dict = compute_signature(signature_function, classifier, blacklist, collapse=True)
            if save_pickle:
                try:
                    # Try to store as precomputed
                    with gzip.open(signature_dict_precomputed, 'wb') as fd:
                        pickle.dump(signature_dict, fd)
                except OSError:
                    logging.debug("Imposible to write precomputed bysample signatures here: {}".format(signature_dict_precomputed))

    elif method == "file":
        if not os.path.exists(path):
            logging.error("Signature file {} not found.".format(path))
            return -1
        else:
            logging.info("Loading signatures")
            signature_probabilities = pd.read_csv(path, sep='\t')
            signature_probabilities.set_index([column_ref, column_alt], inplace=True)
            signature_dict = {classifier: signature_probabilities.to_dict()[column_probability]}
    return signature_dict


def yield_mutations(mutations):
    """
    Yields one mutation each time from
    a list of mutations

    Args:
        mutations (dict): :ref:`mutations <mutations dict>`

    Yields:
        Mutation

    """
    for elem, mutations_list in mutations.items():
        for mutation in mutations_list:
            yield mutation
