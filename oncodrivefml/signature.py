"""
This module contains information related with the signature.

The signature is a way of assigning probabilities to certain mutations that have some
relation amongst them (e.g. cancer type, sample...).

This relation is identified by the **signature_id**.

The ``classifier`` parameter in the :ref:`configuration <project configuration>` of the signature
specifies which column of the mutations file (:data:`~oncodrivefml.load.MUTATIONS_HEADER`) is used as
the identifier for the different signature groups.
If the column does not exist the ``classifier`` itself is used as value for the
*signature_id*.

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

import os
import sys
import gzip
import json
import pickle
import logging
import bgdata
import pandas as pd
from multiprocessing.pool import Pool
from collections import defaultdict, Counter
from bgreference import refseq

from oncodrivefml import __logger_name__
from oncodrivefml.utils import exists_path

logger = logging.getLogger(__logger_name__)

ref_build = 'hg19'
"""
Build of the Reference Genome
"""

__CB = {"A": "T", "T": "A", "G": "C", "C": "G"}


def change_ref_build(build):
    """
    Modify the default build fo the reference genome

    Args:
        build (str): genome reference build

    """
    global ref_build
    ref_build = build
    logger.info('Using %s as reference genome', ref_build.upper())


def get_build():
    return ref_build


def get_ref(chromosome, start, size=1):
    """
    Gets a sequence from the reference genome

    Args:
        chromosome (str): chromosome
        start (int): start position where to look
        size (int): number of bases to retrieve

    Returns:
        str. Sequence from the reference genome

    """
    return refseq(ref_build, chromosome, start, size)


def get_ref_triplet(chromosome, start):
    """

    Args:
        chromosome (str): chromosome identifier
        start (int): starting position

    Returns:
        str: 3 bases from the reference genome

    """
    return get_ref(chromosome, start, size=3)


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


def sum2one_dict(signature_counts):
    """
    Associates to each key (tuple(reference_tripet, altered_triplet)) the value divided by the total amount

    Args:
        signature_counts (dict): pair key-amount {(ref_triplet, alt_triplet): value}

    Returns:
        dict: pair key-(amount/total_amount)

    """
    total = sum([v for v in signature_counts.values()])
    return {k: v/total for k, v in signature_counts.items()}


def compute_signature(signature_function, classifier, collapse=False, include_mnp=False):
    """
    Gets the probability of each substitution that occurs for a certain signature_id.

    Each substitution is identified by the pair (reference_triplet, altered_triplet).

    The signature_id is taken from the mutations field corresponding to the classifier.

    Args:
        signature_function: function that yields one mutation each time
        classifier (str): passed to :func:`~oncodrivefml.load.load_mutations`
            as parameter ``signature_classifier``.
        collapse (bool): consider one substitutions and the complementary one as the same. Defaults to True.
        include_mnp (bool): use MNP mutation in the signature computation or not

    Returns:
        dict: probability of each substitution (measured by the triplets) grouped by the signature_classifier

        .. code-block:: python

            { signature_id:
                {
                    (ref_triplet, alt_triplet): prob
                }
            }

    .. warning::

        Only substitutions (MNP are optional) are taken into account

    """
    total = 0
    mismatches = 0
    signature_count = defaultdict(lambda: defaultdict(int))
    for mut in signature_function():
        pos = mut['POSITION']
        if mut['ALT_TYPE'] == 'snp':
            total += 1
            signature_ref = get_ref_triplet(mut['CHROMOSOME'], mut['POSITION'] - 1)
            signature_alt = signature_ref[0] + mut['ALT'] + signature_ref[2]
            if signature_ref[1] != mut['REF']:
                logger.debug('Discrepancy in substitution at position %d of chr %s', pos, mut['CHROMOSOME'])
                mismatches += 1
                continue

            signature_count[mut.get(classifier, classifier)][(signature_ref, signature_alt)] += 1
        elif include_mnp and mut['ALT_TYPE'] == 'mnp':
            total += 1
            for index, nucleotides in enumerate(zip(mut['REF'], mut['ALT'])):
                ref_nucleotide, alt_nucleotide = nucleotides
                signature_ref = get_ref_triplet(mut['CHROMOSOME'], pos - 1 + index)
                if signature_ref[1] != ref_nucleotide:
                    logger.debug('Discrepancy in MNP at position %d of chr %s', pos + index, mut['CHROMOSOME'])
                    mismatches += 1
                    continue
                signature_alt = signature_ref[0] + alt_nucleotide + signature_ref[2]

                signature_count[mut.get(classifier, classifier)][(signature_ref, signature_alt)] += 1
        else:
            continue

    if mismatches / total > 0.1:
        logger.error('Too many mismatches. You are using %s as reference genome, please check it is right. Program stopped', ref_build)
        sys.exit(-1)
    elif mismatches / total > 0.05:
        logger.warning('There are %d mismatches between your mutations and the reference genome.', mismatches)

    signature = {}
    for k, v in signature_count.items():
        if collapse:
            signature[k] = sum2one_dict(collapse_complementaries(v))
        else:
            signature[k] = sum2one_dict(v)

    return signature


def load_signature(signature_config, signature_function, trinucleotides_counts=None, load_pickle=None, save_pickle=False):
    """
    Computes the probability that certain mutation occurs.

    Args:
        signature_config (dict): information of the signature (see :ref:`configuration <project configuration>`)
        signature_function: function that yields one mutation each time
        trinucleotides_counts (:obj:`dict`, optional): counts of trincleotides used to correct the signature
        load_pickle (:obj:`str`, optional): path to the pickle file
        save_pickle (:obj:`str`, optional): path to pickle file

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
    path = signature_config['path']
    column_ref = signature_config['column_ref']
    column_alt = signature_config['column_alt']
    column_probability = signature_config['column_probability']
    include_mnp = signature_config['include_mnp']

    if path is not None and path.endswith(".pickle.gz"):
        with gzip.open(path, 'rb') as fd:
            return pickle.load(fd)

    signature_dict = None
    if method == "none":
        # We don't use signature
        logger.warning("No signature is being used")

    elif method == "file":
        if not os.path.exists(path):
            logger.error("Signature file %s not found.", path)
            return -1
        else:
            logger.info("Loading signatures")
            signature_probabilities = pd.read_csv(path, sep='\t')
            signature_probabilities.set_index([column_ref, column_alt], inplace=True)
            signature_dict = {classifier: signature_probabilities.to_dict()[column_probability]}
            return signature_dict

    elif method == "full" or method == "complement":
        if method == "complement":
            collapse = True
        else:
            collapse = False

        if exists_path(load_pickle):
            try:
                logger.info("Using precomputed signatures")
                with gzip.open(load_pickle, 'rb') as fd:
                    signature_dict = pickle.load(fd)
            except EOFError:
                logger.error("Loading file %s", load_pickle)
                signature_dict = None

        if signature_dict is None:
            logger.info("Computing signatures")
            signature_dict = compute_signature(signature_function, classifier, collapse, include_mnp)
            if save_pickle is not None:
                try:
                    # Try to store as precomputed
                    with gzip.open(save_pickle, 'wb') as fd:
                        pickle.dump(signature_dict, fd)
                except OSError:
                    logger.debug(
                        "Imposible to write precomputed signature here: {}".format(save_pickle))

        if signature_dict is not None and trinucleotides_counts is not None:  # correct the signature

            if collapse:
                trinucleotides_counts = collapse_complementaries(trinucleotides_counts)

            triplets_probabilities = sum2one_dict(trinucleotides_counts)

            logger.info('Correcting signatures')

            signature_dict = correct_signature_by_triplets_frequencies(signature_dict, triplets_probabilities)

            signature_dict['trinucleotides'] = count_valid_trinucleotides(trinucleotides_counts)

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


def correct_signature_by_triplets_frequencies(signature, triplets_frequencies):
    """
    Normalized de signature by the frequency of the triplets

    Args:
        signature (dict): see :ref:`signature <signature dict>`
        triplets_frequencies (dict): {triplet: frequency}

    Returns:
        dict. Normalized signature

    """
    if signature is None:
        return None
    corrected_signature = {}
    for k,v in signature.items():
        corrected_signature[k] = get_normalized_frequencies(v, triplets_frequencies)

    return corrected_signature


def get_normalized_frequencies(signature, triplets_frequencies):
    """
    Divides the frequency of each triplet alteration by the
    frequency of the reference triplet to get the normalized
    signature

    Args:
        signature (dict): {(ref_triplet, alt_triplet): counts}
        triplets_frequencies (dict): {triplet: frequency}

    Returns:
        dict. Normalized signature

    """
    corrected_signature = {}
    for triplet_pair, frequency in signature.items():
        ref_triplet = triplet_pair[0]
        if ref_triplet not in triplets_frequencies:
            logger.warning('Triplet %s not found', ref_triplet)
        corrected_signature[triplet_pair] = frequency/triplets_frequencies.get(ref_triplet, float("inf"))
    return sum2one_dict(corrected_signature)


def __get_region(region_name):
    if region_name == 'whole_genome' or region_name == 'wgs':
        return 'genome'
    elif region_name == 'whole_exome' or region_name == 'wes' or region_name == 'wxs':
        return 'exome'
    else:
        return region_name


def load_trinucleotides_counts(region):
    """
    Get the trinucleotides counts for a precomputed
    region: whole exome or whole genome

    Args:
        region (str): whole genome or whole exome

    Returns:
        dict. Counts of the different trinucleotides

    """
    region = __get_region(region)
    file = bgdata.get_path('datasets', region+'signature', ref_build)
    with open(file) as fd:
        triplets_counts = json.load(fd)

    return triplets_counts


def triplets(sequence):
    """

    Args:
        sequence (str): sequence of nucleotides

    Yields:
        str. Triplet

    """
    iterator = iter(sequence)

    n1 = next(iterator)
    n2 = next(iterator)

    for n3 in iterator:
        yield n1 + n2 + n3
        n1 = n2
        n2 = n3


def triplet_counter_executor(elements):
    """
    For a list of regions, get all the triplets present
    in all the segments

    Args:
        elements (:obj:`list` of :obj:`list`): list of lists of segments

    Returns:
        :class:`collections.Counter`. Count of each triplet in the regions

    """
    counts = Counter()
    for element in elements:
        for segment in element:
            chrom = segment['CHROMOSOME']
            start = segment['START']
            stop = segment['STOP']
            seq = refseq(ref_build, chrom, start, stop-start+1)
            counts.update(triplets(seq))
    return counts


def chunkizator(iterable, size=1000):
    """
    Creates chunks from an iterable

    Args:
        iterable:
        size (int): elements in the chunk

    Returns:
        list. Chunk

    """
    s = 0
    chunk = []
    for i in iterable:
        if s == size:
            yield chunk
            chunk = []
            s = 0
        chunk.append(i)
        s += 1
    yield chunk


def compute_regions_signature(elements, cores):
    """
    Counts triplets in the elements

    Args:
        elements:
        cores (int): cores to use

    Returns:
        :class:`collections.Counter`. Counts of the triplets in the elements

    """
    counter = Counter()
    pool = Pool(cores)
    for result in pool.imap(triplet_counter_executor, chunkizator(elements, size=500)):
        counter.update(result)
    return counter


def is_valid_trinucleotides(trinucleotide):
    """
    Check if a trinucleotide has a nucleotide distinct than A, C, G, T
    Args:
        trinucleotide (str): triplet

    Returns:
        bool.

    """
    for nucleotide in trinucleotide:
        if nucleotide not in __CB.keys():
            return False
    return True


def count_valid_trinucleotides(trinucleotides_dict):
    """
    Count how many trinucleotides are valid

    Args:
        trinucleotides_dict (dict): trinucleotides counts

    Returns:
        int. Valid trinucleotides

    """
    distinct_trinucleotides = set()
    for k in trinucleotides_dict.keys():
        if is_valid_trinucleotides(k):
            distinct_trinucleotides.add(k)
    return len(k)
