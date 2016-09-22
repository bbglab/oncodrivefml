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
import functools
import gzip
import json
import logging
import os
import mmap
import pickle
from multiprocessing.pool import Pool

import bgdata
import pandas as pd
from os.path import join, exists
from oncodrivefml.load import load_mutations
from collections import defaultdict, Counter

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


def signature_counts(signature_function, classifier):
    """
    Count the number of mutation of each type (ref_triplet, alt_triplet)
    are found

    Args:
        signature_function: function that yields mutations
        classifier (str): identifier of the signature

    Returns:

    """
    signature_count = defaultdict(lambda: defaultdict(int))
    for mut in signature_function():
        if mut['TYPE'] != 'subs':
            continue

        signature_ref = get_ref_triplet(mut['CHROMOSOME'], mut['POSITION'] - 1)
        signature_alt = signature_ref[0] + mut['ALT'] + signature_ref[2]

        signature_count[mut.get(classifier, classifier)][(signature_ref, signature_alt)] += 1
    return signature_count

def compute_signature(signature_function, classifier, collapse=False):
    """
    Gets the probability of each substitution that occurs for a certain signature_id.

    Each substitution is identified by the pair (reference_triplet, altered_triplet).

    The signature_id is taken from the mutations ``SIGNATURE`` field.

    Args:
        signature_function: function that yields one mutation each time
        classifier (str): passed to :func:`oncodrivefml.load.load_mutations`
            as parameter ``signature_classifier``.
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
    signature_count = signature_counts(signature_function, classifier)

    signature = {}
    for k, v in signature_count.items():
        if collapse:
            signature[k] = sum2one_dict(collapse_complementaries(v))
        else:
            signature[k] = sum2one_dict(v)

    return signature


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


def get_signature_function(only_mapped_mutations, mutations, mutations_file, blacklist=None):
    """
    Computes the corresponding signature function depending on the *only_mapped_mutations* flag

    Args:
        only_mapped_mutations (bool): flag indicating if only the mutation that fall into regions are
            used to compute the signature or the whole dataset
        mutations (dict): :ref:`mutations <mutations dict>`
        mutations_file (str): file path
        blacklist (str): file with blacklisted samples

    Returns:
        function. :func:`yield_yield_mutations` if only_mapped_mutations,
            otherwise :func:`oncodrivefml.load.load_mutations`

    """
    if only_mapped_mutations:
        signature_function = functools.partial(yield_mutations, mutations)
    else:
        signature_function = functools.partial(load_mutations, mutations_file, show_warnings=False,
                                               blacklist=blacklist)
    return signature_function


def load_signature(signature_config, regions_file, mutations_file, regions, mutations, save_pickle=False, blacklist=None, cores=None):
    """
    Computes the probability that certain mutation occurs.

    Args:
        signature_config (dict): information of the signature (see :ref:`configuration <project configuration>`)
        regions_file (str): regions file
        mutations_file (str): mutations file
        regions (dict): see :ref:`elements dict`
        mutations (dict): see :ref:`mutations dict`
        save_pickle (:obj:`bool`, optional): save pickle files
        blacklist (str): blacklist used with the mutations file
        cores (int): cores to use

    Returns:
        dict. Signature

    """
    method = signature_config['method']
    classifier = signature_config['classifier']
    limit_signature_to_mapped_mutations = signature_config.get('use_only_mapped_mutations', False)
    correct_signature_by_sites = signature_config.get('correct_signature_by_sites', True)

    if method == "complement":
        collapse = True
    else:
        collapse = False

    signature_dict = None
    if classifier == 'none' or classifier == 'SAMPLE':
        signature_dict_precomputed = mutations_file + '_signature_' + method + '_' + classifier + ".pickle.gz"
    else:
        signature_dict_precomputed = None
        save_pickle = False

    if signature_dict_precomputed is not None and exists(signature_dict_precomputed):
        logging.info("Using precomputed signatures")
        with gzip.open(signature_dict_precomputed, 'rb') as fd:
            signature_dict = pickle.load(fd)
    else:
        logging.info("Computing signatures")

        signature_function = get_signature_function(limit_signature_to_mapped_mutations, mutations, mutations_file, blacklist)

        signature_dict = compute_signature(signature_function, classifier, collapse)
        if save_pickle:
            try:
                # Try to store as precomputed
                with gzip.open(signature_dict_precomputed, 'wb') as fd:
                    pickle.dump(signature_dict, fd)
            except OSError:
                logging.debug(
                    "Imposible to write precomputed signature here: {}".format(signature_dict_precomputed))

    # TODO add checks for the signature

    if signature_dict is not None and correct_signature_by_sites:  # correct the signature

        triplets_probabilities = load_regions_signature(limit_signature_to_mapped_mutations, collapse, regions,
                                                        regions_file, save_pickle, cores)

        logging.info('Correcting signatures')

        signature_dict = correct_signature_by_triplets_frequencies(signature_dict, triplets_probabilities)

    return signature_dict


def get_signature(signature_config, mutations_file, mutations, regions_file, regions, save_pickle=False, blacklist=None, cores=None):
    """
    Computes the probability that certain mutation occurs.

    Args:
        signature_config (dict): information of the signature (see :ref:`configuration <project configuration>`)
        mutations_file (str): mutations file
        mutations (dict): see :ref:`mutations dict`
        regions_file (str): regions file
        regions (dict): see :ref:`elements dict`
        save_pickle (:obj:`bool`, optional): save pickle files
        blacklist (str): blacklist used with the mutations file
        cores (int): cores to use

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
    path = signature_config.get('path', None)

    if path is not None and path.endswith(".pickle.gz"):
        with gzip.open(path, 'rb') as fd:
            return pickle.load(fd)

    signature_dict = None
    if method == "none":
        # We don't use signature
        logging.warning("We are not using any signature")

    elif method == "file":  # read the signature from a file
        if not os.path.exists(path):
            logging.error("Signature file {} not found.".format(path))
            return -1
        else:
            classifier = signature_config['classifier']
            column_ref = signature_config.get('column_ref', None)
            column_alt = signature_config.get('column_alt', None)
            column_probability = signature_config.get('column_probability', None)

            logging.info("Loading signatures")
            signature_probabilities = pd.read_csv(path, sep='\t')
            signature_probabilities.set_index([column_ref, column_alt], inplace=True)
            signature_dict = {classifier: signature_probabilities.to_dict()[column_probability]}

    elif method == "full" or method == "complement":
        signature_dict = load_signature(signature_config, regions_file, mutations_file, regions, mutations, save_pickle, blacklist, cores)

    return signature_dict


def correct_signature_by_triplets_frequencies(signature, triplets_frequencies):
    """
    Normalized de signature by the frequency of the triplets

    Args:
        signature (dict): see :ref:`signature dict`
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


def collapse_complementary_counts(triplets_counts):
    """
    Collapse complementaries

    Args:
        triplets_counts (dict): {triplet: counts}

    Returns:
        dict. Complementary sequences values added together

    """
    comp_sig = defaultdict(int)
    for k, v in triplets_counts.items():
        comp_sig[k] += v
        comp_k = complementary_sequence(k)
        comp_sig[comp_k] += v
    return comp_sig


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
        corrected_signature[triplet_pair] = frequency/triplets_frequencies[ref_triplet]
    return sum2one_dict(corrected_signature)


def triplets_of_sequence(sequence):
    """
    Yields each triplet from a sequence of nucleotides

    Args:
        sequence (str): sequence of nucleotides

    Yields:
        str. Triplet

    """
    iterator = iter(sequence)

    n1 = next(iterator)
    n2 = next(iterator)

    for n3 in iterator:
        yield n1+n2+n3
        n1 = n2
        n2 = n3


def triplet_counter_executor(elements_segments):
    """
    For a list of regions, get all the triplets present
    in all the segments

    Args:
        elements_segments (:obj:`list` of :obj:`list`): list of lists of segments

    Returns:
        :class:`collections.Counter`. Count of each triplet
        in the regions

    """
    counts = Counter()
    for segments in elements_segments:
        for segment in segments:
            chrom = segment['chrom']
            start = segment['start']
            stop = segment['stop']
            seq = get_ref(chrom, start, stop-start+1)
            counts.update(triplets_of_sequence(seq))
    return counts

def triplets_counter(elements, cores=None):
    """
    Counts triplets in the elements

    Args:
        elements (ditc): see :ref:`elements dict`
        cores (int): cores to use

    Returns:
        :class:`collections.Counter`. Counts of the triplets in the elements

    """
    logging.info('Computing the signature of the region')
    if cores is None:
        cores = os.cpu_count()
    counter = Counter()
    pool = Pool(cores)
    for result in pool.imap(triplet_counter_executor, chunkizator(elements.values(), size=500)):
        counter.update(result)
    return counter


def load_regions_signature(only_mapped_mutations, collapse, regions, regions_file, save_pickle=False, cores=None):
    """
    Load the frequency of the triplets

    Args:
        only_mapped_mutations (bool): if we use only the mapped mutations, we must use the
            regions to get the counts of the triplets. If not, we use the whole genome counts.
        collapse (bool): collapse complementaries
        regions (dict): see :ref:`elements dict`
        regions_file (str): path to elements file
        save_pickle (bool): save intermediate pickle files
        cores (int): cores to use

    Returns:
        dict. {triplet: prob}

    """
    if only_mapped_mutations:  # use only mapped regions
        precomputed_regions_counts = regions_file + '_signature.json'
        if exists(precomputed_regions_counts):
            with open(precomputed_regions_counts) as fd:
                triplets_counts = json.load(fd)
        else:
            triplets_counts = triplets_counter(regions, cores)
            if save_pickle:
                with open(precomputed_regions_counts, 'w') as fd:
                    json.dump(triplets_counts, fd)

    else:  # use whole genome
        file = bgdata.get_path('datasets', 'genomesignature', 'hg19')
        with open(file) as fd:
            triplets_counts = json.load(fd)

    if collapse:
        triplets_counts = collapse_complementary_counts(triplets_counts)

    return sum2one_dict(triplets_counts)



def chunkizator(iterable, size=1000):
    """
    Creates chunks from an iterable

    Args:
        iterable:
        size (int): elements in the chunk

    Yields:
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
