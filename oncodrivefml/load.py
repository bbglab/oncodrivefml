"""
This module contains the methods used to
load and parse the input files: elements and mutations

.. _elements dict:

elements (:obj:`dict`)
    contains all the segments related to one element. The information is taken from
    the :file:`elements_file`.
    Basic structure:

    .. code-block:: python

        { element_id:
            [
                {
                'CHROMOSOME': chromosome,
                'START': start_position_of_the_segment,
                'END': end_position_of_the_segment,
                'STRAND': strand (+ -> positive | - -> negative)
                'ELEMENT': element_id,
                'SEGMENT': segment_id,
                'SYMBOL': symbol_id
                }
            ]
        }


.. _mutations dict:

mutations (:obj:`dict`)
    contains all the mutations for each element. Most of the information is taken from
    the mutations_file but the *element_id* and the *segment* that are taken from the **elements**.
    More information is added during the execution.
    Basic structure:

    .. code-block:: python

        { element_id:
            [
                {
                'CHROMOSOME': chromosome,
                'POSITION': position_where_the_mutation_occurs,
                'REF': reference_sequence,
                'ALT': alteration_sequence,
                'SAMPLE': sample_id,
                'ALT_TYPE': type_of_the_mutation,
                'CANCER_TYPE': group to which the mutation belongs to,
                'SIGNATURE': a different grouping category,
                }
            ]
        }

.. _mutations data dict:

mutations_data (:obj:`dict`)
    contains the `mutations dict`_ and some metadata information about the mutations.
    Currently, the number of substitutions and indels.
    Basic structure:

    .. code-block:: python

        {
            'data':
                {
                    `mutations dict`_
                },
            'metadata':
                {
                    'snp': amount of SNP mutations
                    'mnp': amount of MNP mutations
                    'mnp_length': total length of the MNP mutations
                    'indel': amount of indels
                }
        }

"""

import gzip
import logging
import pickle

from bgcache import bgcache
from bgparsers import readers
from collections import defaultdict
from intervaltree import IntervalTree

from oncodrivefml import __logger_name__

logger = logging.getLogger(__logger_name__)


def mutations(file, blacklist=None, metadata_dict=None, indels_max_size=None):
    """
    Parsed the mutations file

    Args:
        file: mutations file (see :class:`~oncodrivefml.main.OncodriveFML`)
        metadata_dict (dict): dict that the function will fill with useful information
        blacklist (optional): file with blacklisted samples (see :class:`~oncodrivefml.main.OncodriveFML`).
            Defaults to None.
        indels_max_size (int, optional): max size of indels. Indels with logner sizes will be discarded.

    Yields:
        One line from the mutations file as a dictionary. Each of the inner elements of
        :ref:`mutations <mutations dict>`

    """

    # Set of samples to blacklist
    samples_blacklisted = set([s.strip() for s in open(blacklist).readlines()]) if blacklist is not None else set()

    snp = 0
    indel = 0
    mnp = 0
    mnp_length = 0

    for row in readers.variants(file, extra=['CANCER_TYPE', 'SIGNATURE'], required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']):
        if row['SAMPLE'] in samples_blacklisted:
            continue

        if row['ALT_TYPE'] == 'snp':
            snp += 1
        elif row['ALT_TYPE'] == 'mnp':
            mnp += 1
            mnp_length += len(row['REF'])
        elif row['ALT_TYPE'] == 'indel':
            # very long indels are discarded
            if indels_max_size and max(len(row['REF']), len(row['ALT'])) > indels_max_size:
                continue
            indel += 1

        yield row

    if metadata_dict is not None:
        metadata_dict['snp'] = snp
        metadata_dict['indel'] = indel
        metadata_dict['mnp'] = mnp
        metadata_dict['mnp_length'] = mnp_length


def snp(file, blacklist=None):
    """Load only SNP"""
    for row in mutations(file, blacklist=blacklist):
        if row['ALT_TYPE'] == 'snp':
            yield row


def build_regions_tree(regions):
    """
    Generates a binary tree with the intervals of the regions

    Args:
        regions (dict): segments grouped by :ref:`elements <elements dict>`.

    Returns:
        dict of :obj:`~intervaltree.IntervalTree`: for each chromosome, it get one :obj:`~intervaltree.IntervalTree` which
        is a binary tree. The leafs are intervals [low limit, high limit) and the value associated with each interval
        is the :obj:`tuple` (element, segment).
        It can be interpreted as:

        .. code-block:: python

            { chromosome:
                (start_position, end_position +1): (element, segment)
            }

    """
    logger.info('Building regions tree')
    regions_tree = {}
    i = 0
    for i, (k, allr) in enumerate(regions.items()):

        if i % 7332 == 0:
            logger.info("[%d of %d]", i+1, len(regions))

        for r in allr:
            tree = regions_tree.get(r['CHROMOSOME'], IntervalTree())
            tree[r['START']:(r['END']+1)] = r['ELEMENT']
            regions_tree[r['CHROMOSOME']] = tree

    logger.info("[%d of %d]", i+1, len(regions))
    return regions_tree


@bgcache
def elements_tree(elements_file):
    elements = readers.elements_dict(elements_file, required=['CHROMOSOME', 'START', 'END', 'ELEMENT'])
    return build_regions_tree(elements)


def mutations_and_elements(variants_file, elements_file, blacklist=None, indels_max_size=None):
    """
    From the elements and variants file, get dictionaries with the segments grouped by element ID and the
    mutations grouped in the same way, as well as some information related to the mutations.

    Args:
        variants_file: mutations file (see :class:`~oncodrivefml.main.OncodriveFML`)
        elements_file: elements file (see :class:`~oncodrivefml.main.OncodriveFML`)
        blacklist (optional): file with blacklisted samples (see :class:`~oncodrivefml.main.OncodriveFML`). Defaults to None.
           If the blacklist option is passed, the mutations are not loaded from a pickle file.
        indels_max_size (int, optional): max size of indels. Indels with logner sizes will be discarded.

    Returns:
        tuple: mutations and elements

        Elements: `elements dict`_

        Mutations: `mutations data dict`_


    The process is done in 3 steps:
       1. :meth:`load_regions`
       #. :meth:`build_regions_tree`.
       #. each mutation (:meth:`mutations`) is associated with the right
          element ID

    """
    # Load elements file
    # TODO add extra and required
    elements = readers.elements_dict(elements_file, required=['CHROMOSOME', 'START', 'END', 'ELEMENT'],
                                     extra=['SEGMENT', 'SYMBOL', 'STRAND'])

    # If the input file is a pickle file do nothing
    if variants_file.endswith(".pickle.gz"):
        with gzip.open(variants_file, 'rb') as fd:
            return pickle.load(fd), elements

    # Loading elements tree
    elements_tree_ = elements_tree(elements_file)

    # Mapping mutations
    variants_dict = defaultdict(list)
    variants_metadata_dict = {}
    logger.info("Mapping mutations")
    i = 0
    show_small_progress_at = 100000
    show_big_progress_at = 1000000
    indels_mapped_multiple_of_3 = 0
    snp_mapped = 0
    mnp_mapped = 0
    indels_mapped = 0
    for i, r in enumerate(mutations(variants_file, metadata_dict=variants_metadata_dict, blacklist=blacklist, indels_max_size=indels_max_size)):

        if r['CHROMOSOME'] not in elements_tree_:
            continue

        if i % show_small_progress_at == 0:
            print('*', end='', flush=True)

        if i % show_big_progress_at == 0:
            print(' [{} muts]'.format(i), flush=True)

        # Get the interval that include that position in the same chromosome
        intervals = elements_tree_[r['CHROMOSOME']][r['POSITION']]

        for interval in intervals:
            element = interval.data
            variants_dict[element].append(r)

        if intervals:
            if r['ALT_TYPE'] == 'snp':
                snp_mapped += 1
            elif r['ALT_TYPE'] == 'mnp':
                mnp_mapped += 1
            else:
                indels_mapped += 1
                if max(len(r['REF']), len(r['ALT'])) % 3 == 0:
                    indels_mapped_multiple_of_3 += 1

    if i > show_small_progress_at:
        print('{} [{} muts]'.format(' '*(((show_big_progress_at-(i % show_big_progress_at)) // show_small_progress_at)+1), i), flush=True)

    variants_metadata_dict['snp_mapped'] = snp_mapped
    variants_metadata_dict['mnp_mapped'] = mnp_mapped
    variants_metadata_dict['indels_mapped'] = indels_mapped
    variants_metadata_dict['indels_mapped_multiple_of_3'] = indels_mapped_multiple_of_3
    mutations_data_dict = {'data': variants_dict, 'metadata': variants_metadata_dict}

    return mutations_data_dict, elements
