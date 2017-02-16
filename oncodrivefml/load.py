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
                'STOP': end_position_of_the_segment,
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
                'TYPE': type_of_the_mutation,
                'CANCER_TYPE': group to which the mutation belongs to,
                'SEGMENT': segment_id
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
                    'subs': amount_of_subs
                    'indels': amount_o_indels
                }
        }

"""

import gzip
import logging
import os
import pickle
import itab
from os.path import exists
from collections import defaultdict
from bgcache import bgcache
from intervaltree import IntervalTree
from bgparsers import readers

from oncodrivefml.config import remove_extension_and_replace_special_characters as get_name


MUTATIONS_HEADER = ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "TYPE", "CANCER_TYPE"]
"""
Headers of the data expected in the mutations file file (see :class:`~oncodrivefml.main.OncodriveFML`).
"""

MUTATIONS_SCHEMA = {
    'fields': {
        'CHROMOSOME':  {'reader': 'str(x)', 'validator': "x in ([str(c) for c in range(1,23)] + ['X', 'Y'])"},
        'POSITION':    {'reader': 'int(x)', 'validator': 'x > 0'},
        'REF':         {'reader': 'str(x).upper()', 'validator': 'match("^[ACTG-]*$",x)'},
        'ALT':         {'reader': 'str(x).upper()', 'validator': 'match("^[ACTG-]*$",x) and r[2]!=x'},
        'SAMPLE':      {'reader': 'str(x)'},
        'TYPE':        {'nullable': 'True', 'validator': 'x in ["subs", "indel", "mnp"]'},
        'CANCER_TYPE': {'reader': 'str(x)', 'nullable': 'True'},
    }
}


def load_mutations(file, show_warnings=True, blacklist=None, metadata_dict=None):
    """
    Parsed the mutations file

    Args:
        file: mutations file (see :class:`~oncodrivefml.main.OncodriveFML`)
        metadata_dict (dict): dict that the function will fill with useful information
        show_warnings (bool, optional): Defaults to True.
        blacklist (optional): file with blacklisted samples (see :class:`~oncodrivefml.main.OncodriveFML`).
            Defaults to None.

    Yields:
        One line from the mutations file as a dictionary. Each of the inner elements of
        :ref:`mutations <mutations dict>`

    """

    # Set of samples to blacklist
    samples_blacklisted = set([s.strip() for s in open(blacklist).readlines()]) if blacklist is not None else set()

    subs = 0
    indels = 0
    mnp = 0
    mnp_length = 0

    reader = itab.DictReader(file, schema=MUTATIONS_SCHEMA) # TODO add the header and switch SAMPLE and TYPE?
    all_errors = []
    for ix, (row, errors) in enumerate(reader, start=1):
        if len(errors) > 0:
            if reader.line_num == 1:
                # Most probable this is a file with a header
                continue
            all_errors += errors
            continue

        if row.get('SAMPLE', None) in samples_blacklisted:
            continue

        if row.get('TYPE', None) is None:
            if '-' in row['REF'] or '-' in row['ALT'] or len(row['REF']) != len(row['ALT']):
                row['TYPE'] = 'indel'
            elif len(row['REF']) == len(row['ALT']) and len(row['REF']) > 1:
                row['TYPE'] = 'mnp'
            else:
                row['TYPE'] = 'subs'

        if metadata_dict is not None:
            # compute the metatada needed
            if row['TYPE'] == 'indel':
                indels += 1
            elif row['TYPE'] == 'subs':
                subs += 1
            else:
                mnp_length += len(row['REF'])
                mnp += 1

        yield row

    if show_warnings and len(all_errors) > 0:
        logging.warning("There are {} errors at {}. {}".format(
            len(all_errors), os.path.basename(file),
            " I show you only the ten first errors." if len(all_errors) > 10 else ""
        ))
        for e in all_errors[:10]:
            logging.warning(e)

    if metadata_dict is not None:
        metadata_dict['subs'] = subs
        metadata_dict['indels'] = indels
        metadata_dict['mnp'] = mnp
        metadata_dict['mnp_length'] = mnp_length

    reader.fd.close()


def build_regions_tree(regions):
    """
    Generates a binary tree with the intervals of the regions

    Args:
        regions (dict): segments grouped by :ref:`elements <elements dict>`.

    Returns:
        :obj:`dict` of :obj:`IntervalTree`: for each chromosome, it get one :obj:`IntervalTree` which
        is a binary tree. The leafs are intervals [low limit, high limit) and the value associated with each interval
        is the :obj:`tuple` (element, segment).
        It can be interpreted as:

        .. code-block:: python

            { chromosome:
                (start_position, stop_position +1): (element, segment)
            }

    """
    regions_tree = {}
    for i, (k, allr) in enumerate(regions.items()):

        if i % 7332 == 0:
            logging.info("[{} of {}]".format(i+1, len(regions)))

        for r in allr:
            tree = regions_tree.get(r['CHROMOSOME'], IntervalTree())
            tree[r['START']:(r['STOP']+1)] = (r['ELEMENT'], r['SEGMENT'])
            regions_tree[r['CHROMOSOME']] = tree

    logging.info("[{} of {}]".format(i+1, len(regions)))
    return regions_tree


@bgcache
def load_elements_tree(elements_file):
    elements = readers.elements(elements_file)
    return build_regions_tree(elements)


def load_and_map_variants(variants_file, elements_file, blacklist=None, save_pickle=False):
    """
    From the elements and variants file, get dictionaries with the segments grouped by element ID and the
    mutations grouped in the same way, as well as some information related to the mutations.

    Args:
        variants_file: mutations file (see :class:`~oncodrivefml.main.OncodriveFML`)
        elements_file: elements file (see :class:`~oncodrivefml.main.OncodriveFML`)
        blacklist (optional): file with blacklisted samples (see :class:`~oncodrivefml.main.OncodriveFML`). Defaults to None.
           If the blacklist option is passed, the mutations are not loaded from/saved to a pickle file.
        save_pickle (:obj:`bool`, optional): save pickle files

    Returns:
        tuple: mutations and elements

        Elements: `elements dict`_

        Mutations: `mutations data dict`_


    The process is done in 3 steps:
       1. :meth:`load_regions`
       #. :meth:`build_regions_tree`.
       #. each mutation (:meth:`load_mutations`) is associated with the right
          element ID

    """
    # Load elements file
    elements = readers.elements(elements_file)

    # If the input file is a pickle file do nothing
    if variants_file.endswith(".pickle.gz"):
        with gzip.open(variants_file, 'rb') as fd:
            return pickle.load(fd), elements

    # Check if it's already done
    variants_dict_precomputed = variants_file + "_mapping_" + get_name(elements_file) + '.pickle.gz'
    if exists(variants_dict_precomputed) and blacklist is None:
        try:
            logging.info("Using precomputed mutations mapping")
            with gzip.open(variants_dict_precomputed, 'rb') as fd:
                return pickle.load(fd), elements
        except EOFError:
            logging.error("Loading file {}".format(variants_dict_precomputed))

    # Loading elements tree
    elements_tree = load_elements_tree(elements_file)

    # Mapping mutations
    variants_dict = defaultdict(list)
    variants_metadata_dict = {}
    logging.info("Mapping mutations")
    i = 0
    show_small_progress_at = 100000
    show_big_progress_at = 1000000
    for i, r in enumerate(load_mutations(variants_file, metadata_dict=variants_metadata_dict, blacklist=blacklist), start=1):

        if r['CHROMOSOME'] not in elements_tree:
            continue

        if i % show_small_progress_at == 0:
            print('*', end='', flush=True)

        if i % show_big_progress_at == 0:
            print(' [{} muts]'.format(i), flush=True)

        position = int(r['POSITION'])
        # Get the interval that include that position in the same chromosome
        intervals = elements_tree[r['CHROMOSOME']][position]

        for interval in intervals:
            element, segment = interval.data
            mutation = r
            mutation['POSITION'] = position
            mutation['SEGMENT'] = segment
            variants_dict[element].append(mutation)

    if i > show_small_progress_at:
        print('{} [{} muts]'.format(' '*(((show_big_progress_at-(i % show_big_progress_at)) // show_small_progress_at)+1), i), flush=True)

    mutations_data_dict = {'data': variants_dict, 'metadata': variants_metadata_dict}

    if save_pickle:
        # Try to store as precomputed
        try:
            with gzip.open(variants_dict_precomputed, 'wb') as fd:
                pickle.dump(mutations_data_dict, fd)
        except OSError:
            logging.debug("Imposible to write precomputed mutations mapping here: {}".format(variants_dict_precomputed))

    return mutations_data_dict, elements
