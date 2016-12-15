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
                'STRAND': strand (+ -> positive | - -> negative | . -> unknown)
                'FEATURE': element_id,
                'SEGMENT': segment_id
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
                'SAMPLE': sample_id,
                'TYPE': type_of_the_mutation,
                'REF': reference_sequence,
                'ALT': alteration_sequence,
                'SIGNATURE': group to which the mutation belongs to,
                'SEGMENT': segment_id
                }
            ]
        }

"""

import gzip
import logging
import os
import pickle
import itab
from os.path import exists
from collections import defaultdict
from intervaltree import IntervalTree

from oncodrivefml.config import remove_extension_and_replace_special_characters as get_name


REGIONS_HEADER = ['CHROMOSOME', 'START', 'STOP', 'STRAND', 'FEATURE', 'SEGMENT']
"""
Headers of the data expected in the elements file (see :class:`~oncodrivefml.main.OncodriveFML`).
"""

REGIONS_SCHEMA = {
    'fields': {
        'CHROMOSOME': {'reader': 'str(x)', 'validator': "x in ([str(c) for c in range(1,23)] + ['X', 'Y'])"},
        'START': {'reader': 'int(x)', 'validator': 'x > 0'},
        'STOP': {'reader': 'int(x)', 'validator': 'x > 0'},
        'STRAND': {'reader': 'str(x)', 'validator': "x in ['.', '+', '-']"},
        'FEATURE': {'reader': 'str(x)'},
        'SEGMENT': {'reader': 'str(x)', 'nullable': 'True'}
}}


MUTATIONS_HEADER = ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "TYPE", 'CANCER_TYPE']
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
        'CANCER_TYPE': {'reader': 'str(x)', 'nullable': 'True'}
    }
}


def load_mutations(file, show_warnings=True, blacklist=None):
    """
    Parsed the mutations file

    Args:
        file: mutations file (see :class:`~oncodrivefml.main.OncodriveFML`)
        show_warnings (bool, optional): Defaults to True.
        blacklist (optional): file with blacklisted samples (see :class:`~oncodrivefml.main.OncodriveFML`).
            Defaults to None.

    Yields:
        One line from the mutations file as a dictionary. Each of the inner elements of
        :ref:`mutations <mutations dict>`

    """

    # Set of samples to blacklist
    samples_blacklisted = set([s.strip() for s in open(blacklist).readlines()]) if blacklist is not None else set()

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

        yield row

    if show_warnings and len(all_errors) > 0:
        logging.warning("There are {} errors at {}. {}".format(
            len(all_errors), os.path.basename(file),
            " I show you only the ten first errors." if len(all_errors) > 10 else ""
        ))
        for e in all_errors[:10]:
            logging.warning(e)

    reader.fd.close()


def load_regions(file):
    """
    Parse an elements file compliant with :attr:`REGIONS_HEADER`

    Args:
        file: elements file. If 'segment' field is not present, the value of the 'feature'
            is used instead.

    Returns:
        dict: elements where the element_id is the 'feature' (see :ref:`elements <elements dict>`).

    """

    regions = defaultdict(list)
    with itab.DictReader(file, header=REGIONS_HEADER, schema=REGIONS_SCHEMA) as reader:
        all_errors = []
        for r, errors in reader:
            # Report errors and continue
            if len(errors) > 0:
                all_errors += errors
                continue

            # If there are no segments use the feature as randomization segment
            if r['SEGMENT'] is None:
                r['SEGMENT'] = r['FEATURE']

            regions[r['FEATURE']].append(r)

        if len(all_errors) > 0:
            logging.warning("There are {} errors at {}. {}".format(
                len(all_errors), os.path.basename(file),
                " I show you only the ten first errors." if len(all_errors) > 10 else ""
            ))
            for e in all_errors[:10]:
                logging.warning(e)
    logging.info("Regions: {}".format(len(regions)))
    return regions


def build_regions_tree(regions):
    """
    Generates a binary tree with the intervals of the regions

    Args:
        regions (dict): segments grouped by :ref:`elements <elements dict>`.

    Returns:
        :obj:`dict` of :obj:`IntervalTree`: for each chromosome, it get one :obj:`IntervalTree` which
        is a binary tree. The leafs are intervals [low limit, high limit) and the value associated with each interval
        is the :obj:`tuple` (feature, segment).
        It can be interpreted as:

        .. code-block:: python

            { chromosome:
                (start_position, stop_position +1): (feature, segment)
            }

    """
    regions_tree = defaultdict(IntervalTree)
    for i, (k, allr) in enumerate(regions.items()):

        if i % 7332 == 0:
            logging.info("[{} of {}]".format(i+1, len(regions)))

        for r in allr:
            regions_tree[r['CHROMOSOME']][r['START']:(r['STOP']+1)] = (r['FEATURE'], r['SEGMENT'])

    logging.info("[{} of {}]".format(i+1, len(regions)))
    return regions_tree


def load_and_map_variants(variants_file, elements_file, blacklist=None, save_pickle=False):
    """
    From the elements and variants files, get dictionaries with the segments grouped by element ID and the
    mutations grouped in the same way.

    Args:
        variants_file: mutations file (see :class:`~oncodrivefml.main.OncodriveFML`)
        elements_file: elements file (see :class:`~oncodrivefml.main.OncodriveFML`)
        blacklist (optional): file with blacklisted samples (see :class:`~oncodrivefml.main.OncodriveFML`). Defaults to None.
        save_pickle (:obj:`bool`, optional): save pickle files

    Returns:
        tuple: mutations and elements

        Elements: `elements dict`_

        Mutations `mutations dict`_


    The process is done in 3 steps:
       1. :meth:`load_regions`
       #. :meth:`build_regions_tree`.
       #. each mutation (:meth:`load_mutations`) is associated with the right
          element ID

    For each of the steps, a precomputed file is search:
        1.
            .. code-block:: python

                elements_file + "_dict.pickle.gz"

        #. Only if the mutations precomputed file is not found

            .. code-block:: python

                elements_file + "_tree.pickle.gz"

        #.
            .. code-block:: python

                variants_file + "_mapping_" + get_name(elements_file) + '.pickle.gz'

    """
    # Load elements file
    elements = None
    elements_precomputed = elements_file + "_dict.pickle.gz"
    if exists(elements_precomputed):
        try:
            logging.info("Using precomputed elements map")
            with gzip.open(elements_precomputed, 'rb') as fd:
                elements = pickle.load(fd)
        except EOFError:
            logging.error("Loading file {}".format(elements_precomputed))
            elements = None

    if elements is None:
        logging.info("Loading elements")
        elements = load_regions(elements_file)
        if save_pickle:
            # Try to store as precomputed
            try:
                with gzip.open(elements_precomputed, 'wb') as fd:
                    pickle.dump(elements, fd)
            except OSError:
                logging.debug("Imposible to write precomputed elements map here: {}".format(elements_precomputed))

    # If the input file is a pickle file do nothing
    if variants_file.endswith(".pickle.gz"):
        with gzip.open(variants_file, 'rb') as fd:
            return pickle.load(fd), elements

    # Check if it's already done
    variants_dict_precomputed = variants_file + "_mapping_" + get_name(elements_file) + '.pickle.gz'
    if exists(variants_dict_precomputed):
        try:
            logging.info("Using precomputed mutations mapping")
            with gzip.open(variants_dict_precomputed, 'rb') as fd:
                return pickle.load(fd), elements
        except EOFError:
            logging.error("Loading file {}".format(variants_dict_precomputed))

    # Loading elements tree
    elements_tree = None
    elements_tree_precomputed = elements_file + "_tree.pickle.gz"
    if exists(elements_tree_precomputed):
        try:
            logging.info("Using precomputed genomic elements tree")
            with gzip.open(elements_tree_precomputed, 'rb') as fd:
                elements_tree = pickle.load(fd)
        except EOFError:
            logging.error("Loading file {}".format(elements_tree_precomputed))
            elements_tree = None

    if elements_tree is None:
        logging.info("Loading genomic elements tree")
        elements_tree = build_regions_tree(elements)
        if save_pickle:
            # Try to store as precomputed
            try:
                with gzip.open(elements_tree_precomputed, 'wb') as fd:
                    pickle.dump(elements_tree, fd)
            except OSError:
                logging.debug("Imposible to write precomputed genomic elements tree here: {}".format(elements_tree_precomputed))

    # Mapping mutations
    variants_dict = defaultdict(list)
    logging.info("Mapping mutations")
    i = 0
    show_small_progress_at = 100000
    show_big_progress_at = 1000000
    for i, r in enumerate(load_mutations(variants_file, blacklist=blacklist), start=1):

        if r['CHROMOSOME'] not in elements_tree:
            continue

        if i % show_small_progress_at == 0:
            print('*', end='', flush=True)

        if i % show_big_progress_at == 0:
            print(' [{} muts]'.format(i), flush=True)

        position = int(r['POSITION'])
        #Get all intervals that include that position in the same chromosome
        intervals = elements_tree[r['CHROMOSOME']][position]

        for interval in intervals:
            element, segment = interval.data
            mutation = r
            mutation['POSITION'] = position
            mutation['SEGMENT'] = segment
            variants_dict[element].append(mutation)

    if i > show_small_progress_at:
        print('{} [{} muts]'.format(' '*(((show_big_progress_at-(i % show_big_progress_at)) // show_small_progress_at)+1), i), flush=True)

    if save_pickle:
        # Try to store as precomputed
        try:
            with gzip.open(variants_dict_precomputed, 'wb') as fd:
                pickle.dump(variants_dict, fd)
        except OSError:
            logging.debug("Imposible to write precomputed mutations mapping here: {}".format(variants_dict_precomputed))

    return variants_dict, elements



def count_mutations(file, show_warnings=True, blacklist=None):
    """
    Count subs and indels in the mutations file
    using :meth:`load_mutations`.

    Args:
        file:
        show_warnings (bool):
        blacklist:

    Returns:
        tuple. Amount of subs and amount of indels

    """
    # TODO add option to not read the file twice
    logging.info('Counting subs and indels')
    subs = 0
    indels = 0
    for mut in load_mutations(file, blacklist=blacklist):
        if mut['TYPE'] == 'indel':
            indels += 1
        else:
            subs += 1
    return subs, indels
