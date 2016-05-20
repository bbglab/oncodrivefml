from collections import defaultdict
import gzip
import logging
import os
import pickle
from intervaltree import IntervalTree
import itab

from os.path import exists

from oncodrivefml.config import file_name

REGIONS_HEADER = ['chrom', 'start', 'stop', 'feature', 'segment']
REGIONS_SCHEMA = {
    'fields': {
        'chrom': {'reader': 'str(x)', 'validator': "x in ([str(c) for c in range(1,23)] + ['X', 'Y'])"},
        'start': {'reader': 'int(x)', 'validator': 'x > 0'},
        'stop': {'reader': 'int(x)', 'validator': 'x > 0'},
        'feature': {'reader': 'str(x)'},
        'segment': {'reader': 'str(x)', 'nullable': 'True'}
}}


MUTATIONS_HEADER = ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "TYPE", "SIGNATURE"]
MUTATIONS_SCHEMA = {
    'fields': {
        'CHROMOSOME': {'reader': 'str(x)', 'validator': "x in ([str(c) for c in range(1,23)] + ['X', 'Y'])"},
        'POSITION':   {'reader': 'int(x)', 'validator': 'x > 0'},
        'REF':        {'reader': 'str(x).upper()', 'validator': 'match("^[ACTG-]*$",x)'},
        'ALT':        {'reader': 'str(x).upper()', 'validator': 'match("^[ACTG-]*$",x) and r[2]!=x'},
        'TYPE':       {'nullable': 'True', 'validator': 'x in ["subs", "indel"]'},
        'SAMPLE':     {'reader': 'str(x)'},
        'SIGNATURE':  {'reader': 'str(x)'}
    }
}


def load_mutations(file, signature=None, show_warnings=True, blacklist=None):
    """
    Yields one line from the mutations file as a dictionary
    :param file: mutations file with format: ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "TYPE", "SIGNATURE"]
    :param signature: only affects if it is bysample (in this case the
    signatures of the sample [local] if modified to be the sample itself) or if
    the signatures of the sample [local] is not present (in this case,
    this value is used)
    :param show_warnings:
    :param blacklist: mutations that are omitted from the list
    :return:
    """

    # Set of samples to blacklist
    samples_blacklisted = set([s.strip() for s in open(blacklist).readlines()]) if blacklist is not None else set()

    reader = itab.DictReader(file, header=MUTATIONS_HEADER, schema=MUTATIONS_SCHEMA)
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
            if '-' in row['REF'] or '-' in row['ALT'] or len(row['REF']) > 1 or len(row['ALT']) > 1:
                row['TYPE'] = 'indel'
            else:
                row['TYPE'] = 'subs'

        if signature == 'bysample':
            row['SIGNATURE'] = row['SAMPLE']
        else:
            if row.get('SIGNATURE', None) is None:
                row['SIGNATURE'] = signature

            if row.get('CANCER_TYPE', None) is not None:
                row['SIGNATURE'] = row['CANCER_TYPE']

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

    :param file: mutation file with format ['chrom', 'start', 'stop', 'feature', 'segment']
    (see :ref: REGIONS_HEADER)
    :return: {feature: [ {chromosome: , start: , stop: , feature:, segment*:
    } ] }
    *: if the field is not present, feature is used for it
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
            if r['segment'] is None:
                r['segment'] = r['feature']

            regions[r['feature']].append(r)

        if len(all_errors) > 0:
            logging.warning("There are {} errors at {}. {}".format(
                len(all_errors), os.path.basename(file),
                " I show you only the ten first errors." if len(all_errors) > 10 else ""
            ))
            for e in all_errors[:10]:
                logging.warning(e)
    return regions


def build_regions_tree(regions):
    """

    :param regions:
    :return: { chromosome: IntrevalTree }
    IntervalTree: binary tree with an interval [low limit, high limit) and a
    value for that interval (can be anything). E.g. [5:(10)]='a string as value'
    In our case [start : (stop+1)] = (feature, segment)

    """
    regions_tree = defaultdict(IntervalTree)
    for i, (k, allr) in enumerate(regions.items()):

        if i % 7332 == 0:
            logging.info("[{} of {}]".format(i+1, len(regions)))

        for r in allr:
            regions_tree[r['chrom']][r['start']:(r['stop']+1)] = (r['feature'], r['segment'])

    logging.info("[{} of {}]".format(i+1, len(regions)))
    return regions_tree


def load_and_map_variants(variants_file, elements_file, signature_name='none', blacklist=None):
    """

    :param variants_file:
    :param elements_file:
    :param signature_name:
    :param blacklist:
    :return: variants_dict, elements
    variants_dict:
    { element* : [ {chromosome: , positon: , sample: , type: , ref: ,
    alt: , signatures: , segment*: } ] }
    *: this value comes from the elements file
    elements:
    {feature: [ {chromosome: , start: , stop: , feature: , segment: } ] }

    Checks if the pickle.gz already exists for mutations and elements files,
    and uses those if possible.
    If not, it loads the elements (see :ref: load_regions), generates the
    elements-tree (see :ref: build_regions_tree) and combines it with the
    variants_file (see :ref: load_mutations) to generate the variants_dict.
    Blacklist and signatures are used when loading the mutations file.
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
    variants_dict_precomputed = variants_file + "_mapping_" + file_name(elements_file) + '.pickle.gz'
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
    for i, r in enumerate(load_mutations(variants_file, signature=signature_name, blacklist=blacklist), start=1):

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
            variants_dict[element].append({
                'CHROMOSOME': r['CHROMOSOME'],
                'POSITION': position,
                'SAMPLE': r['SAMPLE'],
                'TYPE': r['TYPE'],
                'REF': r['REF'],
                'ALT': r['ALT'],
                'SIGNATURE': r['SIGNATURE'],
                'SEGMENT': segment
            })

    if i > show_small_progress_at:
        print('{} [{} muts]'.format(' '*(((show_big_progress_at-(i % show_big_progress_at)) // show_small_progress_at)+1), i), flush=True)

    # Try to store as precomputed
    try:
        with gzip.open(variants_dict_precomputed, 'wb') as fd:
            pickle.dump(variants_dict, fd)
    except OSError:
        logging.debug("Imposible to write precomputed mutations mapping here: {}".format(variants_dict_precomputed))

    return variants_dict, elements
