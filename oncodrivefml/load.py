from collections import defaultdict
import gzip
import logging
import os
import pickle
from intervaltree import IntervalTree
import itab

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
        'ALT':        {'reader': 'str(x).upper()', 'validator': 'match("^[ACTG-]*$",x)'},
        'TYPE':       {'nullable': 'True', 'validator': 'x in ["subs", "indel"]'},
        'SAMPLE':     {'reader': 'str(x)'},
        'SIGNATURE':  {'reader': 'str(x)'}
    }
}


def load_mutations(file, signature=None, show_warnings=True):
    reader = itab.DictReader(file, header=MUTATIONS_HEADER, schema=MUTATIONS_SCHEMA)
    all_errors = []
    for ix, (row, errors) in enumerate(reader, start=1):
        if len(errors) > 0:
            if reader.line_num == 1:
                # Most probable this is a file with a header
                continue
            all_errors += errors
            continue

        if row.get('TYPE', None) is None:
            if '-' in row['REF'] or '-' in row['ALT'] or len(row['REF']) > 1 or len(row['ALT']) > 1:
                row['TYPE'] = 'indel'
            else:
                row['TYPE'] = 'subs'

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
    regions_tree = defaultdict(IntervalTree)
    for i, (k, allr) in enumerate(regions.items()):

        if i % 7332 == 0:
            logging.info("[{} of {}]".format(i+1, len(regions)))

        for r in allr:
            regions_tree[r['chrom']][r['start']:r['stop']] = (r['feature'], r['segment'])

    logging.info("[{} of {}]".format(i+1, len(regions)))
    return regions_tree


def load_variants_dict(variants_file, regions, signature_name='none'):

    if variants_file.endswith(".pickle.gz"):
        with gzip.open(variants_file, 'rb') as fd:
            return pickle.load(fd)

    # Build regions tree
    regions_tree = build_regions_tree(regions)

    # Load mutations
    variants_dict = defaultdict(list)

    # Check the file format
    for r in load_mutations(variants_file, signature=signature_name):

        if r['CHROMOSOME'] not in regions_tree:
            continue

        position = int(r['POSITION'])
        for interval in regions_tree[r['CHROMOSOME']][position]:
            feature, segment = interval.data
            variants_dict[feature].append({
                'CHROMOSOME': r['CHROMOSOME'],
                'POSITION': position,
                'SAMPLE': r['SAMPLE'],
                'TYPE': r['TYPE'],
                'REF': r['REF'],
                'ALT': r['ALT'],
                'SIGNATURE': r['SIGNATURE'],
                'SEGMENT': segment
            })

    return variants_dict
