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
        'ALT':        {'reader': 'str(x).upper()', 'validator': 'match("^[ACTG-]*$",x) and r[2]!=x'},
        'TYPE':       {'nullable': 'True', 'validator': 'x in ["subs", "indel"]'},
        'SAMPLE':     {'reader': 'str(x)'},
        'SIGNATURE':  {'reader': 'str(x)'}
    }
}

INDELS_HEADER = ["CHROM", "POS", "REF", "ALT", "RawScore", "PHRED"]
INDELS_SCHEMA = {
    'fields': {
        'CHROM': {'reader': 'str(x)', 'validator': "x in ([str(c) for c in range(1,23)] + ['X', 'Y'])"},
        'POS': {'reader': 'int(x)', 'validator': 'x > 0'},
        'REF': {'reader': 'str(x).upper()', 'validator': 'match("^[ACTG-]*$",x)'},
        'ALT': {'reader': 'str(x).upper()', 'validator': 'match("^[ACTG-]*$",x)'},
        'RawScore': {'reader': 'float(x)'},
        'PHRED': {'reader': 'float(x)'}
    }
}


def load_mutations(file, signature=None, show_warnings=True, blacklist=None):

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
            regions_tree[r['chrom']][r['start']:(r['stop']+1)] = (r['feature'], r['segment'])

    logging.info("[{} of {}]".format(i+1, len(regions)))
    return regions_tree


def load_indels_dict(indels_file):

    if indels_file is None:
        return None

    indels = {}
    with itab.DictReader(indels_file, schema=INDELS_SCHEMA) as reader:
        all_errors = []
        for r, errors in reader:

            # Report errors and continue
            if len(errors) > 0:
                all_errors += errors
                continue

            indel_key = "{}|{}|{}|{}".format(r['CHROM'], r['POS'], r['REF'], r['ALT'])
            indels[indel_key] = r['PHRED']

        if len(all_errors) > 0:
            logging.warning("There are {} errors at {}. {}".format(
                len(all_errors), os.path.basename(indels_file),
                " I show you only the ten first errors." if len(all_errors) > 10 else ""
            ))
            for e in all_errors[:10]:
                logging.warning(e)

    return indels


def load_variants_dict(variants_file, regions, indels=None, signature_name='none', blacklist=None):

    if variants_file.endswith(".pickle.gz"):
        with gzip.open(variants_file, 'rb') as fd:
            return pickle.load(fd)

    # Build regions tree
    logging.info("Building regions tree")
    regions_tree = build_regions_tree(regions)

    # Load mutations
    variants_dict = defaultdict(list)

    # Check the file format
    indels_skip = []
    logging.info("Mapping mutations")
    i = 0
    show_small_progress_at = 100000
    show_big_progress_at = 1000000
    for i, r in enumerate(load_mutations(variants_file, signature=signature_name, blacklist=blacklist), start=1):

        if r['CHROMOSOME'] not in regions_tree:
            continue

        if i % show_small_progress_at == 0:
            print('*', end='', flush=True)

        if i % show_big_progress_at == 0:
            print(' [{} muts]'.format(i), flush=True)

        position = int(r['POSITION'])
        intervals = regions_tree[r['CHROMOSOME']][position]

        if len(intervals) > 0:
            if r['TYPE'] == 'indel':
                if indels is None:
                    continue
                else:
                    score = indels.get("{}|{}|{}|{}".format(r['CHROMOSOME'], r['POSITION'], r['REF'], r['ALT']), None)
                    if score is None:
                        # Skip indels without score (may be are very long indels)
                        indels_skip.append("At {}:{} of length {}.".format(r['CHROMOSOME'], r['POSITION'], max(len(r['REF']), len(r['ALT']))))
                        continue
            else:
                score = None

        for interval in intervals:
            feature, segment = interval.data
            variants_dict[feature].append({
                'CHROMOSOME': r['CHROMOSOME'],
                'POSITION': position,
                'SAMPLE': r['SAMPLE'],
                'TYPE': r['TYPE'],
                'REF': r['REF'],
                'ALT': r['ALT'],
                'SCORE': score,
                'SIGNATURE': r['SIGNATURE'],
                'SEGMENT': segment
            })

    if i > show_small_progress_at:
        print('{} [{} muts]'.format(' '*(((show_big_progress_at-(i % show_big_progress_at)) // show_small_progress_at)+1), i), flush=True)

    if len(indels_skip) > 0:
        logging.warning("There are {} indels without score. {}".format(
            len(indels_skip),
            " I show you only the ten first." if len(indels_skip) > 10 else ""
        ))
        for e in indels_skip[:10]:
            logging.warning(e)

    return variants_dict
