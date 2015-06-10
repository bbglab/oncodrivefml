import csv
import logging
import os
import stat
import mmap
import itab
from collections import Counter, defaultdict
import subprocess

import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt
from intervaltree import IntervalTree


#TODO This must be configurable
TABIX_LIST = [
    "/usr/bin/tabix",
    os.path.expanduser('~/bin/tabix'),
    os.path.expanduser(os.path.expandvars("$TABIX_PATH")),
    "/soft/bio/sequence/tabix-0.2.3/tabix"
]

TABIX = None
for t in TABIX_LIST:
    if os.path.exists(t):
        TABIX = t
        break


HG19_LIST = [
    "/projects_bg/bg/soft/intogen_home/gencluster/software/mutsigCV/reffiles/chr_files_hg19",
    os.path.expanduser(os.path.expandvars("$FULL_GENOME_PATH"))
]
HG19 = None
for t in HG19_LIST:
    if os.path.exists(t):
        HG19 = t
        break

HG19_MMAP_FILES = {}
SCORE_CONF = {'chr': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5}
SCORES = {
    'whole_genome_SNVs.tsv.gz': {
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5, 'element': None
    },
    'hg19_wg_score.tsv.gz': {
        'chr': 0, 'chr_prefix': 'chr', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'hg19_rnasnp_scores.txt.gz': {
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 5, 'element': 6
    },
    'tfbs_creation.tsv.gz': {
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'tfbs_disruption.tsv.gz': {
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    }
}

REGIONS_HEADER = ['chrom', 'start', 'stop', 'feature']
REGIONS_SCHEMA = {
    'fields': {
        'chrom': {'PARSER': 'str(x)', 'VALIDATOR': "x in ([str(c) for c in range(1,23)] + ['X', 'Y'])"},
        'start': {'PARSER': 'int(x)', 'VALIDATOR': 'x > 0'},
        'stop': {'PARSER': 'int(x)', 'VALIDATOR': 'x > 0'}
}}

MUTATIONS_HEADER = ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "TYPE"]
MUTATIONS_SCHEMA = {
    'fields': {
        'CHROMOSOME': {'PARSER': 'str(x)', 'VALIDATOR': "x in ([str(c) for c in range(1,23)] + ['X', 'Y'])"},
        'POSITION':   {'PARSER': 'int(x)', 'VALIDATOR': 'x > 0'},
        'REF':        {'PARSER': 'str(x).upper()', 'VALIDATOR': 'match("^[ACTG-]*$",x)'},
        'ALT':        {'PARSER': 'str(x).upper()', 'VALIDATOR': 'match("^[ACTG-]*$",x)'},
        'TYPE':       {'NULLABLE': 'True', 'VALIDATOR': 'x in ["subs", "indel"]'}
    }
}


def _compress_format(file):
    if file.endswith(".gz"):
        return "gzip"
    if file.endswith(".bz2"):
        return "bz2"
    if file.endswith(".xz"):
        return "xz"
    return None


def _load_regions(regions_file):

    if not os.path.exists(regions_file):
        raise RuntimeError("Feature file '{}' not found".format(regions_file))

    regions = defaultdict(IntervalTree)
    elements = []
    with itab.open(regions_file, header=REGIONS_HEADER, schema=REGIONS_SCHEMA) as reader:
        for r, errors in reader:
            # Report errors and continue
            if len(errors) > 0:
                for e in errors:
                    logging.error(e)
                continue
            elements.append((r[0], r[1], r[2], r[3]))

    for i, r in enumerate(elements):
        if i % 15632 == 0:
            logging.info("[{} of {}]".format(i+1, len(elements)))
        chrom, start, stop, feature = r[0], r[1], r[2], r[3]
        regions[chrom][start:stop] = feature
    logging.info("[{} of {}]".format(i+1, len(elements)))

    return regions


def _load_variants_dict(variants_file, regions_file, signature_name='none'):

    # Load regions
    logging.info("Loading regions")
    regions = _load_regions(regions_file)

    # Load mutations
    variants_dict = defaultdict(list)
    logging.info("Loading and mapping mutations")

    # Check the file format
    for r in load_mutations(variants_file):

        if r.get('TYPE', None) is None:
            if '-' in r['REF'] or '-' in r['ALT'] or len(r['REF']) > 1 or len(r['ALT']) > 1:
                r['TYPE'] = 'indel'
            else:
                r['TYPE'] = 'subs'

        if r['CHROMOSOME'] not in regions:
            continue

        position = int(r['POSITION'])
        for interval in regions[r['CHROMOSOME']][position]:
            variants_dict[interval.data].append({
                'CHROMOSOME': r['CHROMOSOME'],
                'POSITION': position,
                'SAMPLE': r['SAMPLE'],
                'TYPE': r['TYPE'],
                'REF': r['REF'],
                'ALT': r['ALT'],
                'SIGNATURE': signature_name
            })

    return variants_dict


def _create_background_tsv(background_tsv_file, element, regions_file, scores_file):

    if os.path.exists(background_tsv_file):
        logging.debug("%s - background.tsv [skip]", element)
    else:
        if not os.path.exists(regions_file):
            logging.error("%s - background.tsv [error]", element)
            raise RuntimeError("File {} not found".format(regions_file))

        # Load regions at this feature
        all_regions = pd.read_csv(regions_file,
                                  names=['chr', 'start', 'stop', 'feature'],
                                  dtype={'chr': object, 'start': np.int, 'stop': np.int, 'feature': object},
                                  sep='\t')
        regions = all_regions[all_regions.feature == element]

        # Build the command to query with Tabix
        cmd = [TABIX, scores_file]
        cmd += ['{}:{}-{}'.format(region.chr, region.start, region.stop) for idx, region in regions.iterrows()]
        cmd += ['>', background_tsv_file]

        # Execute it
        command = ' '.join(cmd)
        logging.debug("%s", command)
        ret_code = subprocess.call(command, shell=True)

        if ret_code != 0:
            logging.error("%s - background.tsv [error]", element)
            raise RuntimeError("Running {}".format(command))

        _group_writable(background_tsv_file)
        logging.debug("%s - background.tsv [done]", element)


def _get_hg19_mmap(chromosome):
    if chromosome not in HG19_MMAP_FILES:
        fd = open(os.path.join(HG19, "chr{0}.txt".format(chromosome)), 'rb')
        HG19_MMAP_FILES[chromosome] = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    return HG19_MMAP_FILES[chromosome]


def _get_ref_triplet(chromosome, start):
    mm_file = _get_hg19_mmap(chromosome)
    mm_file.seek(start-1)
    return mm_file.read(3).decode().upper()


def _random_scores(num_samples, sampling_size, background, signature):

    values = None

    result = None
    if num_samples > 0:

        if len(background) == 0:
            return None

        if signature is None:
            # Subs sampling without signature
            p_normalized = None

        else:
            # Subs sampling with signature
            p_normalized = signature / sum(signature)

        to_pick = sampling_size - len(values) if values is not None else sampling_size
        if to_pick > 0:
            if num_samples == 1:
                result = np.array(
                    np.random.choice(background, size=to_pick, p=p_normalized, replace=True),
                    dtype='float32'
                )
            else:
                result = np.array(
                    [np.mean(np.random.choice(background, size=num_samples, p=p_normalized, replace=False)) for a in range(to_pick)],
                    dtype='float32'
                )

        # Concat results with previous sampling
        if values is not None:
            if result is not None:
                result = np.concatenate((values, result), axis=0)
            else:
                result = values

    return result[:sampling_size]


def _sampling(sampling_size, scores, signature, item):

    e, m = item[0], item[1]

    values_mean = None
    values_mean_count = 0
    all_scores = []
    for m_tissue in m['muts_by_tissue']['subs'].keys():

        m_scores = m['muts_by_tissue']['subs'][m_tissue]
        m_count = len(m_scores)

        # TODO Use per tissue signature
        values = _random_scores(m_count, sampling_size, scores, signature)

        if values is None:
            logging.warning("There are no scores at {}-{}-{}".format(e, m_count, sampling_size))
            continue

        values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count],
                                 axis=0) if values_mean is not None else values
        values_mean_count += m_count

        all_scores += m_scores

    obs = len(values_mean[values_mean >= np.mean(all_scores)]) if len(all_scores) > 0 else float(sampling_size)

    logging.debug(" %s - sampling.bin [done]", e)

    return e, obs


def _create_background_signature(cache_folder, signature_dict, e):

    cache_path = os.path.join(cache_folder, e[len(e) - 2:], e.upper())
    element_scores_tsv_file = os.path.join(cache_path, 'background.tsv')

    probabilities = []
    scores = []

    with open(element_scores_tsv_file, 'r') as fd:
        reader = csv.reader(fd, delimiter='\t')
        for row in reader:
            # Read score
            value = float(row[SCORE_CONF['score']])
            alt = row[SCORE_CONF['alt']]

            # Expand refseq2 dots
            if alt == '.':
                scores.append(value)
                scores.append(value)
                scores.append(value)
            else:
                scores.append(value)

            # Compute signature
            if signature_dict is not None:
                ref_triplet = _get_ref_triplet(row[SCORE_CONF['chr']], int(row[SCORE_CONF['pos']]) - 1)
                ref = row[SCORE_CONF['ref']]
                alt = row[SCORE_CONF['alt']]

                if ref_triplet[1] != ref:
                    logging.warning("Background mismatch at position %d at '%s'", int(row[SCORE_CONF['pos']]), element_scores_tsv_file)

                # Expand funseq2 dots
                alts = alt if alt != '.' else 'ACGT'.replace(ref, '')

                for a in alts:
                    alt_triplet = ref_triplet[0] + a + ref_triplet[2]
                    try:
                        signature = signature_dict[(ref_triplet, alt_triplet)]
                        probabilities.append(signature)
                    except KeyError:
                        logging.warning("Triplet without probability ref: '%s' alt: '%s' (%s:%s)", ref_triplet, alt_triplet, row[SCORE_CONF['chr']], row[SCORE_CONF['pos']])
                        probabilities.append(0.0)

    signature = np.array(probabilities) if len(probabilities) > 0 else None
    scores = np.array(scores)
    return scores, signature


def _scores_by_position(e, regions_file, scores_file, cache_folder):

    background_folder = os.path.join(cache_folder, e[len(e) - 2:], e.upper())
    background_tsv_file = os.path.join(background_folder, 'background.tsv')

    if not os.path.exists(background_tsv_file):
        _silent_mkdir(background_folder)
        _create_background_tsv(background_tsv_file, e, regions_file, scores_file)

    scores = defaultdict(list)
    with open(background_tsv_file, "r") as fd:
        reader = csv.reader(fd, delimiter='\t')
        for row in reader:
            value = float(row[SCORE_CONF['score']])
            ref = row[SCORE_CONF['ref']]
            alt = row[SCORE_CONF['alt']]
            pos = row[SCORE_CONF['pos']]

            scores[pos].append({'ref': ref, 'alt': alt, 'value': value})

    return scores


def _compute_element(regions_file, scores_file, cache_folder, signature_dict, min_randomizations, max_randomizations, element):

    element, item = _compute_score_means(regions_file, scores_file, cache_folder, element)

    if item is None:
        return element, "The element {} has mutations but not scores.".format(element)

    scores, signature = _create_background_signature(cache_folder, signature_dict, element)

    obs = 0
    randomizations = min_randomizations
    while obs <= 5:
        element, obs = _sampling(randomizations, scores, signature, (element, item))
        if randomizations >= max_randomizations:
            break
        randomizations = min(max_randomizations, randomizations*2)

    item['pvalue'] = max(1, obs) / float(randomizations)
    return element, item


def _compute_score_means(regions_file, scores_file, cache_folder, element):

    element, muts = element[0], element[1]

    # Skip features without mutations
    if len(muts) == 0:
        return element, None

    # Read element scores
    scores = _scores_by_position(element, regions_file, scores_file, cache_folder)

    # Add scores to the element mutations
    muts_by_tissue = {'subs': defaultdict(list)}
    scores_by_sample = {}
    scores_list = []
    scores_subs_list = []
    scores_indels_list = []
    total_subs = 0
    total_subs_score = 0
    positions = []
    mutations = []
    for m in muts:

        if m['TYPE'] == "subs":
            total_subs += 1
            values = scores.get(str(m['POSITION']), [])
            for v in values:
                if v['ref'] == m['REF'] and (v['alt'] == m['ALT'] or v['alt'] == '.'):
                    m['SCORE'] = v['value']
                    total_subs_score += 1
                    break

        # Update scores
        if 'SCORE' in m:

            sample = m['SAMPLE']
            if sample not in scores_by_sample:
                scores_by_sample[sample] = []

            scores_by_sample[sample].append(m['SCORE'])
            scores_list.append(m['SCORE'])

            if m['TYPE'] == "subs":
                scores_subs_list.append(m['SCORE'])
                muts_by_tissue['subs'][m['SIGNATURE']].append(m['SCORE'])

            positions.append(m['POSITION'])
            mutations.append(m)

    if len(scores_list) == 0:
        return element, None

    # Aggregate scores
    num_samples = len(scores_by_sample)

    item = {
        'samples_mut': num_samples,
        'all_mean': np.mean(scores_list),
        'muts': len(scores_list),
        'muts_recurrence': len(set(positions)),
        'samples_max_recurrence': max(Counter(positions).values()),
        'subs': total_subs,
        'subs_score': total_subs_score,
        'scores': scores_list,
        'scores_subs': scores_subs_list,
        'scores_indels': scores_indels_list,
        'positions': positions,
        'muts_by_tissue': muts_by_tissue,
        'mutations': mutations
    }

    return element, item


def _multiple_test_correction(results, num_significant_samples=2):

    results_all = pd.DataFrame.from_dict(results, orient='index')

    # Filter minimum samples
    results_good = results_all[results_all['samples_mut'] >= num_significant_samples].copy()
    results_masked = results_all[results_all['samples_mut'] < num_significant_samples].copy()

    # Multiple test correction
    if len(results_good) > 1:
        results_good['qvalue'] = mlpt(results_good['pvalue'], alpha=0.05, method='fdr_bh')[1]
    else:
        results_good['qvalue'] = np.nan

    # Concat results
    results_concat = pd.concat([results_good, results_masked])
    return results_concat

def load_mutations(file):
    reader = itab.DictReader(file, header=MUTATIONS_HEADER, schema=MUTATIONS_SCHEMA)
    for ix, (row, errors) in enumerate(reader, start=1):
        if len(errors) > 0:
            if reader.line_num == 1:
                # Most probable this is a file with a header
                continue
            for e in errors:
                logging.warning(e)
            continue

        yield row
    reader.fd.close()

def _get_reference_signature(line):
    return _get_ref_triplet(line['CHROMOSOME'], line['POSITION'] - 1)

def _get_alternate_signature(line):
    return line['Signature_reference'][0] + line['ALT'] + line['Signature_reference'][2]

def _compute_signature(variants_file):
    mutations = pd.DataFrame.from_dict([r for r in load_mutations(variants_file) if r['TYPE'] == 'subs'])
    mutations = mutations.groupby(['CHROMOSOME', 'POSITION', 'REF', 'ALT']).count()
    mutations.reset_index(inplace=True)
    mutations['Signature_reference'] = mutations.apply(_get_reference_signature, axis=1)
    mutations['Signature_alternate'] = mutations.apply(_get_alternate_signature, axis=1)
    result = mutations.groupby(['Signature_reference', 'Signature_alternate']).agg({'SAMPLE': 'sum'})
    result.columns = ['count']
    result['probability'] = result['count'] / result['count'].sum()
    return result.to_dict()['probability']

def _file_name(file):
    if file is None:
        return None
    return os.path.splitext(os.path.basename(file))[0]


def _silent_mkdir(folder):
    if not os.path.exists(folder):
        try:
            os.makedirs(folder, mode=0o774)
        except FileExistsError:
            pass


def _group_writable(file):
    st = os.stat(file)
    if not bool(st.st_mode & stat.S_IWGRP):
        try:
            os.chmod(file, 0o664)
        except PermissionError:
            pass