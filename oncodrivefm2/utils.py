import csv
import logging
import os
import stat
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt
import subprocess

TABIX = "/home/jdeu/programs/tabix-0.2.6/tabix"
SCORE_CONF = {'chr': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5}


def _load_variants(variants_file, signature=None):

    if not os.path.exists(variants_file):
        raise RuntimeError("Input dataset not found '{}'".format(variants_file))

    muts = pd.read_csv(variants_file, sep='\t',
                       dtype={'CHROMOSOME': object, 'POSITION': np.int, 'REF': object, 'ALT': object, 'SAMPLE': object, 'TYPE': object})

    if signature is not None:
        muts['SIGNATURE'] = signature

    muts = muts[muts.TYPE.isin(['subs', 'indel'])]
    muts.set_index(['CHROMOSOME', 'POSITION'], inplace=True)
    muts.sortlevel(level=0, axis=0, inplace=True)

    return muts


def _load_regions(regions_file):
    if not os.path.exists(regions_file):
        raise RuntimeError("Regions file '{}' not found".format(regions_file))
    regions = pd.read_csv(regions_file, names=['chr', 'start', 'stop', 'feature'],
                          dtype={'chr': object, 'start': np.int, 'stop': np.int, 'feature': object}, sep='\t')
    regions = regions.groupby(['feature'])
    return regions


def _intersect_mutations(elements, muts, chunk):

    result = []

    i = 1
    total = len(elements)
    valid_chromosomes = set(muts.index.levels[0])

    for element, segments in elements:

        if i % 100 == 0:
            print("[process-{}] {} of {}".format(chunk, i, total))
        i += 1

        # Lookup the mutations at this feature
        element_muts = []
        for r in segments.iterrows():
            chrom = r[1].chr
            if chrom not in valid_chromosomes:
                continue
            start = r[1].start
            stop = r[1].stop
            muts_region = muts.loc(axis=0)[chrom, start:stop]
            if len(muts_region) > 0:
                element_muts += muts_region.reset_index().to_dict(orient='record')

        # Skip features without mutations
        if len(element_muts) == 0:
            continue

        result.append((element, element_muts))

    return result


def _create_background_tsv(background_tsv_file, element, regions_file, scores_file):
    if os.path.exists(background_tsv_file):
        logging.info("%s - background.tsv [skip]", element)
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
        logging.info("%s - background.tsv [done]", element)


def _scores_by_position(e, regions_file, scores_file, output_folder):

    background_folder = os.path.join(output_folder, "cache", e[len(e) - 2:], e.upper())
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


def _compute_score_means(elements, regions_file, scores_file, output_folder, chunk):

    result = []

    i = 1
    total = len(elements)

    for element, element_muts in elements:

        if i % 100 == 0:
            print("[process-{}] {} of {}".format(chunk, i, total))
        i += 1

        # Skip features without mutations
        if len(element_muts) == 0:
            continue

        # Read element scores
        scores = _scores_by_position(element, regions_file, scores_file, output_folder)

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
        for m in element_muts:

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
            continue

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

        result.append((element, item))

    return result


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


def _file_name(file):
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