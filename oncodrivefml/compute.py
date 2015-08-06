import gzip
import logging
import os
import pickle
import pandas as pd
import numpy as np
import tabix
from scipy import stats

from statsmodels.sandbox.stats.multicomp import multipletests as mlpt
from collections import Counter, defaultdict
from oncodrivefml.signature import get_ref_triplet


def gmean(a):
    return stats.gmean(np.array(a) + 1.0) - 1.0


def read_score(row, score_conf, element):
    value_str = row[score_conf['score']]
    if value_str is None or value_str == '':
        if 'extra' in score_conf:
            value_str = row[score_conf['extra']].split(',')
            value = 0
            for val in value_str:
                elm, vals = val.split(':')
                if elm == element:
                    value = float(vals)
                    break
        else:
            value = 0
    else:
        value = float(value_str)

    return value


def load_scores(element, regions, signature_dict, score_conf):
    tb = tabix.open(score_conf['file'])
    scores_by_pos = defaultdict(list)
    scores_by_segment = defaultdict(list)
    signature_by_segment = defaultdict(lambda: defaultdict(list))
    missing_signatures = {}

    for region in regions:
        for row in tb.query("{}{}".format(score_conf['chr_prefix'], region['chrom']), region['start'], region['stop']):
            value = read_score(row, score_conf, element)
            ref = row[score_conf['ref']]
            alt = row[score_conf['alt']]
            pos = row[score_conf['pos']]

            scores_by_pos[pos].append({'ref': ref, 'alt': alt, 'value': value})

            # Signature
            # Expand refseq2 dots
            if alt == '.':
                scores_by_segment[region['segment']].append(value)
                scores_by_segment[region['segment']].append(value)
                scores_by_segment[region['segment']].append(value)
            else:
                scores_by_segment[region['segment']].append(value)

            # Compute signature
            if signature_dict is not None:
                ref_triplet = get_ref_triplet(row[score_conf['chr']].replace(score_conf['chr_prefix'], ''), int(row[score_conf['pos']]) - 1)
                ref = row[score_conf['ref']]
                alt = row[score_conf['alt']]

                if ref_triplet[1] != ref:
                    logging.warning("Background mismatch at position %d at '%s'", int(row[score_conf['pos']]), element)

                # Expand funseq2 dots
                alts = alt if alt != '.' else 'ACGT'.replace(ref, '')

                for a in alts:
                    alt_triplet = ref_triplet[0] + a + ref_triplet[2]
                    try:
                        for k in signature_dict.keys():
                            signature = signature_dict[k][(ref_triplet, alt_triplet)]
                            signature_by_segment[region['segment']][k].append(signature)
                    except KeyError:
                        missing_signatures[ref_triplet] = alt_triplet
                        for k in signature_dict.keys():
                            signature_by_segment[region['segment']][k].append(0.0)

    return scores_by_pos, scores_by_segment, signature_by_segment, [(k, v) for k, v in missing_signatures.items()]


def random_scores(num_samples, sampling_size, background, signature, geometric):

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

                # Select mean
                mean = gmean if geometric else np.mean

                result = np.array(
                    [mean(np.random.choice(background, size=num_samples, p=p_normalized, replace=False)) for a in range(to_pick)],
                    dtype='float32'
                )

        # Concat results with previous sampling
        if values is not None:
            if result is not None:
                result = np.concatenate((values, result), axis=0)
            else:
                result = values

    return result[:sampling_size]


def sampling(sampling_size, scores_by_segment, signature_by_segment, e, m, geometric):

    values_mean = None
    values_mean_count = 0
    all_scores = []
    for m_tissue, muts_by_segment in m['muts_by_tissue']['subs'].items():

        # Do randomization per segment
        for segment, m_scores in muts_by_segment.items():
            m_count = len(m_scores)

            signature = signature_by_segment[segment][m_tissue]
            scores = scores_by_segment[segment]
            values = random_scores(m_count, sampling_size, scores, signature, geometric)

            if values is None:
                logging.warning("There are no scores at {}-{}-{}".format(e, m_count, sampling_size))
                continue

            values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count],
                                     axis=0) if values_mean is not None else values
            values_mean_count += m_count
            all_scores += m_scores

    # Select mean
    mean = gmean if geometric else np.mean

    obs = len(values_mean[values_mean >= mean(all_scores)]) if len(all_scores) > 0 else float(sampling_size)

    return e, obs


def compute_element(signature_dict, min_randomizations, max_randomizations, geometric, score_conf, input_data):

    element, muts, regions, trace = input_data

    # Load all element scores
    scores_by_position, scores_by_segment, signature_by_segment, missing_signatures = load_scores(
        element,
        regions,
        signature_dict,
        score_conf
    )

    # Compute elements statistics
    item = compute_muts_statistics(muts, scores_by_position, geometric)

    if item is None:
        return element, "The element {} has mutations but not scores.".format(element), []

    obs = 0
    randomizations = min_randomizations
    while obs <= 5:
        element, obs = sampling(randomizations, scores_by_segment, signature_by_segment, element, item, geometric)
        if randomizations >= max_randomizations or obs > 5:
            break
        randomizations = min(max_randomizations, randomizations*2)

    item['pvalue'] = max(1, obs) / float(randomizations)

    if trace is not None:
        with gzip.open(trace, 'wb') as fd:

            # Remove defaultdict lambdas
            item['muts_by_tissue'] = {k: dict(v) for k, v in item['muts_by_tissue'].items()}
            signature_by_segment = {k: dict(v) for k, v in signature_by_segment.items()}

            trace_object = {
                'item': item,
                'scores_by_segment': scores_by_segment,
                'signature_by_segment': signature_by_segment
            }

            pickle.dump(trace_object, fd)

    # Remove this key to simplify serialization
    del item['muts_by_tissue']

    return element, item, missing_signatures


def compute_muts_statistics(muts, scores, geometric):

    # Skip features without mutations
    if len(muts) == 0:
        return None

    # Add scores to the element mutations
    muts_by_tissue = {'subs': defaultdict(lambda: defaultdict(list))}
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
                muts_by_tissue['subs'][m['SIGNATURE']][m['SEGMENT']].append(m['SCORE'])

            positions.append(m['POSITION'])
            mutations.append(m)

    if len(scores_list) == 0:
        return None

    # Aggregate scores
    num_samples = len(scores_by_sample)

    # Select mean
    mean = gmean if geometric else np.mean

    item = {
        'samples_mut': num_samples,
        'all_mean': mean(scores_list),
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

    return item


def multiple_test_correction(results, num_significant_samples=2):

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


def file_name(file):
    if file is None:
        return None
    return os.path.splitext(os.path.basename(file))[0]


def silent_mkdir(folder):
    if not os.path.exists(folder):
        try:
            os.makedirs(folder, mode=0o774)
        except FileExistsError:
            pass