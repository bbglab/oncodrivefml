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


def gmean_weighted(vectors, weights):
    v_a = [np.array(i) + 1 for i in vectors]
    v_b = [np.power(i, wi) for i, wi in zip(v_a,weights)]
    total_weight = np.sum(weights)
    v_c = [(np.prod([j[i] for j in v_b])**(1/total_weight)) - 1.0 for i in range(len(vectors[0]))]
    return np.array(v_c)


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


def load_scores(element, regions, signature_dict, score_conf, signature_ratio):
    tb = tabix.open(score_conf['file'])
    scores_by_pos = defaultdict(list)
    scores_by_segment = defaultdict(list)
    signature_by_segment = defaultdict(lambda: defaultdict(list))
    missing_signatures = {}

    # Load signature fold change
    foldchange_vector = None
    if signature_ratio is not None:
        signature_foldchange_file = os.path.join(signature_ratio, element[-2:], "{}.bin".format(element))
        if os.path.exists(signature_foldchange_file):
            foldchange_vector = np.fromfile(signature_foldchange_file, dtype='float32')
        else:
            logging.warning("Missing signature fold change vector at '{}'".format(element))

    current_base = 0
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

                    foldchange = foldchange_vector[int(current_base / 3)] if foldchange_vector is not None else 1
                    current_base += 1

                    try:
                        for k in signature_dict.keys():
                            signature = signature_dict[k][(ref_triplet, alt_triplet)]
                            signature_by_segment[region['segment']][k].append(signature*foldchange)
                    except KeyError:
                        missing_signatures[ref_triplet] = alt_triplet
                        for k in signature_dict.keys():
                            signature_by_segment[region['segment']][k].append(0.0)


            else:
                signature_by_segment = None

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


def sampling(sampling_size, scores_by_segment, signature_by_segment, e, m, geometric, indels_background, trace=False):

    values_mean = None
    values_mean_count = 0

    trace_dict = {}

    # Substitutions sampling
    for m_tissue, muts_by_segment in m['muts_by_tissue']['subs'].items():

        # Do randomization per segment
        for segment, m_scores in muts_by_segment.items():
            m_count = len(m_scores)

            signature = signature_by_segment[segment][m_tissue] if signature_by_segment is not None else None
            scores = scores_by_segment[segment]
            values = random_scores(m_count, sampling_size, scores, signature, geometric)

            if trace:
                trace_dict["subs-{}-{}-{}".format(m_tissue, segment, m_count)] = [s for s in values]

            if values is None:
                logging.warning("There are no scores at {}-{}-{}".format(e, m_count, sampling_size))
                continue

            if values_mean is None:
                values_mean = values
            elif geometric:
                values_mean = gmean_weighted([values_mean, values], [values_mean_count, m_count])
            else:
                values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count], axis=0)

            values_mean_count += m_count

    # Select mean
    mean = gmean if geometric else np.mean

    # Indels sampling
    if indels_background is not None:
        for m_tissue, muts_by_segment in m['muts_by_tissue']['indel'].items():

            # Count total indels
            m_count = 0
            for segment, m_scores in muts_by_segment.items():
                m_count += len(m_scores)

            # Load random scores
            if m_count > 0:
                indels_file = os.path.join(indels_background, e[-2:], e, "{}_{}.bin".format(e, m_tissue))
                if os.path.exists(indels_file):
                    values = np.fromfile(indels_file, dtype='float32')

                    if m_count > 1:
                        values = np.array([mean(values[i:i+m_count]) for i in range(0, len(values), m_count)])

                    remaining_values = 0 if values_mean is None else (len(values_mean) - len(values))
                    if remaining_values > 0:
                        # Select uniformly randomly the needed values
                        values_random = np.random.choice(values, size=remaining_values, replace=True)
                        values = np.concatenate((values, values_random))

                    # Add random scores to the mean
                    if values_mean is None:
                        values_mean = values
                    elif geometric:
                        values_mean = gmean_weighted([values_mean, values], [values_mean_count, m_count])
                    else:
                        values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count], axis=0)

                    values_mean_count += m_count

                    if trace:
                        trace_dict["indel-{}-{}-{}".format(m_tissue, e, m_count)] = [s for s in values]

                else:
                    raise RuntimeError("Indels background file '{}' not found.".format(indels_file))

    if values_mean is None:
        return e, None, trace_dict

    obs = len(values_mean[values_mean >= mean(m['scores'])]) if len(m['scores']) > 0 else float(sampling_size)

    return e, obs, trace_dict


def compute_element(signature_dict, min_randomizations, max_randomizations, geometric, score_conf, indels_background, signature_ratio, input_data):

    element, muts, regions, trace_file = input_data

    # Check if we need to trace this element
    trace = trace_file is not None

    # Load all element scores
    scores_by_position, scores_by_segment, signature_by_segment, missing_signatures = load_scores(
        element,
        regions,
        signature_dict,
        score_conf,
        signature_ratio
    )

    # Compute elements statistics
    item = compute_muts_statistics(muts, scores_by_position, geometric)

    if item is None:
        return element, "The element {} has mutations but not scores.".format(element), []

    obs = 0
    randomizations = min_randomizations
    sampling_values = {}
    while obs <= 5:
        element, obs, sampling_values = sampling(randomizations, scores_by_segment, signature_by_segment, element, item, geometric, indels_background, trace=trace)
        if randomizations >= max_randomizations or obs is None or obs > 5:
            break
        randomizations = min(max_randomizations, randomizations*2)

    item['pvalue'] = max(1, obs) / float(randomizations) if obs is not None else None

    if trace:
        with gzip.open(trace_file, 'wb') as fd:

            # Remove defaultdict lambdas
            item['muts_by_tissue'] = None if item['muts_by_tissue'] is None else {k: dict(v) for k, v in item['muts_by_tissue'].items()}
            signature_by_segment = None if signature_by_segment is None else {k: dict(v) for k, v in signature_by_segment.items()}

            trace_object = {
                'item': item,
                'scores_by_segment': scores_by_segment,
                'signature_by_segment': signature_by_segment,
                'sampling': sampling_values
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
    muts_by_tissue = {'subs': defaultdict(lambda: defaultdict(list)), 'indel': defaultdict(lambda: defaultdict(list))}
    scores_by_sample = {}
    scores_list = []
    scores_subs_list = []
    scores_indels_list = []
    total_subs = 0
    total_subs_score = 0
    positions = []
    mutations = []
    for m in muts:

        # Get substitutions scores
        if m['TYPE'] == "subs":
            total_subs += 1
            values = scores.get(str(m['POSITION']), [])
            for v in values:
                if v['ref'] == m['REF'] and (v['alt'] == m['ALT'] or v['alt'] == '.'):
                    m['SCORE'] = v['value']
                    total_subs_score += 1
                    break

        # Update scores
        if m.get('SCORE', None) is not None:

            sample = m['SAMPLE']
            if sample not in scores_by_sample:
                scores_by_sample[sample] = []

            scores_by_sample[sample].append(m['SCORE'])
            scores_list.append(m['SCORE'])

            if m['TYPE'] == "subs":
                scores_subs_list.append(m['SCORE'])
                muts_by_tissue['subs'][m['SIGNATURE']][m['SEGMENT']].append(m['SCORE'])
            elif m['TYPE'] == "indel":
                scores_indels_list.append(m['SCORE'])
                muts_by_tissue['indel'][m['SIGNATURE']][m['SEGMENT']].append(m['SCORE'])

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
    results_good = results_all[(results_all['samples_mut'] >= num_significant_samples) & (~results_all['pvalue'].isnull())].copy()
    results_masked = results_all[(results_all['samples_mut'] < num_significant_samples) | (results_all['pvalue'].isnull())].copy()

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
    return os.path.basename(file).split('.')[0]


def silent_mkdir(folder):
    if not os.path.exists(folder):
        try:
            os.makedirs(folder, mode=0o774)
        except FileExistsError:
            pass