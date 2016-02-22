import gzip
import logging
import os
import numpy as np
import tabix
from scipy import stats

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


def rmean(a):
    n_a = np.array(a)
    return np.mean(n_a) - 0.37687190270063464


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

    # Load signature fold change
    foldchange_vector = None

    current_base = 0
    for region in regions:
        try:
            for row in tb.query("{}{}".format(score_conf['chr_prefix'], region['chrom']), region['start']-1, region['stop']):
                value = read_score(row, score_conf, element)

                ref = row[score_conf['ref']] if 'ref' in score_conf else None
                alt = row[score_conf['alt']] if 'alt' in score_conf else None
                pos = row[score_conf['pos']]

                if score_conf.get('element', None) is not None:
                    if row[score_conf['element']] != element:
                        continue

                scores_by_pos[pos].append({'ref': ref, 'alt': alt, 'value': value})

                # Signature
                # Expand refseq2 dots and fitcons like scores (only one score per postion)
                if alt is None or alt == '.':
                    scores_by_segment[region['segment']].append(value)
                    scores_by_segment[region['segment']].append(value)
                    scores_by_segment[region['segment']].append(value)
                else:
                    scores_by_segment[region['segment']].append(value)

                # Compute signature
                if signature_dict is not None:
                    ref_triplet = get_ref_triplet(row[score_conf['chr']].replace(score_conf['chr_prefix'], ''), int(row[score_conf['pos']]) - 1)
                    ref = row[score_conf['ref']] if 'ref' in score_conf else ref_triplet[1]
                    alt = row[score_conf['alt']] if 'alt' in score_conf else None

                    if ref is not None and ref_triplet[1] != ref:
                        logging.warning("Background mismatch at position %d at '%s'", int(row[score_conf['pos']]), element)

                    # Expand funseq2 dots
                    alts = alt if alt is not None and alt != '.' else 'ACGT'.replace(ref, '')

                    for a in alts:
                        alt_triplet = ref_triplet[0] + a + ref_triplet[2]

                        foldchange = foldchange_vector[int(current_base / 3)] if foldchange_vector is not None else 1
                        current_base += 1

                        try:
                            for k in signature_dict.keys():
                                signature = signature_dict[k].get((ref_triplet, alt_triplet), 0.0)
                                signature_by_segment[region['segment']][k].append(signature*foldchange)
                        except KeyError:
                            missing_signatures[ref_triplet] = alt_triplet
                            for k in signature_dict.keys():
                                signature_by_segment[region['segment']][k].append(0.0)

                else:
                    signature_by_segment = None
        except tabix.TabixError:
            logging.warning("Tabix error at {}='{}{}:{}-{}'".format(element, score_conf['chr_prefix'], region['chrom'], region['start']-1, region['stop']))
            continue

    return scores_by_pos, scores_by_segment, signature_by_segment, [(k, v) for k, v in missing_signatures.items()]


def random_scores(num_samples, sampling_size, background, signature, statistic_name):

    result = None
    if num_samples > 0:

        if len(background) == 0:
            return None

        if signature is None:
            # Subs sampling without signature
            p_normalized = None

        else:
            # Subs sampling with signature
            p_normalized = np.array(signature[:len(background)]) / sum(signature)

        if num_samples == 1:
            result = np.array(
                np.random.choice(background, size=sampling_size, p=p_normalized, replace=True),
                dtype='float32'
            )
        else:

            # Select mean
            if statistic_name == 'gmean':
                statistic_test = gmean
            elif statistic_name == 'rmean':
                statistic_test = rmean
            elif statistic_name == 'max':
                statistic_test = np.max
            else:
                statistic_test = np.mean

            result = np.array(
                [statistic_test(np.random.choice(background, size=num_samples, p=p_normalized, replace=True)) for a in range(sampling_size)],
                dtype='float32'
            )

    return result[:sampling_size]


def sampling(sampling_size, scores_by_segment, signature_by_segment, e, m, statistic_name):

    values_mean = None
    values_mean_count = 0

    # Substitutions sampling
    for m_tissue, muts_by_segment in m['muts_by_tissue']['subs'].items():

        # Do randomization per segment
        for segment, m_scores in muts_by_segment.items():
            m_count = len(m_scores)

            signature = signature_by_segment[segment][m_tissue] if signature_by_segment is not None else None
            scores = scores_by_segment[segment]
            values = random_scores(m_count, sampling_size, scores, signature, statistic_name)

            if values is None:
                logging.warning("There are no scores at {}-{}-{}".format(e, m_count, sampling_size))
                continue

            if values_mean is None:
                values_mean = values
            elif statistic_name == 'gmean':
                values_mean = gmean_weighted([values_mean, values], [values_mean_count, m_count])
            elif statistic_name == 'rmean':
                values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count], axis=0)
            elif statistic_name == 'max':
                values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count], axis=0)
            else:
                values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count], axis=0)

            values_mean_count += m_count

    # Select mean
    if statistic_name == 'gmean':
        statistic_test = gmean
    elif statistic_name == 'rmean':
        statistic_test = np.mean
    elif statistic_name == 'max':
        statistic_test = np.mean
    else:
        statistic_test = np.mean

    if statistic_name == 'max':
        max_by_sample = [max(sum(v.values(), [])) for v in m['muts_by_tissue']['subs'].values()]
        weight_by_sample = [len(sum(v.values(), [])) for v in m['muts_by_tissue']['subs'].values()]
        obs_mean = np.average(max_by_sample, weights=weight_by_sample, axis=0)
    else:
        obs_mean = statistic_test(m['scores'])

    obs = len(values_mean[values_mean >= obs_mean]) if len(m['scores']) > 0 else float(sampling_size)
    neg_obs = len(values_mean[values_mean <= obs_mean]) if len(m['scores']) > 0 else float(sampling_size)

    return obs, neg_obs


def compute_muts_statistics(muts, scores, statistic_name, recurrence):

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

            if recurrence or m['POSITION'] not in positions:

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

    # Aggregate scores
    num_samples = len(scores_by_sample)

    item = {
        'samples_mut': num_samples,
        'muts': len(scores_list),
        'muts_recurrence': len(set(positions)),
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