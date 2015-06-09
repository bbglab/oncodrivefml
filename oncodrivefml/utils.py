import csv
import logging
import os
import stat
import mmap
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt
from intervaltree import IntervalTree
import subprocess

TABIX = "/soft/bio/sequence/tabix-0.2.3/tabix"
SCORE_CONF = {'chr': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5}
HG19_DIR = "/projects_bg/bg/soft/intogen_home/gencluster/software/mutsigCV/reffiles/chr_files_hg19"
HG19_MMAP_FILES = {}

def _compress_format(file):
    if file.endswith("gz"):
        return "gzip"
    if file.endswith("bz"):
        return "bz2"
    return None


def _load_regions(regions_file):
    if not os.path.exists(regions_file):
        raise RuntimeError("Feature file '{}' not found".format(regions_file))

    regions = {}
    with open(regions_file, 'rt') as fd:
        reader = csv.reader(fd, delimiter='\t')
        for r in reader:
            chrom, start, stop, feature = r[0], r[1], r[2], r[3]
            if chrom not in regions:
                regions[chrom] = IntervalTree()
            regions[chrom][int(start):int(stop)] = feature
    return regions


def _load_variants_dict(variants_file, regions_file, signature_name='none'):

    # Load regions
    logging.info("Loading regions")
    regions = _load_regions(regions_file)

    # Load mutations
    variants_dict = defaultdict(list)
    logging.info("Loading mutations")
    with open(variants_file) as fd:
        reader = csv.DictReader(fd, delimiter='\t')
        for r in reader:
            if r['TYPE'] not in ['subs', 'indel']:
                continue

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
        fd = open(os.path.join(HG19_DIR, "chr{0}.txt".format(chromosome)), 'rb')
        HG19_MMAP_FILES[chromosome] = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    return HG19_MMAP_FILES[chromosome]


def _get_ref_triplet(chromosome, start):
    mm_file = _get_hg19_mmap(chromosome)
    mm_file.seek(start-1)
    return mm_file.read(3).decode().upper()


def _random_scores(num_samples, sampling_size, scores_file, signature_file, sampling_file):

    values = None
    if os.path.exists(sampling_file):
        values = np.fromfile(sampling_file, dtype='float32').astype('float64')

    result = None
    if num_samples > 0:
        background = np.fromfile(scores_file, dtype='float32').astype('float64')

        if len(background) == 0:
            return None

        if signature_file is None:

            # Subs sampling without signature
            p_normalized = None

        else:
            # Subs sampling with signature
            probabilities = np.fromfile(signature_file, dtype='float32').astype('float64')
            p_normalized = probabilities / sum(probabilities)

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

        # Store the sampling
        result.tofile(sampling_file, format='float32')
        _group_writable(sampling_file)

    return result[:sampling_size]


def _sampling(chunks, sampling_size, cache_folder, process):

    result = []
    total = len(chunks)
    for i, (e, m) in enumerate(chunks):

        if i % 50 == 0:
            print("[process-{}] {} of {}".format(process, i, total))
        i += 1

        cache_path = os.path.join(cache_folder, e[len(e) - 2:], e.upper())
        scores_file = os.path.join(cache_path, 'background.bin')
        signature_file = os.path.join(cache_path, 'signature.bin')
        if not os.path.exists(signature_file):
            signature_file = None
        sampling_file = os.path.join(cache_path, 'sampling.bin')

        values_mean = None
        values_mean_count = 0
        all_scores = []
        for m_tissue in m['muts_by_tissue']['subs'].keys():

            m_scores = m['muts_by_tissue']['subs'][m_tissue]
            m_count = len(m_scores)

            # TODO Use per tissue signature
            values = _random_scores(m_count, sampling_size, scores_file, signature_file, sampling_file)

            if values is None:
                logging.warning("There are no scores at {}-{}-{}".format(e, m_count, sampling_size))
                continue

            values_mean = np.average([values_mean, values], weights=[values_mean_count, m_count],
                                     axis=0) if values_mean is not None else values
            values_mean_count += m_count

            all_scores += m_scores

        m['obs'] = len(values_mean[values_mean >= np.mean(all_scores)]) if len(all_scores) > 0 else float(sampling_size)
        m['pvalue'] = max(1, m['obs']) / float(sampling_size)
        result.append((e, m))

        logging.debug(" %s - sampling.bin [done]", e)

    return result


def _create_background_signature(elements, cache_folder, signature_dict, chunk):

    i = 0
    total = len(elements)
    for e in elements:

        if i % 50 == 0:
            print("[process-{}] {} of {}".format(chunk, i, total))
        i += 1

        cache_path = os.path.join(cache_folder, e[len(e) - 2:], e.upper())
        element_scores_tsv_file = os.path.join(cache_path, 'background.tsv')
        element_scores_bin_file = os.path.join(cache_path, 'background.bin')
        element_signature_bin_file = os.path.join(cache_path, 'signature.bin')

        if os.path.exists(element_scores_bin_file) and os.path.exists(element_signature_bin_file):
            logging.debug("%s - background.bin signature.bin [skip]", e)
        else:

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
                            logging.warning("Triplet without probability ref: '%s' alt: '%s'", ref_triplet, alt_triplet)
                            probabilities.append(0.0)

            # Save background.bin
            np.array(scores, dtype='float32').tofile(element_scores_bin_file, format='float32')
            _group_writable(element_scores_bin_file)

            # Save signature.bin
            c_array = np.array([p for p in probabilities], dtype='float32')
            c_array.tofile(element_signature_bin_file, format='float32')
            _group_writable(element_signature_bin_file)

            logging.debug("%s - signature.bin [done]", e)

    return True


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


def _compute_score_means(elements, regions_file, scores_file, cache_folder, process):

    result = []

    i = 1
    total = len(elements)

    for element, muts in elements:

        if i % 50 == 0:
            print("[process-{}] {} of {}".format(process, i, total))
        i += 1

        # Skip features without mutations
        if len(muts) == 0:
            continue

        ## SCORES PHASE

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