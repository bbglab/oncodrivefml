import argparse
import json
import logging
import signaturesampling as sampling
import os
import numpy as np
import csv
import gzip

from multiprocessing import Pool
from collections import Counter
from oncodrivefm2_step2 import _signature

STEP0_VERSION = "v01"
STEP1_VERSION = "v03"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Constants
NC_PATH = "/projects_bg/bg/shared/datasets/"
ENSEMBLE_FILE = os.path.join(NC_PATH, "non-coding/ensemble_genes_75.txt")
CHROMOSOMES = [str(i) for i in range(1, 23)]
CHROMOSOMES.extend(['X', 'Y'])
GENE_CONVERSION = {line.split("\t")[0]: line.strip().split("\t")[-1]
                   for line in open(ENSEMBLE_FILE, 'r').readlines()
                   if line.split("\t")[1] in CHROMOSOMES}

def _load_indels(project, score, tissue):
    score_conf = sampling.SCORES[score]
    indels_file = os.path.join(sampling.INPUT_FOLDER, "indels", score, project, "{}.txt.gz".format(tissue))
    indels = {}
    if os.path.exists(indels_file):
        with gzip.open(indels_file, 'rt') as fd:
            reader = csv.reader(fd, delimiter='\t')

            # Skip header
            next(reader)

            for row in reader:
                value = float(row[score_conf['score']])
                ref = row[score_conf['ref']]
                alt = row[score_conf['alt']]
                chrom = row[score_conf['chr']]
                pos = row[score_conf['pos']]
                key = "{}-{}".format(chrom, pos)

                if key not in indels:
                    indels[key] = []

                item = {'ref': ref, 'alt': alt, 'value': value}
                if score_conf['element'] is not None:
                    item['element'] = row[score_conf['element']]
                indels[key].append(item)
    return indels


def _scores(e, feature, score):
    score_conf = sampling.SCORES[score]

    background_folder = os.path.join(sampling.SAMPLING_CACHE_FOLDER, score, feature, e[len(e) - 2:], e.upper())
    background_tsv_file = os.path.join(background_folder, 'background.tsv')

    if not os.path.exists(background_tsv_file):
        sampling._silent_mkdir(background_folder)
        sampling.create_background_tsv(score, background_tsv_file, "none", feature, e)

    scores = {}
    with open(background_tsv_file, "r") as fd:
        reader = csv.reader(fd, delimiter='\t')
        for row in reader:
            value = float(row[score_conf['score']])
            ref = row[score_conf['ref']]
            alt = row[score_conf['alt']]
            pos = row[score_conf['pos']]

            # Skip scores of other elements
            if score_conf['element'] is not None:
                back_element = row[score_conf['element']]
                if back_element != e:
                    continue

            if pos not in scores:
                scores[pos] = []

            scores[pos].append({'ref': ref, 'alt': alt, 'value': value})

    return scores


def _add_by_tissue(scores_by_tissue, tissue, score):

    if tissue not in scores_by_tissue:
        scores_by_tissue[tissue] = []

    scores_by_tissue[tissue].append(score)


def _compute_score_means(elements, feature, indels, score, chunk):

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
        scores = _scores(element, feature, score)

        # Add scores to the element mutations
        muts_by_tissue = {'subs': {}, 'indel': {}}
        scores_by_sample = {}
        scores_list = []
        scores_subs_list = []
        scores_indels_list = []
        total_subs = 0
        total_subs_score = 0
        total_indels = 0
        total_indels_score = 0
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

            elif m['TYPE'] == "indel":
                total_indels += 1
                #if score.startswith('cadd') and len(indels) > 0:
                #    key = '{}-{}'.format(m['CHROMOSOME'], m['POSITION'])
                #    values = indels.get(key, [])
                #    for v in values:
                #        if v['ref'] == m['REF'] and (v['alt'] == m['ALT'] or v['alt'] == '.'):
                #            m['SCORE'] = v['value']
                #            total_indels_score += 1
                #            break
                #elif score.startswith('funseq') or score.startswith('cadd'):
                pos = m['POSITION']
                seq = m['REF'] if m['ALT'] == '-' else m['ALT']
                seq = seq.upper()

                # Skip indels longer than 50
                if len(seq) < 50:
                    indels_scores = []
                    for pos, n in enumerate(seq, start=pos):
                        if n not in "ACTG":
                            continue

                        values = scores.get(str(pos), [])

                        for v in values:
                            indels_scores.append(v['value'])

                    if len(indels_scores) > 0:
                        m['SCORE'] = max(indels_scores)
                        m['INDEL_SIZE'] = len(indels_scores)
                        total_indels_score += 1

            # Update scores
            if 'SCORE' in m:

                sample = m['SAMPLE']
                if sample not in scores_by_sample:
                    scores_by_sample[sample] = []

                scores_by_sample[sample].append(m['SCORE'])
                scores_list.append(m['SCORE'])

                if m['TYPE'] == "subs":
                    scores_subs_list.append(m['SCORE'])
                    _add_by_tissue(muts_by_tissue['subs'], m['PROJECT'], m['SCORE'])
                elif m['TYPE'] == "indel":
                    scores_indels_list.append(m['SCORE'])
                    _add_by_tissue(muts_by_tissue['indel'], m.get('INDEL_SIZE', 0), m['SCORE'])

                positions.append(m['POSITION'])
                mutations.append(m)

        if len(scores_list) == 0:
            continue

        # Aggregate scores
        scores_by_sample_maxs = []
        for s in scores_by_sample:
            scores_by_sample_maxs.append(max(scores_by_sample[s]))
        num_samples = len(scores_by_sample)

        item = {
            'samples_mut': num_samples,
            'all_mean': np.mean(scores_list),
            'max_mean': np.mean(scores_by_sample_maxs),
            'muts': len(scores_list),
            'muts_recurrence': len(set(positions)),
            'samples_max_recurrence': max(Counter(positions).values()),
            'subs': total_subs,
            'subs_score': total_subs_score,
            'indels': total_indels,
            'indels_score': total_indels_score,
            'symbol': GENE_CONVERSION.get(element, "UNKNOWN"),
            'scores': scores_list,
            'scores_subs': scores_subs_list,
            'scores_indels': scores_indels_list,
            'positions': positions,
            'muts_by_tissue': muts_by_tissue,
            'mutations': mutations
        }

        result.append((element, item))

    return result


def _prepare_background(feature, pool, prepare_arguments, elements, score):
    to_background = []
    for e in elements:
        background_folder = os.path.join(sampling.SAMPLING_CACHE_FOLDER, score, feature, e[len(e) - 2:], e.upper())
        background_tsv_file = os.path.join(background_folder, 'background.tsv')
        if not os.path.exists(background_tsv_file):
            if not os.path.exists(os.path.join(background_folder, '.lock')):
                to_background.append((e, 0))
    if len(to_background) > 0:
        logger.info("Send {} backgrounds".format(len(to_background)))
        return pool.apply(sampling.sampling_prepare, (to_background,), prepare_arguments)
    else:
        return None


def _run(project, tissue, score, feature):

    # Skip if done
    output_folder = os.path.join(sampling.OUTPUT_FOLDER, project, tissue, feature)
    results_file = os.path.join(output_folder, "results_{}_{}.json.gz".format(score, STEP1_VERSION))
    if os.path.exists(results_file):
        logger.info("{} - {} - {} - {} [skip]".format(project, tissue, score, feature))
        return

    # Loading mutations
    logger.info("Loading mutations")
    output_folder = os.path.join(sampling.OUTPUT_FOLDER, project, tissue, feature)
    mutations_file = os.path.join(output_folder, "mutations_{}.json.gz".format(STEP0_VERSION))
    if not os.path.exists(mutations_file):
        raise RuntimeError("Mutations file '{}' not found".format(mutations_file))
    with gzip.open(mutations_file, 'rt') as fd:
        elements = json.load(fd)

    # Initialize
    cores = 12
    pool = Pool(cores)
    signature = _signature(project, tissue) if tissue != 'pan' else 'none'
    prepare_arguments = dict(feature=feature, signature=signature, score=score, verbose=True, max_jobs=50)
    if score not in sampling.SCORES:
        raise RuntimeError("Unknown score '{}'".format(score))

    # Load indels scores
    #logger.info("Loading indels scores")
    #indels = _load_indels(project, score, tissue)
    indels = {}

    # Prepare background
    _prepare_background(feature, pool, prepare_arguments, elements, score)

    # Compute elements statistics
    logger.info("Computing score statistics by element")
    results = {}
    elements_all = list(elements.items())
    elements_chunks = np.array_split(elements_all, cores)
    compute_arguments = ((chunk, feature, indels, score, num) for num, chunk in enumerate(elements_chunks, start=1))
    c = 0
    for done in pool.starmap(_compute_score_means, compute_arguments):
        c += 1
        logger.info("Chunk {} of {} [done]".format(c, cores))
        for element, means in done:
            results[element] = means

    # Store results
    logger.info("Store results")
    sampling._silent_mkdir(output_folder)
    with gzip.open(results_file, 'wt') as fd:
        json.dump(results, fd)

    logger.info("{} - {} - {} - {} [done]".format(project, tissue, score, feature))


if __name__ == "__main__":

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')

    # Parse the arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("project")
    parser.add_argument("tissue")
    parser.add_argument("score")
    parser.add_argument("feature")

    args = parser.parse_args()

    _run(args.project, args.tissue, args.score, args.feature)
