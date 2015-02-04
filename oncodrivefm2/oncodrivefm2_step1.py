import json
import logging
import math
from configuration import Configuration, SCORES
import os
import numpy as np
import csv
import gzip

from multiprocessing import Pool
from collections import Counter
from signaturesampling import SignatureSampling, silent_mkdir

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Constants
CHROMOSOMES = [str(i) for i in range(1, 23)]
CHROMOSOMES.extend(['X', 'Y'])


class OncodriveFM2Error(RuntimeError):
    pass


class OncodriveFM2Scoring(object):
    def __init__(self, config: Configuration, sampling: SignatureSampling):
        self._discarded = []
        self._conf = config
        self._sampling = sampling
        self.gene_conversion = {line.split("\t")[0]: line.strip().split("\t")[-1]
                   for line in open(self._conf.ensembl_file, 'r').readlines()
                   if line.split("\t")[1] in CHROMOSOMES}
        self._output_file = None


    def get_output_file(self):
        return self._output_file

    def get_discarded(self):
        return self._discarded

    def _load_indels(self):
        indels_file = self._conf.indel_file
        logger.info("Loading indels scores from {}".format(indels_file))
        score_conf = SCORES[self._conf.score]
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


    def _scores(self, e):
        score_conf = SCORES[self._conf.score]

        background_folder = os.path.join(self._conf.cache_dir, self._conf.score, self._conf.feature, e[len(e) - 2:], e.upper())
        background_tsv_file = os.path.join(background_folder, 'background.tsv')

        if not os.path.exists(background_tsv_file):
            silent_mkdir(background_folder)
            self._sampling.create_background_tsv(background_tsv_file, "none", e)

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

    def _indel_score_percentiles(self, element, m, total_indels_score):

        if self._conf.feature != 'cds':
            raise OncodriveFM2Error('The "percentile" indel strategy is only available for the coding sequence (cds) feature')
        background_file = os.path.join(self._sampling.get_cache_dir(element), "background.bin")
        background = np.fromfile(background_file, dtype='float32')
        if len(background) == 0:
            return total_indels_score

        ref = m['REF'].replace("-", "")
        alt = m['ALT'].replace("-", "")
        is_inframe = math.fabs(len(ref) - len(alt)) % 3 == 0

        if is_inframe:
            m['SCORE'] = np.percentile(background, 50)
        else:
            m['SCORE'] = np.percentile(background, 75)

        total_indels_score += 1

        return total_indels_score

    def _indel_score_max(self, m, scores, total_indels_score):
        pos = m['POSITION']
        seq = m['REF'] if m['ALT'] == '-' else m['ALT']
        seq = seq.upper()
        # Skip indels longer than 50
        if len(seq) < 50:
            indels_scores = []
            for pos, n in enumerate(seq, start=pos):
                if n not in "ACTG":
                    break
                values = scores.get(str(pos), [])
                for v in values:
                    if v['alt'] == n:
                        indels_scores.append(v['value'])

            if len(indels_scores) > 0:
                m['SCORE'] = max(indels_scores)
                total_indels_score += 1

        return total_indels_score

    def _indel_score_user_supplied(self, indels, m, total_indels_score):
        if self._conf.score.startswith('cadd') and len(indels) > 0:
            key = '{}-{}'.format(m['CHROMOSOME'], m['POSITION'])
            values = indels.get(key, [])
            for v in values:
                if v['ref'] == m['REF'] and (v['alt'] == m['ALT'] or v['alt'] == '.'):
                    m['SCORE'] = v['value']
                    total_indels_score += 1
                    break
        return total_indels_score

    def _compute_score_means(self, elements, indels, chunk):

        result = []

        i = 1
        total = len(elements)

        logger.info(self._discarded)

        for element, element_muts in elements:

            if i % 100 == 0:
                print("[process-{}] {} of {}".format(chunk, i, total))
            i += 1

            # Skip features without mutations
            if len(element_muts) == 0:
                continue

            # Skip features that were discarded:
            if element in self._discarded:
                continue

            # Read element scores
            scores = self._scores(element)

            # Add scores to the element mutations
            scores_by_sample = {}
            scores_list = []
            # count subs & indels) and subs & indels with score
            total_subs = 0
            total_subs_score = 0
            total_indels = 0
            total_indels_score = 0
            positions = []
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
                    if self._conf.indel_strategy == self._conf.INDEL_STRATEGY_USER_INPUT:
                        total_indels_score = self._indel_score_user_supplied(indels, m, total_indels_score)

                    elif self._conf.indel_strategy == self._conf.INDEL_STRATEGY_MAX:
                        total_indels_score = self._indel_score_max(m, scores, total_indels_score)

                    elif self._conf.indel_strategy == self._conf.INDEL_STRATEGY_PERCENTILES:
                        total_indels_score = self._indel_score_percentiles(element, m, total_indels_score)

                # Update scores
                if 'SCORE' in m:

                    sample = m['SAMPLE']
                    if sample not in scores_by_sample:
                        scores_by_sample[sample] = []

                    scores_by_sample[sample].append(m['SCORE'])
                    scores_list.append(m['SCORE'])
                    positions.append(m['POSITION'])

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
                'symbol': self.gene_conversion.get(element, "UNKNOWN"),
                'scores': scores_list,
                'positions': positions
            }

            result.append((element, item))

        return result


    def _prepare_background(self, pool, elements):
        to_background = []
        signature = self._conf.signature

        for e in elements:

            if not self._sampling.has_regions(e):
                self._discarded.append(e)
                continue

            background_folder = os.path.join(self._conf.cache_dir, self._conf.score, self._conf.feature, e[len(e) - 2:], e.upper())
            background_tsv_file = os.path.join(background_folder, 'background.tsv')
            if not os.path.exists(background_tsv_file):
                if not os.path.exists(os.path.join(background_folder, '.lock')):
                    to_background.append((e, 0))
        if len(to_background) > 0:
            logger.info("Send {} backgrounds".format(len(to_background)))
            if self._sampling.on_cluster:
                prepare_arguments = dict(signature=signature, verbose=True, max_jobs=50)
                return pool.apply(self._sampling.sampling_prepare, (to_background,), prepare_arguments)
            else:
                compute_arguments = ((signature, e, n) for e, n in to_background)
                return pool.starmap(self._sampling.sampling, compute_arguments)
        else:
            return None


    def run(self):

        # Skip if done
        output_folder = self._conf.output_dir
        results_file = os.path.join(output_folder, self._conf.SAMPLING_RESULTS_FILE_NAME.format(self._conf.score))
        if os.path.exists(results_file):
            logger.info("{} - {} - {} - {} [skip]".format(self._conf.project, self._conf.tissue, self._conf.score, self._conf.feature))
            self._output_file = results_file
            return

        # Loading mutations
        mutations_file = os.path.join(self._conf.input, self._conf.MUTS_FILE_NAME)
        if not os.path.exists(mutations_file):
            raise RuntimeError("Mutations file '{}' not found".format(mutations_file))
        with gzip.open(mutations_file, 'rt') as fd:
            elements = json.load(fd)

        if self._conf.testing:
            smallset = sorted(list(elements.keys()))[10:45]
            elements = {k: elements[k] for k in smallset}

        # Initialize
        pool = Pool(self._conf.cores)
        if self._conf.score not in SCORES:
            raise RuntimeError("Unknown score '{}'".format(self._conf.score))

        # Load indels scores if supplied by user
        indels = self._load_indels() if self._conf.indel_strategy == self._conf.INDEL_STRATEGY_USER_INPUT else None

        # Prepare background
        self._prepare_background(pool, elements)

        # Compute elements statistics
        logger.info("Computing score statistics by element")
        results = {}
        elements_all = list(elements.items())
        elements_chunks = np.array_split(elements_all, self._conf.cores)
        compute_arguments = ((chunk, indels, num) for num, chunk in enumerate(elements_chunks, start=1))
        c = 0
        for done in pool.starmap(self._compute_score_means, compute_arguments):
            c += 1
            logger.info("Chunk {} of {} [done]".format(c, self._conf.cores))
            for element, means in done:
                results[element] = means

        # Store results
        logger.info("Store score sampling results")
        silent_mkdir(output_folder)
        with gzip.open(results_file, 'wt') as fd:
            json.dump(results, fd)
            self._output_file = results_file

        logger.info("{} - {} - {} - {} [done]".format(self._conf.project, self._conf.tissue, self._conf.score, self._conf.feature))