import argparse
import csv
import logging
import stat
import subprocess
import os
import tempfile
from bgcore.multiprocess.qmap import QMapExecutor
import pandas as pd
import numpy as np
import shutil
from os.path import join

try:
    import drmaa
    SAMPLING_CACHE_FOLDER = "/data/users/jdeu/projects/fmdrivers/cache"
    OUTPUT_FOLDER = "/data/users/jdeu/projects/fmdrivers/output"
except RuntimeError:
    SAMPLING_CACHE_FOLDER = "/shared/projects/fmdrivers/cache"
    OUTPUT_FOLDER = "/shared/projects/fmdrivers/output"

DEFAULT_SAMPLING_SIZE = 10000
TABIX = "/soft/bio/sequence/tabix-0.2.3/tabix"
INPUT_FOLDER = "/projects_bg/bg/shared/projects/fmdrivers/input"

SCORES = {
    'cadd': {
        'file': join(INPUT_FOLDER, 'scores', 'cadd', 'whole_genome_SNVs.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5, 'element': None
    },
    'cadd_12': {
        'file': join(INPUT_FOLDER, 'scores', 'cadd_12', 'whole_genome_SNVs.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5, 'element': None
    },
    'funseq2': {
        'file': join(INPUT_FOLDER, 'scores', 'funseq2', 'hg19_wg_score.tsv.gz'),
        'chr': 0, 'chr_prefix': 'chr', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'rnasnp': {
        'file': join(INPUT_FOLDER, 'scores', 'rnasnp', 'hg19_rnasnp_scores.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 5, 'element': 6
    },
    'rnasnp_lncrna': {
        'file': join(INPUT_FOLDER, 'scores', 'rnasnp_lncrna', 'hg19_rnasnp_lncRNA_scores.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 5, 'element': 6
    },
    'rnasnp_mirna': {
        'file': join(INPUT_FOLDER, 'scores', 'rnasnp_mirna', 'hg19_rnasnp_preMirna_scores.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 5, 'element': 6
    },
    'tfbs_creation': {
        'file': join(INPUT_FOLDER, 'scores', 'tfbs_creation', 'tfbs_creation.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'tfbs_disruption': {
        'file': join(INPUT_FOLDER, 'scores', 'tfbs_disruption', 'tfbs_disruption.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'mirnats_bladder': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 6, 'element': 5
    },
    'mirnats_blood': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 7, 'element': 5
    },
    'mirnats_brain': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 8, 'element': 5
    },
    'mirnats_breast': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 9, 'element': 5
    },
    'mirnats_colon': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 10, 'element': 5
    },
    'mirnats_heart': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 11, 'element': 5
    },
    'mirnats_kidney': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 12, 'element': 5
    },
    'mirnats_liver': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 13, 'element': 5
    },
    'mirnats_lung': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 14, 'element': 5
    },
    'mirnats_ovary': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 15, 'element': 5
    },
    'mirnats_pancreas': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 16, 'element': 5
    },
    'mirnats_prostate': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 17, 'element': 5
    },
    'mirnats_spleen': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 18, 'element': 5
    },
    'mirnats_stomach': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 19, 'element': 5
    },
    'mirnats_testicle': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 20, 'element': 5
    },
    'mirnats_uterus': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 21, 'element': 5
    },
    'mirnats_all': {
        'file': join(INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 22, 'element': 5
    }
}

SIGNATURES = {
    'pilot': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_pilotproject.tsv',
    'tcga': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_tcga.tsv',
    'intogen': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_intogen.tsv',
    'cll': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_cll.tsv',
    'pwcag': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_pwcag.tsv',
    'simulation': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_simulation.tsv',
    'pilot_simulation_1kb': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_pilotproject.tsv',
    'pilot_simulation_10kb': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_pilotproject.tsv',
    'pilot_simulation_50kb': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_pilotproject.tsv',
}

HG19_DIR = "/projects_bg/bg/soft/intogen_home/gencluster/software/mutsigCV/reffiles/chr_files_hg19"
COMMAND = "python /projects_bg/bg/shared/projects/fmdrivers/scripts/signaturesampling.py"

logger = logging.getLogger(__name__)


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


def _get_ref_triplet(chromosome, start):

    # Normalize chromosome
    chromosome = chromosome.replace('chr', '')

    with open(os.path.join(HG19_DIR, "chr{0}.txt".format(chromosome)), 'rb') as hg19:
        hg19.seek(start-1)
        return hg19.read(3).decode().upper()


def create_background_tsv(score, background_tsv_file, signature, feature, element):
    if os.path.exists(background_tsv_file):
        logger.info("%s - %s - %s - %s - background.tsv [skip]", score, signature, feature, element)
    else:
        feature_regions_file = os.path.join(INPUT_FOLDER, 'features_regions', '{}.regions'.format(feature))
        if not os.path.exists(feature_regions_file):
            logger.error("%s - %s - %s - %s - background.tsv [error]", score, signature, feature, element)
            raise RuntimeError("File {} not found".format(feature_regions_file))

        # Load regions at this feature
        all_regions = pd.read_csv(feature_regions_file,
                                  names=['chr', 'start', 'stop', 'feature'],
                                  dtype={'chr': object, 'start': np.int, 'stop': np.int, 'feature': object},
                                  sep='\t')
        regions = all_regions[all_regions.feature == element]

        score_conf = SCORES[score]

        # Build the command to query with Tabix
        cmd = [TABIX, score_conf['file']]
        cmd += ['{}{}:{}-{}'.format(score_conf['chr_prefix'], region.chr, region.start, region.stop) for idx, region in regions.iterrows()]
        cmd += ['>', background_tsv_file]

        # Execute it
        command = ' '.join(cmd)
        logger.debug("%s", command)
        ret_code = subprocess.call(command, shell=True)

        if ret_code != 0:
            logger.error("%s - %s - %s - %s - background.tsv [error]", score, signature, feature, element)
            raise RuntimeError("Running {}".format(command))

        _group_writable(background_tsv_file)
        logger.info("%s - %s - %s - %s - background.tsv [done]", score, signature, feature, element)


def _create_background_bin(score, background_tsv_file, background_bin_file, signature, feature, element):
    if os.path.exists(background_bin_file):
        logger.info("%s - %s - %s - %s - background.bin [skip]", score, signature, feature, element)
    else:
        with open(background_tsv_file, 'r') as fd:
            reader = csv.reader(fd, delimiter='\t')
            scores = []
            score_conf = SCORES[score]
            for row in reader:
                value = float(row[score_conf['score']])
                alt = row[score_conf['alt']]

                # Skip scores of other elements
                if score_conf['element'] is not None:
                    back_element = row[score_conf['element']]
                    if back_element != element:
                        continue

                # Expand refseq2 dots
                if alt == '.':
                    scores.append(value)
                    scores.append(value)
                    scores.append(value)
                else:
                    scores.append(value)

            np.array(scores, dtype='float32').tofile(background_bin_file, format='float32')
            _group_writable(background_bin_file)

    logger.info("%s - %s - %s - %s - background.bin [done]", score, signature, feature, element)


def _create_background_signature(score, background_tsv_file, background_signature_bin_file, out_folder, signature, feature, element):

    # Do nothing if output file is None
    if background_signature_bin_file is None:
        return

    # Indels don't use signature
    if signature.startswith("indels"):
        return

    if os.path.exists(background_signature_bin_file):
        logger.info("%s - %s - %s - %s - background-%s.bin [skip]", score, signature, feature, element, signature)
    else:

        # Load signature probability
        if '_' not in signature:
            logger.error("%s - %s - %s - %s - background-%s.bin [error]", score, signature, feature, element, signature)
            raise RuntimeError("Invalid signature '{}' name. It must contain a '_' character.".format(signature))

        signature_type = signature.split('_')[0]
        if signature_type not in SIGNATURES:
            logger.error("%s - %s - %s - %s - background-%s.bin [error]", score, signature, feature, element, signature)
            raise RuntimeError("No signature file found for the signature {}".format(signature))

        signature_probabilities = pd.read_csv(SIGNATURES[signature_type], sep='\t')
        signature_probabilities.set_index(['Signature_reference', 'Signature_alternate'], inplace=True)

        probabilities = []
        signature_columns = [f for f in signature_probabilities.columns if f.startswith('Probability_')]
        signatures = ["{}_{}".format(signature_type, c.replace('Probability_', '').lower()) for c in signature_columns]
        score_conf = SCORES[score]
        with open(background_tsv_file, 'r') as fd:
            reader = csv.reader(fd, delimiter='\t')
            for row in reader:
                ref_triplet = _get_ref_triplet(row[score_conf['chr']], int(row[score_conf['pos']]) - 1)
                ref = row[score_conf['ref']]
                alt = row[score_conf['alt']]

                # Skip scores of other elements
                if score_conf['element'] is not None:
                    back_element = row[score_conf['element']]
                    if back_element != element:
                        continue

                if ref_triplet[1] != ref:
                    logger.warning("Background mismatch at position %d at '%s'", int(row[score_conf['pos']]), background_tsv_file)

                # Expand funseq2 dots
                alts = alt if alt != '.' else 'ACGT'.replace(ref, '')

                for a in alts:
                    alt_triplet = ref_triplet[0] + a + ref_triplet[2]
                    try:
                        signature = signature_probabilities.loc[ref_triplet, alt_triplet]
                        probabilities.append(dict(zip(signatures, signature[signature_columns])))
                    except KeyError:
                        logger.warning("Triplet without probability ref: '%s' alt: '%s'", ref_triplet, alt_triplet)
                        probabilities.append(dict(zip(signatures, [0]*len(signature_columns))))

        for s in signatures:
            c_array = np.array([p[s] for p in probabilities], dtype='float32')

            if len(c_array) != len(probabilities):
                raise RuntimeError("Missing signature probabilities at {}".format(s))

            c_file = os.path.join(out_folder, 'background-{}.bin'.format(s))
            c_array.tofile(c_file, format='float32')
            _group_writable(c_file)

            logger.info("%s - %s - %s - %s - background-%s.bin [done]", score, s, feature, element, s)


def _run(score, signature, feature, element, num_samples, sampling_size=DEFAULT_SAMPLING_SIZE, verbose=True):

    # Force num_samples to be an integer
    num_samples = int(num_samples)

    # Configure logging
    logger.setLevel(logging.INFO if verbose else logging.ERROR)

    # Use 'none' if we don't want signature
    signature = 'none' if signature is None else signature

    # Output folder
    if not os.path.exists(SAMPLING_CACHE_FOLDER):
        raise RuntimeError("Output folder not found")

    # Score folder
    out_folder = os.path.join(SAMPLING_CACHE_FOLDER, score)
    _silent_mkdir(out_folder)

    # Feature folder
    out_folder = os.path.join(out_folder, feature)
    _silent_mkdir(out_folder)

    # Element hash folder
    out_folder = os.path.join(out_folder, element[len(element)-2:])
    _silent_mkdir(out_folder)

    # Element folder
    out_folder = os.path.join(out_folder, element.upper())
    _silent_mkdir(out_folder)

    # Fail if the folder is lock
    lock_file = os.path.join(out_folder, '.lock')
    if os.path.exists(lock_file):
        logger.error("%s - %s - %s - %s - [error]", score, signature, feature, element)
        raise RuntimeError("Folder {} is lock".format(out_folder))

    # Output file
    sampling_file = os.path.join(out_folder, 'sampling-{}-{}.bin'.format(signature, num_samples))

    # Check if it's already precomputed
    values = None
    if os.path.exists(sampling_file):
        values = np.fromfile(sampling_file, dtype='float32')
        # Check that we have enough precomputed samples
        if len(values) >= sampling_size:
            logger.info("%s - %s - %s - sampling-%s-%s.bin [skip]", signature, feature, element, signature, num_samples)
            return values

    # Files that must exists to be able to do the sampling
    background_tsv_file = os.path.join(out_folder, 'background.tsv')
    background_bin_file = os.path.join(out_folder, 'background.bin')
    background_signature_bin_file = os.path.join(out_folder, 'background-{}.bin'.format(signature.lower())) if signature != "none" else None

    # Check if all background files needed exists, otherwise compute them
    if not(os.path.exists(background_tsv_file) and os.path.exists(background_bin_file) and (signature == "none" or os.path.exists(background_signature_bin_file))):
        try:
            # Lock the folder
            open(lock_file, 'a').close()

            # Create background file for the feature with tabix
            create_background_tsv(score, background_tsv_file, signature, feature, element)

            # Create background scores array file
            _create_background_bin(score, background_tsv_file, background_bin_file, signature, feature, element)

            # Create probabilities arrays
            _create_background_signature(score, background_tsv_file, background_signature_bin_file, out_folder, signature, feature, element)

        finally:
            # Unlock the folder
            os.remove(lock_file)

    # Sampling
    result = None
    if num_samples > 0:
        background = np.fromfile(background_bin_file, dtype='float32')

        if len(background) == 0:
            return None

        if signature == 'none':

            # Subs sampling without signature
            p_normalized = None
        elif signature.startswith("indels"):

            # Indels sampling
            indels_size = int(signature.split("_")[1]) * 3

            # Compute a new background that is the maximum score of all indels_size posible windows
            background = np.array([np.max(background[i:i+indels_size]) for i in range(len(background))])

        else:
            # Subs sampling with signature
            probabilities = np.fromfile(background_signature_bin_file, dtype='float32')
            p_normalized = probabilities / sum(probabilities)


            # Force the float32 array to sum 1 when using float64
            p_normalized[0] += (1 - sum(p_normalized))
            for i in range(0,len(p_normalized)-1):
                if p_normalized[i] < 0:
                    p_normalized[i+1] += p_normalized[i]
                    p_normalized[i] = 0
                else:
                    break

        to_pick = sampling_size - len(values) if values is not None else sampling_size
        if num_samples == 1:
            result = np.array(
                np.random.choice(background, size=to_pick, p=p_normalized, replace=True),
                dtype='float32'
            )
        else:
            # TODO do this in parallel if num_samples * sampling_size is bigger than 10 millions
            result = np.array(
                [np.mean(np.random.choice(background, size=num_samples, p=p_normalized, replace=False)) for a in range(to_pick)],
                dtype='float32'
            )

        # Concat results with previous sampling
        if values is not None:
            result = np.concatenate((values, result), axis=0)

        # Store the sampling
        result.tofile(sampling_file, format='float32')
        _group_writable(sampling_file)

        logger.info("%s - %s - %s - sampling-%s-%s.bin [done]", signature, feature, element, signature, num_samples)

    return result


def fitting(score, feature):

    if score not in SCORES:
        raise RuntimeError("Unknown score '{}'".format(score))

    fitting_file = os.path.join(SAMPLING_CACHE_FOLDER, score, feature, 'fitting.txt')
    return pd.read_csv(fitting_file, sep='\t')


def sampling(score, signature, feature, element, num_samples, sampling_size=DEFAULT_SAMPLING_SIZE, verbose=False):

    if score not in SCORES:
        raise RuntimeError("Unknown score '{}'".format(score))

    values = _run(score, signature, feature, element, num_samples, sampling_size=sampling_size, verbose=verbose)

    if values is None:
        return None

    return values[:sampling_size]


def sampling_prepare(values, feature='cds', signature='none', score='cadd', sampling_size=DEFAULT_SAMPLING_SIZE, verbose=False, max_jobs=200):

    if score not in SCORES:
        raise RuntimeError("Unknown score '{}'".format(score))

    # Configure logging
    logger.setLevel(logging.INFO if verbose else logging.ERROR)

    jobs_len = len(values)

    output_dir = tempfile.mkdtemp(prefix="sampling_", dir=os.path.join(SAMPLING_CACHE_FOLDER, '..', 'logs'))
    _silent_mkdir(output_dir)

    if signature is None:
        signature = 'none'

    values = sorted(values, key=lambda i: i[1], reverse=True)

    def arguments():
        for a in values:
            m_signature = signature if len(a) == 2 else a[2]
            yield "{0} {1} {2} {3} {4} {5}".format(score, m_signature, feature, a[0], a[1], sampling_size)

    retry = 1
    jobs_done = 0
    jobs_fail = 0
    while retry <= 5:

        # Create the qmap executor
        executor = QMapExecutor(
            ['normal', 'long', 'short-high', 'short-low', 'bigmem'],
            min(jobs_len, max_jobs),
            0,
            output_folder=output_dir,
            interactive=False,
            adaptative=False
        )

        # Run all
        jobs_done, jobs_fail, jobs_skip = executor.run(
            COMMAND,
            arguments(),
            len(values),
            job_name="sp"
        )

        # Close the executor
        executor.exit()

        if jobs_fail == 0:
            shutil.rmtree(output_dir)
            return jobs_done, jobs_fail
        else:
            retry += 1
            logger.info("Some jobs fail, retry {} of maximum 5".format(retry))

    logging.error("%d jobs fail. Check the logs at '%s'.", jobs_fail, output_dir)

    return jobs_done, jobs_fail


if __name__ == "__main__":

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')

    # Parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("score")
    parser.add_argument("signature")
    parser.add_argument("feature")
    parser.add_argument("element")
    parser.add_argument("num_samples", type=float)
    parser.add_argument("sampling_size", nargs='?', type=int, default=DEFAULT_SAMPLING_SIZE)
    args = parser.parse_args()

    # Execute the sampling
    _run(args.score, args.signature, args.feature, args.element, args.num_samples, sampling_size=args.sampling_size)

