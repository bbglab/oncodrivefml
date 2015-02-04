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
from configuration import SAMPLING_INPUT_FOLDER, TABIX, SCORES, DEFAULT_SAMPLING_SIZE, HG19_DIR, SIGNATURES, COMMAND, Configuration

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

def silent_mkdir(folder):
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


class SignatureSampling(object):
    def __init__(self, conf: Configuration):

        # execution configuration
        self._conf = conf

        # sampling on cluster
        self.on_cluster = False
        try:
            import drmaa
            self.on_cluster = True
        except RuntimeError:
            pass

        # feature regions coordinates
        self._region_elements = None
        self._all_regions = None
        feature_regions_file = os.path.join(SAMPLING_INPUT_FOLDER, 'features_regions', '{}.regions'.format(conf.feature))
        if not os.path.exists(feature_regions_file):
            logger.error("%s - %s - %s - %s - background.tsv [error]", conf.score, conf.signature, conf.feature, conf.element)
            raise RuntimeError("File {} not found".format(feature_regions_file))
        # Load regions at this feature
        self._all_regions = pd.read_csv(feature_regions_file,
                                  names=['chr', 'start', 'stop', 'feature'],
                                  dtype={'chr': object, 'start': np.int, 'stop': np.int, 'feature': object},
                                  sep='\t')

    def _get_ref_triplet(self, chromosome, start):

        # Normalize chromosome
        chromosome = chromosome.replace('chr', '')

        with open(os.path.join(HG19_DIR, "chr{0}.txt".format(chromosome)), 'rb') as hg19:
            hg19.seek(start-1)
            return hg19.read(3).decode().upper()


    def create_background_tsv(self, background_tsv_file, signature, element):

        score = self._conf.score
        feature = self._conf.feature


        if os.path.exists(background_tsv_file):
            logger.info("%s - %s - %s - %s - background.tsv [skip]", score, signature, feature, element)
        else:
            regions = self._all_regions[self._all_regions.feature == element]

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


    def _create_background_bin(self, background_tsv_file, background_bin_file, signature, element):

        score = self._conf.score
        feature = self._conf.feature

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

    def _create_background_signature(self, background_tsv_file, background_signature_bin_file, out_folder, signature, element):

        score = self._conf.score
        feature = self._conf.feature

        # Do nothing if output file is None
        if background_signature_bin_file is None:
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
                    ref_triplet = self._get_ref_triplet(row[score_conf['chr']], int(row[score_conf['pos']]) - 1)
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


    def get_cache_dir(self, element):

        score = self._conf.score
        feature = self._conf.feature
        cache_dir = self._conf.cache_dir

        # Output folder
        if not os.path.exists(cache_dir):
            raise RuntimeError("Output folder not found")

        # Score folder
        out_folder = os.path.join(cache_dir, score)
        silent_mkdir(out_folder)
        # Feature folder
        out_folder = os.path.join(out_folder, feature)
        silent_mkdir(out_folder)
        # Element hash folder
        out_folder = os.path.join(out_folder, element[len(element) - 2:])
        silent_mkdir(out_folder)
        # Element folder
        out_folder = os.path.join(out_folder, element.upper())
        silent_mkdir(out_folder)
        return out_folder

    def run(self, signature, element, num_samples, sampling_size=DEFAULT_SAMPLING_SIZE, verbose=True):

        score = self._conf.score
        feature = self._conf.feature

        # Force num_samples to be an integer
        num_samples = int(num_samples)

        # Configure logging
        logger.setLevel(logging.INFO if verbose else logging.ERROR)

        # Use 'none' if we don't want signature
        signature = 'none' if signature is None else signature

        out_folder = self.get_cache_dir(element)

        # Fail if the folder is lock
        lock_file = os.path.join(out_folder, '.lock')
        if os.path.exists(lock_file):
            logger.error("%s - %s - %s - %s - [error]", score, signature, feature, element)
            raise RuntimeError("Folder {} is locked by .lock file".format(out_folder))

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
                self.create_background_tsv(background_tsv_file, signature, element)

                # Create background scores array file
                self._create_background_bin(background_tsv_file, background_bin_file, signature, element)

                # Create probabilities arrays
                self._create_background_signature(background_tsv_file, background_signature_bin_file, out_folder, signature, element)

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
                p_normalized = None
            else:
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

    def sampling(self, signature, element, num_samples, sampling_size=DEFAULT_SAMPLING_SIZE, verbose=False):
        score = self._conf.score

        if score not in SCORES:
            raise RuntimeError("Unknown score '{}'".format(score))

        values = self.run(signature, element, num_samples, sampling_size=sampling_size, verbose=verbose)

        if values is None:
            return None

        return values[:sampling_size]


    def sampling_prepare(self, values, signature='none', sampling_size=DEFAULT_SAMPLING_SIZE, verbose=False, max_jobs=450):

        # values partitioned by pool.apply

        score = self._conf.score
        feature = self._conf.feature
        cache_dir = self._conf.cache_dir

        if score not in SCORES:
            raise RuntimeError("Unknown score '{}'".format(score))

        # Configure logging
        logger.setLevel(logging.INFO if verbose else logging.ERROR)

        jobs_len = len(values)

        output_dir = tempfile.mkdtemp(prefix="sampling_", dir=os.path.join(cache_dir, '..', 'logs'))
        silent_mkdir(output_dir)

        if signature is None:
            signature = 'none'

        values = sorted(values, key=lambda i: i[1], reverse=True)

        def arguments():
            for a in values:
                yield "{0} {1} {2} {3} {4} {5}".format(score, signature, feature, a[0], a[1], sampling_size)

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
                job_name="sam_prep"
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

    def has_regions(self, element):
        if self._region_elements == None:
            self._region_elements = self._all_regions.feature.tolist()
        return element in self._region_elements

#
# if __name__ == "__main__":
#
#     # Configure the logging
#     logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')
#
#     # Parse the arguments
#     parser = argparse.ArgumentParser()
#     parser.add_argument("score")
#     parser.add_argument("signature")
#     parser.add_argument("feature")
#     parser.add_argument("element")
#     parser.add_argument("num_samples", type=float)
#     parser.add_argument("sampling_size", nargs='?', type=int, default=DEFAULT_SAMPLING_SIZE)
#     parser.add_argument("cache_folder")
#     parser.add_argument("input_folder")
#     args = parser.parse_args()
#
#     # Execute the sampling
#     run(args.score, args.signature, args.feature, args.element, args.num_samples, args.input_folder,  args.cache_folder, sampling_size=args.sampling_size)

