import argparse
import glob
import json
import logging
import oncodrivefm2
import signaturesampling as sampling
import os
import pandas as pd
import numpy as np
import gzip

from multiprocessing import Pool


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


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


def _load_mutations(project, tissue):
    if tissue != 'pan':
        mutations_file = os.path.join(sampling.SAMPLING_INPUT_FOLDER, "datasets", project, "{}.txt".format(tissue))
        if not os.path.exists(mutations_file):
            raise RuntimeError("Input dataset not found '{}/{}.txt'".format(project, tissue))
        muts = pd.read_csv(mutations_file, sep='\t',
                           dtype={'CHROMOSOME': object, 'POSITION': np.int, 'REF': object, 'ALT': object,
                                  'SAMPLE': object, 'TYPE': object})
        muts = muts[muts.TYPE.isin(['subs', 'indel'])]
        muts.set_index(['CHROMOSOME', 'POSITION'], inplace=True)
        muts.sortlevel(level=0, axis=0, inplace=True)
    else:
        mutations_files = glob.glob(os.path.join(sampling.SAMPLING_INPUT_FOLDER, "datasets", project, "*.txt"))
        muts = pd.read_csv(mutations_files[0], sep='\t',
                           dtype={'CHROMOSOME': object, 'POSITION': np.int, 'REF': object, 'ALT': object,
                                  'SAMPLE': object, 'TYPE': object})
        for mutations_file in mutations_files[1:]:
            more_muts = pd.read_csv(mutations_file, sep='\t',
                                    dtype={'CHROMOSOME': object, 'POSITION': np.int, 'REF': object, 'ALT': object,
                                           'SAMPLE': object, 'TYPE': object})
            muts = pd.concat([muts, more_muts])
        muts = muts[muts.TYPE.isin(['subs', 'indel'])]
        muts.set_index(['CHROMOSOME', 'POSITION'], inplace=True)
        muts.sortlevel(level=0, axis=0, inplace=True)

    return muts


def _load_features(feature):
    regions_file = os.path.join(sampling.SAMPLING_INPUT_FOLDER, "features_regions", "{}.regions".format(feature))
    if not os.path.exists(regions_file):
        raise RuntimeError("Feature file '{}' not found".format(regions_file))
    regions = pd.read_csv(regions_file, names=['chr', 'start', 'stop', 'feature'],
                          dtype={'chr': object, 'start': np.int, 'stop': np.int, 'feature': object}, sep='\t')
    regions = regions.groupby(['feature'])
    return regions


def _run(cores, project, tissue, feature):

    # Skip if done
    output_folder = os.path.join(sampling.OUTPUT_FOLDER, project, tissue, feature)
    mutations_file = os.path.join(output_folder, oncodrivefm2.Configuration.MUTS_FILE_NAME)
    if os.path.exists(mutations_file):
        logger.info("{} - {} - {} [skip]".format(project, tissue, feature))
        return

    # Initialize pool
    pool = Pool(cores)

    # Loading mutations
    logger.info("Loading mutations")
    muts = _load_mutations(project, tissue)

    # Load features
    logger.info("Loading feature regions")
    regions = _load_features(feature)

    # Aggregate mutations
    logger.info("Aggregating mutations by feature element")
    mutations = {}
    elements_all = [(e, segments) for e, segments in regions]
    elements_chunks = np.array_split(elements_all, cores)
    compute_arguments = ((chunk, muts, num) for num, chunk in enumerate(elements_chunks, start=1))
    c = 0
    for done in pool.starmap(_intersect_mutations, compute_arguments):
        c += 1
        logger.info("Chunk {} of {} [done]".format(c, cores))
        for element, variants in done:
            mutations[element] = variants

    # Store results
    logger.info("Store results")
    sampling.silent_mkdir(output_folder)
    with gzip.open(mutations_file, 'wt') as fd:
        json.dump(mutations, fd)

    logger.info("{} - {} - {} [done]".format(project, tissue, feature))


if __name__ == "__main__":

    # Configure the logging
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')

    # Parse the arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("cores", type=int)
    parser.add_argument("project")
    parser.add_argument("tissue")
    parser.add_argument("feature")

    args = parser.parse_args()

    _run(args.cores, args.project, args.tissue, args.feature)