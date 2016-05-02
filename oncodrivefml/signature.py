import gzip
import logging
import os
import mmap
import pickle

import bgdata
import pandas as pd
from os.path import join, exists
from oncodrivefml.load import load_mutations
from collections import defaultdict

HG19 = None
HG19_MMAP_FILES = {}
__CB = {"A": "T", "T": "A", "G": "C", "C": "G"}


def get_hg19_dataset():
    global HG19

    if HG19 is None:
        HG19 = bgdata.get_path('datasets', 'genomereference', 'hg19')

    return HG19


def get_hg19_mmap(chromosome):
    if chromosome not in HG19_MMAP_FILES:
        fd = open(join(get_hg19_dataset(), "chr{0}.txt".format(chromosome)), 'rb')
        HG19_MMAP_FILES[chromosome] = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    return HG19_MMAP_FILES[chromosome]


def get_ref_triplet(chromosome, start):
    mm_file = get_hg19_mmap(chromosome)
    try:
        mm_file.seek(start-1)
    except ValueError:
        return "ERR"
    return mm_file.read(3).decode().upper()


def get_ref(chromosome, start, size=1):
    mm_file = get_hg19_mmap(chromosome)
    try:
        mm_file.seek(start - 1)
    except ValueError:
        return "ERR"
    return mm_file.read(size).decode().upper()


def get_reference_signature(line):
    return get_ref_triplet(line['CHROMOSOME'], line['POSITION'] - 1)


def get_alternate_signature(line):
    return line['Signature_reference'][0] + line['ALT'] + line['Signature_reference'][2]


def complementary_sequence(seq):
    return "".join([__CB[base] if base in __CB else base for base in seq.upper()])


def collapse_complementaries(signature):
    comp_sig = defaultdict(int)
    for k, v in signature.items():
        comp_sig[k] += v
        comp_k = (complementary_sequence(k[0]), complementary_sequence(k[1]))
        comp_sig[comp_k] += v
    return comp_sig


def signature_probability(signature_counts):
    total = sum([v for v in signature_counts.values()])
    return {k: v/total for k, v in signature_counts.items()}


def compute_signature(variants_file, signature_name, blacklist):
    signature_count = defaultdict(lambda: defaultdict(int))
    for mut in load_mutations(variants_file, signature=signature_name, show_warnings=False, blacklist=blacklist):
        if mut['TYPE'] != 'subs':
            continue

        signature_ref = get_ref_triplet(mut['CHROMOSOME'], mut['POSITION'] - 1)
        signature_alt = signature_ref[0] + mut['ALT'] + signature_ref[2]

        signature_count[mut['SIGNATURE']][(signature_ref, signature_alt)] += 1

    signature = {}
    for k, v in signature_count.items():
        signature[k] = signature_probability(v)

    return signature


def compute_signature_by_sample(variants_file, blacklist, collapse=True):
    signature_count = defaultdict(lambda: defaultdict(int))
    for mut in load_mutations(variants_file, show_warnings=False, blacklist=blacklist):
        if mut['TYPE'] != 'subs':
            continue

        signature_ref = get_ref_triplet(mut['CHROMOSOME'], mut['POSITION'] - 1)
        signature_alt = signature_ref[0] + mut['ALT'] + signature_ref[2]

        signature_count[mut['SAMPLE']][(signature_ref, signature_alt)] += 1

    signature = {}
    for k, v in signature_count.items():
        if collapse:
            signature[k] = signature_probability(collapse_complementaries(v))
        else:
            signature[k] = signature_probability(v)

    return signature


def load_signature(variants_file, signature_config, blacklist=None, signature_name="none"):

    method = signature_config['method']
    path = signature_config.get('path', None)
    column_ref = signature_config.get('column_ref', None)
    column_alt = signature_config.get('column_alt', None)
    column_probability = signature_config.get('column_probability', None)

    if path is not None and path.endswith(".pickle.gz"):
        with gzip.open(path, 'rb') as fd:
            return pickle.load(fd)

    signature_dict = None
    if method == "none":
        # We don't use signature
        logging.warning("We are not using any signature")

    elif method == "full" or method == "complement":

        signature_dict_precomputed = variants_file + "_signature_full.pickle.gz"
        if exists(signature_dict_precomputed):
            logging.info("Using precomputed signature")
            with gzip.open(signature_dict_precomputed, 'rb') as fd:
                signature_dict = pickle.load(fd)
        else:
            logging.info("Computing full global signature")
            signature_dict = compute_signature(variants_file, signature_name, blacklist)
            try:
                # Try to store as precomputed
                with gzip.open(signature_dict_precomputed, 'wb') as fd:
                    pickle.dump(signature_dict, fd)
            except OSError:
                logging.debug("Imposible to write precomputed full signature here: {}".format(signature_dict_precomputed))

        if method == "complement":
            signature_dict = collapse_complementaries(signature_dict)

    elif method == "bysample":
        signature_dict_precomputed = variants_file + "_signature_bysample.pickle.gz"
        if exists(signature_dict_precomputed):
            logging.info("Using precomputed per sample signature")
            with gzip.open(signature_dict_precomputed, 'rb') as fd:
                signature_dict = pickle.load(fd)
        else:
            logging.info("Computing signature per sample")
            signature_dict = compute_signature_by_sample(variants_file, blacklist)
            try:
                # Try to store as precomputed
                with gzip.open(signature_dict_precomputed, 'wb') as fd:
                    pickle.dump(signature_dict, fd)
            except OSError:
                logging.debug("Imposible to write precomputed bysample signature here: {}".format(signature_dict_precomputed))

    elif method == "file":
        if not os.path.exists(path):
            logging.error("Signature file {} not found.".format(path))
            return -1
        else:
            logging.info("Loading signature")
            signature_probabilities = pd.read_csv(path, sep='\t')
            signature_probabilities.set_index([column_ref, column_alt], inplace=True)
            signature_dict = {signature_name: signature_probabilities.to_dict()[column_probability]}
    return signature_dict
