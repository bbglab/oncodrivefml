import gzip
import logging
import os
import mmap
import pickle
import pandas as pd
from oncodrivefml.load import load_mutations
from collections import defaultdict

HG19 = None
HG19_MMAP_FILES = {}
__CB = {"A": "T", "T": "A", "G": "C", "C": "G"}


def get_hg19_mmap(chromosome):
    if chromosome not in HG19_MMAP_FILES:
        fd = open(os.path.join(HG19, "chr{0}.txt".format(chromosome)), 'rb')
        HG19_MMAP_FILES[chromosome] = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    return HG19_MMAP_FILES[chromosome]


def get_ref_triplet(chromosome, start):
    mm_file = get_hg19_mmap(chromosome)
    try:
        mm_file.seek(start-1)
    except ValueError:
        return "ERR"
    return mm_file.read(3).decode().upper()


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


def load_signature(variants_file, signature_file, signature_field, signature_type, signature_name, blacklist=None):

    if signature_file is not None and signature_file.endswith(".pickle.gz"):
        with gzip.open(signature_file, 'rb') as fd:
            return pickle.load(fd)

    signature_dict = None
    if signature_type == "none":
        # We don't use signature
        logging.warning("We are not using any signature")
    elif signature_type == "compute":
        logging.info("Computing global signature")
        signature_dict = compute_signature(variants_file, signature_name, blacklist)
    elif signature_type == "bysample":
        logging.info("Computing signature per sample")
        signature_dict = compute_signature_by_sample(variants_file, blacklist)
    else:
        if not os.path.exists(signature_file):
            logging.error("Signature file {} not found.".format(signature_file))
            return -1
        else:
            logging.info("Loading signature")
            signature_probabilities = pd.read_csv(signature_file, sep='\t')
            signature_probabilities.set_index(['Signature_reference', 'Signature_alternate'], inplace=True)
            signature_dict = {signature_name: signature_probabilities.to_dict()[signature_field]}
    return signature_dict
