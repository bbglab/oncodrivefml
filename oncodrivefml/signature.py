import logging
import os
import mmap
import pandas as pd
from oncodrivefml.load import load_mutations

HG19_LIST = [
    "/projects_bg/bg/soft/intogen_home/gencluster/software/mutsigCV/reffiles/chr_files_hg19",
    os.path.expanduser(os.path.expandvars("$FULL_GENOME_PATH"))
]
HG19 = None
for t in HG19_LIST:
    if os.path.exists(t):
        HG19 = t
        break

HG19_MMAP_FILES = {}

def get_hg19_mmap(chromosome):
    if chromosome not in HG19_MMAP_FILES:
        fd = open(os.path.join(HG19, "chr{0}.txt".format(chromosome)), 'rb')
        HG19_MMAP_FILES[chromosome] = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    return HG19_MMAP_FILES[chromosome]


def get_ref_triplet(chromosome, start):
    mm_file = get_hg19_mmap(chromosome)
    mm_file.seek(start-1)
    return mm_file.read(3).decode().upper()

def get_reference_signature(line):
    return get_ref_triplet(line['CHROMOSOME'], line['POSITION'] - 1)


def get_alternate_signature(line):
    return line['Signature_reference'][0] + line['ALT'] + line['Signature_reference'][2]


def compute_signature(variants_file):
    mutations = pd.DataFrame.from_dict([r for r in load_mutations(variants_file, show_warnings=False) if r['TYPE'] == 'subs'])
    mutations = mutations.groupby(['CHROMOSOME', 'POSITION', 'REF', 'ALT']).count()
    mutations.reset_index(inplace=True)
    mutations['Signature_reference'] = mutations.apply(get_reference_signature, axis=1)
    mutations['Signature_alternate'] = mutations.apply(get_alternate_signature, axis=1)
    result = mutations.groupby(['Signature_reference', 'Signature_alternate']).agg({'SAMPLE': 'sum'})
    result.columns = ['count']
    result['probability'] = result['count'] / result['count'].sum()
    return result.to_dict()['probability']

def load_signature(variants_file, signature_file, signature_field, signature_type):
    signature_dict = None
    if signature_type == "none":
        # We don't use signature
        logging.warning("We are not using any signature")
    elif signature_type == "compute":
        logging.info("Computing signature")
        signature_dict = compute_signature(variants_file)
    else:
        if not os.path.exists(signature_file):
            logging.error("Signature file {} not found.".format(signature_file))
            return -1
        else:
            logging.info("Loading signature")
            signature_probabilities = pd.read_csv(signature_file, sep='\t')
            signature_probabilities.set_index(['Signature_reference', 'Signature_alternate'], inplace=True)
            signature_dict = signature_probabilities.to_dict()[signature_field]
    return signature_dict
