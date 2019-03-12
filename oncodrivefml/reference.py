"""
This module contains information related to the reference genome.
"""

import logging

from collections import  Counter
from bgreference import refseq

from oncodrivefml import __logger_name__

logger = logging.getLogger(__logger_name__)

ref_build = 'hg38'
"""
Build of the Reference Genome
"""

__CB = {"A": "T", "T": "A", "G": "C", "C": "G"}


def change_build(build):
    """
    Modify the default build fo the reference genome

    Args:
        build (str): genome reference build

    """
    global ref_build
    ref_build = build
    logger.info('Using %s as reference genome', ref_build.upper())


def get_build():
    return ref_build


def get_ref(chromosome, start, size=1):
    """
    Gets a sequence from the reference genome

    Args:
        chromosome (str): chromosome
        start (int): start position where to look
        size (int): number of bases to retrieve

    Returns:
        str. Sequence from the reference genome

    """
    return refseq(ref_build, chromosome, start, size)


def get_ref_triplet(chromosome, start):
    """

    Args:
        chromosome (str): chromosome identifier
        start (int): starting position

    Returns:
        str: 3 bases from the reference genome

    """
    return get_ref(chromosome, start, size=3)


__CB = {"A": "T", "T": "A", "G": "C", "C": "G"}


def reverse_complementary_sequence(seq):
    return "".join([__CB[base] if base in __CB else base for base in seq.upper()[::-1]])


def triplets(sequence):
    """

    Args:
        sequence (str): sequence of nucleotides

    Yields:
        str. Triplet

    """
    iterator = iter(sequence)

    n1 = next(iterator)
    n2 = next(iterator)

    for n3 in iterator:
        yield n1 + n2 + n3
        n1 = n2
        n2 = n3


def triplet_counter_executor(elements):
    """
    For a list of regions, get all the triplets present
    in all the segments

    Args:
        elements (:obj:`list` of :obj:`list`): list of lists of segments

    Returns:
        :class:`collections.Counter`. Count of each triplet in the regions

    """
    counts = Counter()
    for element in elements:
        for segment in element:
            chrom = segment['CHROMOSOME']
            start = segment['START']
            end = segment['END']
            seq = refseq(ref_build, chrom, start, end-start+1)
            counts.update(triplets(seq))
    return counts


def is_valid_trinucleotides(trinucleotide):
    """
    Check if a trinucleotide has a nucleotide distinct than A, C, G, T
    Args:
        trinucleotide (str): triplet

    Returns:
        bool.

    """
    for nucleotide in trinucleotide:
        if nucleotide not in __CB.keys():
            return False
    return True


def count_valid_trinucleotides(trinucleotides_dict):
    """
    Count how many trinucleotides are valid

    Args:
        trinucleotides_dict (dict): trinucleotides counts

    Returns:
        int. Valid trinucleotides

    """
    distinct_trinucleotides = set()
    for k in trinucleotides_dict.keys():
        if is_valid_trinucleotides(k):
            distinct_trinucleotides.add(k)
    return len(k)
