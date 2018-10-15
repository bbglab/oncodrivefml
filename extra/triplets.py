
import os
import click
import json
import logging
import pandas as pd
from collections import Counter
from multiprocessing.pool import Pool
from bgreference import refseq, _get_dataset


build = 'hg19'
CHR = 'CHROMOSOME'
START = 'START'
STOP = 'STOP'


CHROMOSOMES = ['chr{}'.format(c) for c in range(1,23)] + ['chrX', 'chrY']


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


def chunkizator(iterable, size=1000):
    """
    Creates chunks from an iterable

    Args:
        iterable:
        size (int): elements in the chunk

    Returns:
        list. Chunk

    """
    s = 0
    chunk = []
    for i in iterable:
        if s == size:
            yield chunk
            chunk = []
            s = 0
        chunk.append(i)
        s += 1
    yield chunk


def element_counter_executor(elements):
    """
    For a list of regions, get all the triplets present
    in all the segments

    Args:
        elements (:obj:`list` of :obj:`list` or :obj:`str`): list of lists of segments or a chomosome

    Returns:
        :class:`collections.Counter`. Count of each triplet in the regions

    """
    counts = Counter()
    for segment in elements:
        chrom = segment[CHR]
        start = segment[START]
        stop = segment[STOP]
        try:
            seq = refseq(build, chrom, start, stop-start+1)
        except (ValueError, RuntimeError):
            logging.warning('Error in ref for CHR: {} positions: {}:{}'.format(chrom, start, stop))
            continue
        counts.update(triplets(seq))
    return counts


def element_counter(elements, cores=None):
    """
    Counts triplets in the elements

    Args:
        elements:
        cores (int): cores to use

    Returns:
        :class:`collections.Counter`. Counts of the triplets in the elements

    """
    if cores is None:
        cores = os.cpu_count()
    counter = Counter()
    pool = Pool(cores)

    for result in pool.imap(element_counter_executor, chunkizator(elements, size=8000)):
        counter.update(result)
    return counter


def chromosome_counter_executor(chr):
    """
    Count triplets in chromosome

    Args:
        chr (str): chromosome ID

    Returns:
        :class:`collections.Counter`. Count of each triplet in the regions

    """
    counts = Counter()
    with open(os.path.join(_get_dataset(build), "{}.txt".format(chr)), 'rt') as fd:
        counts.update(triplets(fd.read().upper()))
    return counts


def chromosome_counter(cores=None):
    if cores is None:
        cores = os.cpu_count()
    counter = Counter()
    pool = Pool(cores)

    for result in pool.imap(chromosome_counter_executor, CHROMOSOMES):
        counter.update(result)
    return counter


def save_counts(counts, file):
    filtered_signature = {k: v for k, v in counts.items() if 'N' not in k}
    with open(file, 'wt') as fd:
        json.dump(filtered_signature, fd)
    logging.info('Saved to {}'.format(file))


def load_counts(file):
    with open(file, 'rt') as fd:
        return json.load(fd)


def load_regions(file):
    df = pd.read_csv(file, sep='\t', header=None, names=[CHR, START, STOP, 'strand', 'gene', 'transcript', 'symbol'])
    regions = df[[CHR, START, STOP]]
    return regions.to_dict(orient='records')


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(short_help='Count in specific regions')
@click.argument('regions_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
@click.option('--cores')
def regions(regions_file, output_file, cores):
    elements = load_regions(regions_file)
    counts = element_counter(elements, cores)
    save_counts(counts, output_file)


@click.command(short_help='Count in whole genome')
@click.argument('output_file', type=click.Path())
@click.option('--cores', type=int)
def genome(output_file, cores):
    """Count trinucleotes in genome and save to file"""
    counts = chromosome_counter(cores)
    save_counts(counts,output_file)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option('-g', '--genome_build', type=click.Choice(['hg19', 'hg38', 'hg18', 'mm10', 'c3h']), help='Reference genome build')
def count(genome_build):
    """
    Get the trinucleotides counts
    """
    global build
    if genome_build is not None:
        build = genome_build
        logging.info('Using genome build {}'.format(build))


count.add_command(genome)
count.add_command(regions)


@click.command(short_help='Substract the counts from two files', context_settings=CONTEXT_SETTINGS)
@click.argument('file1', type=click.Path(exists=True))
@click.argument('file2', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
def substract(file1, file2, output_file):
    """Substract the counts from two files"""
    counts1 = Counter(load_counts(file1))
    counts2 = Counter(load_counts(file2))
    counts1.subtract(counts2)
    save_counts(counts1, output_file)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option('--debug', is_flag=True, help='Enable debugging')
@click.version_option(version='0.1')
def cli(debug):
    if debug:
        fmt = logging.Formatter('%(asctime)s %(message)s', datefmt='%H:%M:%S')
        logging.basicConfig(level=logging.DEBUG, format=fmt)
        logging.debug('Debug mode enabled')
    else:
        logging.basicConfig(level=logging.INFO)


cli.add_command(count)
cli.add_command(substract)


if __name__ == "__main__":
    cli()
