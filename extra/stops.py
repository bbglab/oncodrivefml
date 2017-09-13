import click
import logging
import pandas as pd
from bgparsers.readers import elements
from bgreference import refseq

genome_builds = ['hg19', 'hg38', 'hg18', 'mm10', 'c3h']
build = 'hg19'
nucleotides = 'ACTG'

CHR = 'CHROMOSOME'
START = 'START'
STOP = 'STOP'
POS= 'POSITION'
REF = 'REF'
ALT = 'ALT'
GENE_ID = 'GENE'


def generate_vep_mutations(regions_file):
    """
    Generate all possible SNP for all positions in a regions file

    Args:
        regions_file (str): path to a regions file

    Returns:
        :class:`~pandas.DataFrame`.

    """
    muts_row = []
    for element, segments in elements(regions_file).items():
        for segment in segments:
            chr = segment[CHR]
            start = segment[START]
            stop = segment[STOP]
            strand = segment['STRAND']
            if strand == '.':  # unknown strand treated as forward  TODO check is fine
                strand = '+'
            for index, ref in enumerate(refseq(build, chr, start, stop - start + 1)):
                for n in nucleotides:
                    if n != ref:
                        muts_row.append([chr, start+index, start+index, ref+'/'+n, strand])  # TODO check that position is fine
    return pd.DataFrame(muts_row)


def get_ref(row):
    chr = row[CHR]
    pos = int(row[POS])
    return refseq(build, chr, pos)


def get_stops_from_vep(vep_file):
    df = pd.read_csv(vep_file, sep='\t', na_values='-')
    df = df[df['Consequence'].str.contains("stop_gained")]  # keep only stops
    df[[CHR, POS]] = df['Location'].str.split(':', expand=True)
    df.rename(columns={'Allele': ALT, 'Gene': GENE_ID}, inplace=True)
    df = df[[CHR, POS, ALT, GENE_ID]]
    df.drop_duplicates(inplace=True)
    df[REF] = df.apply(get_ref, axis=1)
    return df.sort_values([CHR, POS])


def set_build(b):
    global build
    if b is not None:
        build = b
        logging.info('Using genome build {}'.format(build))


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(short_help='Generate mutations')
@click.argument('regions_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
@click.option('-g', '--genome_build', 'build', type=click.Choice(genome_builds), help='Reference genome build')
def muts4vep(regions_file, output_file, build):
    """Generate a mutations file for VEP with all possible mutations for a given set of regions"""
    set_build(build)
    df = generate_vep_mutations(regions_file)
    df.to_csv(output_file, sep='\t', index=False)
    logging.info('Saving to {}'.format(output_file))


@click.command(short_help='Get stops')
@click.argument('vep_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
@click.option('-g', '--genome_build', 'build', type=click.Choice(genome_builds), help='Reference genome build')
def vep2stops(vep_file, output_file, build):
    """Generate a mutations file with the stops from a VEP output file"""
    set_build(build)
    df = get_stops_from_vep(vep_file)
    df.to_csv(output_file, sep='\t', index=False, header=False)
    logging.info('Saving to {}'.format(output_file))


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


cli.add_command(muts4vep)
cli.add_command(vep2stops)

if __name__ == "__main__":
    cli()
