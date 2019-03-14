"""
This module contains the methods associated with the
scores that are assigned to the mutations.

The scores are read from a file.


Information about the stop scores.

As of December 2016, we have only measured
the stops using CADD1.0.

The stops of a gene retrieved only if there are
ast least 3 stops in the regions being analysed.
If not, a formula is applied to derived the
value of the stops from the rest of the
values.

.. note::

    This formula was obtained using the CADD scores
    of the coding regions. Using a different regions
    or scores files will make the function to return
    totally nonsense values.

"""
import logging
import gzip
import mmap
import json
import os
import struct
from collections import defaultdict, namedtuple
from typing import List

import bgdata
import numpy as np
import tabix

from oncodrivefml import __logger_name__
from oncodrivefml.reference import get_ref_triplet, get_build

logger = logging.getLogger(__logger_name__)


ScoreValue = namedtuple('ScoreValue', ['ref', 'alt', 'value', 'change'])
"""
Tuple that contains the reference, the alteration, the score value and the triplets

Parameters:
    ref (str): reference base
    alt (str): altered base
    value (float): score value of that substitution
    change (str): reference triplet > alt (e.g. ACG>T)
"""


def null(x):
    return x

stop_function = null
min_stops = 1
stops_file = None
scores_reader = None


class ReaderError(Exception):

    def __init__(self, msg):
        self.message = msg


class ReaderGetError(ReaderError):

    def __init__(self, chr, start, end):
        self.message = 'Error reading chr: {} start: {} end: {}'.format(chr, start, end)


class PackScoresReader:

    STRUCT_SIZE = 6
    SCORE_ALT = {'T': 'ACG', 'A': 'CGT', 'C': 'AGT', 'G': 'ACT'}
    BIT_TO_REF = {
        (0, 0, 0): '?',
        (0, 0, 1): 'T',
        (0, 1, 0): 'A',
        (0, 1, 1): 'C',
        (1, 0, 0): 'G'
    }
    SCORE_ORDER = {
        'T': {'A': 0, 'C': 1, 'G': 2},
        'A': {'C': 0, 'G': 1, 'T': 2},
        'C': {'A': 0, 'G': 1, 'T': 2},
        'G': {'A': 0, 'C': 1, 'T': 2}
    }

    def __init__(self, conf):

        self.file_path = os.path.join(conf['file'], 'whole_genome_SNVs.fml')

        with gzip.open('{}.idx'.format(self.file_path), 'rt') as fd:
            index = json.load(fd)

        self.scores = index['scores']
        self.metadata = index['metadata']
        self.fd = None

    def unpack(self, block):
        values = struct.unpack('3H', block)
        bits = tuple([v % 2 for v in values])
        values = [v // 2 for v in values]
        ref = PackScoresReader.BIT_TO_REF[bits]
        return ref, values

    def __enter__(self):
        self.fd = open(self.file_path, 'rb')
        self.mm = mmap.mmap(self.fd.fileno(), 0, access=mmap.ACCESS_READ)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.mm.close()
        self.fd.close()

    def get(self, chromosome, start, stop, *args, **kwargs):

        if chromosome not in self.metadata['positions']:
            raise ReaderError("Chromosome '{}' not found".format(chromosome))

        c_ini, c_end, c_start, c_stop = self.metadata['positions'][chromosome]

        q_ini = c_ini + (start - c_start)
        q_end = q_ini + (stop - start)

        # Check boundaries
        q_ini = c_ini if q_ini < c_ini else q_ini
        q_end = c_end if q_end > c_end else q_end
        if q_end < q_ini:
            raise ReaderError("Position out of boundaries. Chr:{} start:{} stop:{}".format(chromosome, start, stop))

        b_now = q_ini
        p_now = start

        self.mm.seek(b_now*PackScoresReader.STRUCT_SIZE)
        while b_now <= q_end:
            ref, values = self.unpack(self.mm.read(PackScoresReader.STRUCT_SIZE))

            if ref != '?':
                yield self.scores[values[0]], ref, PackScoresReader.SCORE_ALT[ref][0], p_now
                yield self.scores[values[1]], ref, PackScoresReader.SCORE_ALT[ref][1], p_now
                yield self.scores[values[2]], ref, PackScoresReader.SCORE_ALT[ref][2], p_now

            p_now += 1
            b_now += 1


class ScoresTabixReader:

    def __init__(self, conf):
        self.file = conf['file']
        self.conf_chr_prefix = conf['chr_prefix']
        self.ref_pos = conf['ref']
        self.alt_pos = conf['alt']
        self.pos_pos = conf['pos']
        self.score_pos = conf['score']
        self.element_pos = conf['element']

    def __enter__(self):
        self.tb = tabix.open(self.file)
        self.index_errors = 0
        self.elements_errors = 0
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.index_errors > 0 or self.elements_errors > 0:
            raise ReaderError('{} index errors and {} discrepancies between the expected and retreived element'.format(self.index_errors, self.elements_errors))
        return True

    def _read_row(self, row):
        score = float(row[self.score_pos])
        ref = None if self.ref_pos is None else row[self.ref_pos]
        alt = None if self.alt_pos is None else row[self.alt_pos]
        pos = None if self.pos_pos is None else int(row[self.pos_pos])
        element = None if self.element_pos is None else row[self.element_pos]

        return (score, ref, alt, pos), element

    def get(self, chromosome, start, stop, element=None):
        try:
            for row in self.tb.query("{}{}".format(self.conf_chr_prefix, chromosome), start, stop):
                try:
                    r = self._read_row(row)
                except IndexError:
                    self.index_errors += 1
                    continue
                else:
                    if self.element_pos is not None and element is not None and r[1] != element:
                        self.elements_errors += 1
                        continue
                    yield r[0]
        except tabix.TabixError:
            raise ReaderGetError(chromosome, start, stop)


def init_scores_module(conf, stops_required=False):
    global min_stops, stops_file, scores_reader

    if stops_required:
        stops_file = bgdata.get_path('datasets', 'genestops', get_build())

        min_stops = conf.get('minimum_number_of_stops', min_stops)
        logger.debug('Below %d stops in the element the function for stops will be used', min_stops)
        if conf['mean_to_stop_function'] is not None:
            exec("def stop_function(x): return {}".format(conf['mean_to_stop_function']), globals())
        else:
            logger.warning('You have not provided any function for computing the stops')

    if conf['format'] == 'tabix':
        scores_reader = ScoresTabixReader(conf)
    elif conf['format'] == 'pack':
        scores_reader = PackScoresReader(conf)


class Scores(object):
    """

    Args:
        element (str): element ID
        segments (list): list of the segments associated to the element
        config (dict): configuration

    Attributes:
        scores_by_pos (dict): for each positions get all possible changes, and for each change the triplets

            .. code-block:: python

                    { position:
                        [
                            ScoreValue(
                                ref,
                                alt_1,
                                value,
                                change
                            ),
                            ScoreValue(
                                ref,
                                alt_2,
                                value,
                                change
                            ),
                            ScoreValue(
                                ref,
                                alt_3,
                                value,
                                change
                            )
                        ]
                    }
    """

    def __init__(self, element: str, segments: list, config: dict):

        self.element = element
        self.segments = segments

        # Score configuration
        self.conf_file = config['file']
        self.conf_score = config['score']
        self.conf_chr = config['chr']
        self.conf_chr_prefix = config['chr_prefix']
        self.conf_ref = config['ref']
        self.conf_alt = config['alt']
        self.conf_pos = config['pos']
        self.conf_element = config['element']
        self.conf_extra = config['extra']

        # Scores to load
        self.scores_by_pos = defaultdict(list)

        # Initialize background scores
        self._load_scores()

    def get_score_by_position(self, position: int) -> List[ScoreValue]:
        """
        Get all ScoreValue objects that are asocated with that position

        Args:
            position (int): position

        Returns:
            :obj:`list` of :obj:`ScoreValue`: list of all ScoreValue related to that positon

        """
        return self.scores_by_pos.get(position, [])

    def get_all_positions(self) -> List[int]:
        """
        Get all positions in the element

        Returns:
            :obj:`list` of :obj:`int`: list of positions

        """
        return self.scores_by_pos.keys()

    def _load_scores(self):
        """
        For each position get all possible substitutions and for each
        obtatins the assigned score

        Returns:
            dict: for each positions get a list of ScoreValue
            (see :attr:`scores_by_pos`)
        """

        try:
            with scores_reader as reader:
                for region in self.segments:
                    try:
                        for row in reader.get(region['CHROMOSOME'], region['START'], region['END'], self.element):
                            score, ref, alt, pos = row
                            ref_triplet = get_ref_triplet(region['CHROMOSOME'], pos - 1)
                            ref = ref_triplet[1] if ref is None else ref

                            if ref_triplet[1] != ref:
                                logger.warning("Background mismatch at position %d at '%s'", pos, self.element)

                            # Expand funseq2 dots
                            alts = alt if alt is not None and alt != '.' else 'ACGT'.replace(ref, '')

                            for a in alts:
                                self.scores_by_pos[pos].append(ScoreValue(ref, a, score, ref_triplet+'>'+a))

                    except ReaderError as e:
                        logger.warning(e.message)
                        continue
        except ReaderError as e:
            logger.warning("Reader error: %s. Regions being analysed %s", e.message, self.segments)

    def get_stop_scores(self):
        """
        Get the scores of the stops in a gene that fall in the regions
        being analyzed
        """
        stops = defaultdict(list)

        tb = tabix.open(stops_file)
        for region in self.segments:
            try:
                for row in tb.query(region['CHROMOSOME'], region['START'], region['END']):
                    pos = int(row[1])
                    ref = row[2]
                    alt = row[3]
                    # Here we can add a check for the elements

                    stops[pos].append(alt)

            except tabix.TabixError:
                logger.warning(
                    "Tabix error at %s='%s:%d-%d'", self.element, region['CHROMOSOME'], region['START'], region['END'])
                continue

        self.stop_scores = []
        for pos, alts in stops.items():
            for s in self.get_score_by_position(pos):
                if s.alt in alts:
                    self.stop_scores.append(s.value)
        if len(self.stop_scores) < min_stops:
            all_scores = []
            positions = self.get_all_positions()
            for pos in positions:
                for s in self.get_score_by_position(pos):
                    all_scores.append(s.value)
            self.stop_scores = [stop_function(np.mean(all_scores))]
