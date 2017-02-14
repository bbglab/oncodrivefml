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
import tabix
import bgdata
import numpy as np
from typing import List
from collections import defaultdict, namedtuple

from oncodrivefml.signature import get_ref_triplet, get_build

ScoreValue = namedtuple('ScoreValue', ['ref', 'alt', 'value', 'ref_triplet', 'alt_triplet'])
"""
Tuple that contains the reference, the alteration, the score value and the triplets

Parameters:
    ref (str): reference base
    alt (str): altered base
    value (float): score value of that substitution
    ref_triplet (str): reference triplet
    alt_triplet (str): altered triplet
"""


def null(x):
    return x

stop_function = null
min_stops = 3
stops_file = None
scores_reader = None


class ReaderError(Exception):

    def __init__(self, msg):
        self.message = msg


class ReaderGetError(ReaderError):

    def __init__(self, chr, start, stop):
        self.message = 'Error reading chr: {} start: {} stop: {}'.format(chr, start, stop)


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
            for row in self.tb.query("{}{}".format(self.conf_chr_prefix, chromosome), start - 1, stop):
                try:
                    r = self._read_row(row)
                except IndexError:
                    self.index_errors += 1
                    continue
                else:
                    # TODO is this check useful???
                    if self.element_pos is not None and r[1] != element:
                        self.elements_errors += 1
                        continue
                    yield r[0]
        except tabix.TabixError:
            raise ReaderGetError(chromosome, start, stop)


def init_scores_module(conf):
    global stop_function, min_stops, stops_file, scores_reader

    min_stops = conf.get('limit', min_stops)
    logging.debug('Below {} stops in the element the function for stops will be used'.format(min_stops))
    if 'function' in conf:
        exec("def stops_function(x): return {}".format(conf['function']), globals())
        stop_function = stops_function
    else:
        logging.warning('You have not provided any function for computing the stops')
    # stops_file = bgdata.get_path('datasets', 'genestops', get_build())  # TODO what to do in case of error
    stops_file = bgdata.get_path('datasets', 'genestops', 'cds')
    scores_reader = ScoresTabixReader(conf)


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
                                ref_triplet,
                                alt_triple
                            ),
                            ScoreValue(
                                ref,
                                alt_2,
                                value,
                                ref_triplet,
                                alt_triple
                            ),
                            ScoreValue(
                                ref,
                                alt_3,
                                value,
                                ref_triplet,
                                alt_triple
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
                        for row in reader.get(region['CHROMOSOME'], region['START'], region['STOP'], self.element):
                            score, ref, alt, pos = row
                            ref_triplet = get_ref_triplet(region['CHROMOSOME'], pos - 1)
                            ref = ref_triplet[1] if ref is None else ref

                            if ref_triplet[1] != ref:
                                logging.warning("Background mismatch at position %d at '%s'", pos, self.element)

                            # Expand funseq2 dots
                            alts = alt if alt is not None and alt != '.' else 'ACGT'.replace(ref, '')

                            for a in alts:
                                alt_triplet = ref_triplet[0] + a + ref_triplet[2]
                                self.scores_by_pos[pos].append(ScoreValue(ref, a, score, ref_triplet, alt_triplet))

                    except ReaderError as e:
                        logging.warning(e.message)
                        continue
        except ReaderError as e:
            logging.warning("Reader error: {}. Regions being analysed {}".format(e.message, self.segments))

    def get_stop_scores(self):
        """
        Get the scores of the stops in a gene that fall in the regions
        being analyzed
        """
        stops = defaultdict(list)

        tb = tabix.open(stops_file)
        for region in self.segments:
            try:
                for row in tb.query(region['CHROMOSOME'], region['START'] - 1, region['STOP']):
                    pos = int(row[1])
                    ref = row[2]
                    alt = row[3]

                    # TODO remove this check?
                    element_id = row[4]
                    if element_id != self.element:
                        continue

                    stops[pos].append(alt)

            except tabix.TabixError:
                logging.warning(
                    "Tabix error at {}='{}:{}-{}'".format(self.element, region['CHROMOSOME'], region['START'] - 1, region['STOP']))
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
            return [stop_function(np.mean(all_scores))]
