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

from oncodrivefml.signature import get_ref_triplet

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

    def _read_score(self, row: list) -> float:
        """
        Parses a score line and returns the score value

        Args:
            row (list): row from the scores file

        Returns:
            float: score value

        """
        value_str = row[self.conf_score]
        if value_str is None or value_str == '':
            if self.conf_extra is not None:
                value_str = row[self.conf_extra].split(',')
                value = 0
                for val in value_str:
                    elm, vals = val.split(':')
                    if elm == self.element:
                        value = float(vals)
                        break
            else:
                value = 0
        else:
            value = float(value_str)

        return value

    def _load_scores(self):
        """
        For each position get all possible substitutions and for each
        obtatins the assigned score

        Returns:
            dict: for each positions get a list of ScoreValue
            (see :attr:`scores_by_pos`)
        """
        tb = tabix.open(self.conf_file)  # conf_file is the file with the scores

        #Loop through a list of dictionaries from the elements dictionary
        for region in self.segments:
            try:
                # get all rows with certain chromosome and start and stop
                # between the element start -1 and stop
                for row in tb.query("{}{}".format(self.conf_chr_prefix, region['CHROMOSOME']), region['START']-1, region['STOP']):

                    if self.conf_element is not None:
                        # Check that is the element we want
                        if row[self.conf_element] != self.element:
                            continue

                    try:
                        value = self._read_score(row)
                    except IndexError:
                        continue

                    if self.conf_alt is None:
                        alt = None
                    else:
                        alt = row[self.conf_alt]

                    pos = int(row[self.conf_pos])

                    ref_triplet = get_ref_triplet(row[self.conf_chr].replace(self.conf_chr_prefix, ''), int(row[self.conf_pos]) - 1)

                    if self.conf_ref is None:
                        ref = ref_triplet[1]
                    else:
                        ref = row[self.conf_ref]

                    if ref is not None and ref_triplet[1] != ref:
                        logging.warning("Background mismatch at position %d at '%s'", int(row[self.conf_pos]), self.element)

                    # Expand funseq2 dots
                    alts = alt if alt is not None and alt != '.' else 'ACGT'.replace(ref, '')

                    for a in alts:
                        alt_triplet = ref_triplet[0] + a + ref_triplet[2]
                        self.scores_by_pos[pos].append(ScoreValue(ref, a, value, ref_triplet, alt_triplet))


            except tabix.TabixError:
                logging.warning("Tabix error at {}='{}{}:{}-{}'".format(self.element, self.conf_chr_prefix, region['CHROMOSOME'], region['START']-1, region['STOP']))
                continue


    def get_stop_scores(self):
        """
        Get the scores of the stops in a gene that fall in the regions
        being analized
        """
        # TODO add other scores
        # TODO note that is only for coding
        stops = defaultdict(list)
        # stops_file = '/home/iker/Desktop/cds_stop/cds_stops.bgz'
        stops_file = bgdata.get_path('datasets', 'genestops', 'cds')

        tb = tabix.open(stops_file)
        for region in self.segments:
            try:
                for row in tb.query(region['CHROMOSOME'], region['START'] - 1, region['STOP']):
                    pos = int(row[1])
                    ref = row[2]
                    alt = row[3]
                    element_id = row[4]

                    if element_id != self.element:
                        continue

                    stops[pos].append(alt)

            except tabix.TabixError:
                logging.warning(
                    "Tabix error at {}='{}:{}-{}'".format(self.element, region['CHROMOSOME'], region['START'] - 1, region['STOP']))
                continue

        self.stop_scores = []
        if len(stops) > 3:
            # if more than 3 positions have stops we get the stop value from those
            for pos, alts in stops.items():
                for s in self.get_score_by_position(pos):
                    if s.alt in alts:
                        self.stop_scores.append(s.value)
        if len(self.stop_scores) < 3:
            A = 8.9168668946147314
            B = 0.082688007694096191
            all_scores = []
            positions = self.get_all_positions()
            for pos in positions:
                for s in self.get_score_by_position(pos):
                    all_scores.append(s.value)
            mean = np.mean(all_scores)
            self.stop_scores = [A * np.exp(B * mean)]
