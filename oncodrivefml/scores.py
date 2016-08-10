"""
This module contains the methods associated with the
scores that are assigned to the mutations.

The scores are read from a file.
"""


import logging
import tabix

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
        segments (list): list of the segmenst associated to the element
        signature (dict): probabilities {signature_key: { (ref, alt): prob }}
        config (dict): configuration

    Attributes:
        scores_by_pos (dict): for each positions get all possible changes, and for each change the probability
            according to the different signature IDs

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

    def __init__(self, element: str, segments: list, signature: dict, config: dict):

        self.element = element
        self.segments = segments
        self.signature = signature

        # Score configuration
        self.conf_file = config['file']
        self.conf_score = config['score']
        self.conf_chr = config['chr']
        self.conf_chr_prefix = config['chr_prefix']
        self.conf_ref = config['ref']
        self.conf_alt = config['alt']
        self.conf_pos = config['pos']
        self.conf_element = config.get('element', None)
        self.conf_extra = config.get('extra', None)

        # Scores to load
        self.scores_by_pos = defaultdict(list)
        self.missing_signatures = {}

        # Initialize background scores and signatures
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

        Returns:
            dict: for each positions get a list of ScoreValue with all signatures for that triplet
            (see :attr:`scores_by_pos`)
        """
        tb = tabix.open(self.conf_file)#conf_file is the file with the scores

        #Loop through a list of dictionaries from the elements dictionary
        for region in self.segments:
            try:
                #get all rows with certain chromosome and startcompute_muts_statistics and stop
                # between the element start -1 and stop
                for row in tb.query("{}{}".format(self.conf_chr_prefix, region['chrom']), region['start']-1, region['stop']):

                    if self.conf_element is not None:
                        # Check that is the element we want
                        if row[self.conf_element] != self.element:
                            continue

                    value = self._read_score(row)

                    ref = row[self.conf_ref]
                    alt = row[self.conf_alt]
                    pos = int(row[self.conf_pos])

                    ref_triplet = get_ref_triplet(row[self.conf_chr].replace(self.conf_chr_prefix, ''), int(row[self.conf_pos]) - 1)

                    if ref is not None and ref_triplet[1] != ref:
                        logging.warning("Background mismatch at position %d at '%s'", int(row[self.conf_pos]), self.element)

                    # Expand funseq2 dots
                    alts = alt if alt is not None and alt != '.' else 'ACGT'.replace(ref, '')

                    for a in alts:
                        alt_triplet = ref_triplet[0] + a + ref_triplet[2]
                        self.scores_by_pos[pos].append(ScoreValue(ref, a, value, ref_triplet, alt_triplet))

            except tabix.TabixError:
                logging.warning("Tabix error at {}='{}{}:{}-{}'".format(self.element, self.conf_chr_prefix, region['chrom'], region['start']-1, region['stop']))
                continue
