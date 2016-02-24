import logging
import tabix

from typing import List
from collections import defaultdict, namedtuple
from oncodrivefml.signature import get_ref_triplet

ScoreValue = namedtuple('ScoreValue', ['ref', 'alt', 'value', 'signature'])


class Scores(object):

    def __init__(self, element: str, segments: list, signature: dict, config: dict):
        """

        :param element: The element id
        :param segments: A list with the element segments definition
        :param signature: A dict with the signature probability {'signature_key' : { (ref, alt): probability }}
        :param config: Dictionary with the score configuration
        :return:
        """

        self.element = element
        self.segments = segments
        self.signature = signature
        self.config = config

        self.scores_by_pos = defaultdict(list)
        self.missing_signatures = {}

        # Initialize background scores and signature
        self._load_scores()

    def get_score_by_position(self, position: int) -> List[ScoreValue]:
        """
        Get all the posible scores at the given position

        :param position: Genomic position in the gene
        :return: A list of dicts like this {'ref': 'A', 'alt': 'T', 'value': 3.23, 'signature': {'sample1': 0.1, 'sample2': 0.4}}
        """
        return self.scores_by_pos.get(position, [])

    def _read_score(self, row: list) -> float:
        """
        Parses one score line and returns the score value

        :param row: A row from the score file
        :return: The score parsed
        """
        value_str = row[self.config['score']]
        if value_str is None or value_str == '':
            if 'extra' in self.config:
                value_str = row[self.config['extra']].split(',')
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

        tb = tabix.open(self.config['file'])

        for region in self.segments:
            try:
                for row in tb.query("{}{}".format(self.config['chr_prefix'], region['chrom']), region['start']-1, region['stop']):
                    value = self._read_score(row)

                    ref = row[self.config['ref']] if 'ref' in self.config else None
                    alt = row[self.config['alt']] if 'alt' in self.config else None
                    pos = int(row[self.config['pos']])

                    if self.config.get('element', None) is not None:
                        if row[self.config['element']] != self.element:
                            continue

                    if self.signature is not None:
                        ref_triplet = get_ref_triplet(row[self.config['chr']].replace(self.config['chr_prefix'], ''), int(row[self.config['pos']]) - 1)
                        ref = row[self.config['ref']] if 'ref' in self.config else ref_triplet[1]
                        alt = row[self.config['alt']] if 'alt' in self.config else None

                        if ref is not None and ref_triplet[1] != ref:
                            logging.warning("Background mismatch at position %d at '%s'", int(row[self.config['pos']]), self.element)

                    # Expand funseq2 dots
                    alts = alt if alt is not None and alt != '.' else 'ACGT'.replace(ref, '')

                    for a in alts:

                        pos_signature = {}
                        if self.signature is not None:
                            alt_triplet = ref_triplet[0] + a + ref_triplet[2]
                            try:
                                for k in self.signature.keys():
                                    signature_value = self.signature[k].get((ref_triplet, alt_triplet), 0.0)
                                    pos_signature[k] = signature_value
                            except KeyError:
                                self.missing_signatures[ref_triplet] = alt_triplet

                        self.scores_by_pos[pos].append(ScoreValue(ref, a, value, pos_signature))

            except tabix.TabixError:
                logging.warning("Tabix error at {}='{}{}:{}-{}'".format(self.element, self.config['chr_prefix'], region['chrom'], region['start']-1, region['stop']))
                continue
