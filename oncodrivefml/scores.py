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

        # Score configuration
        self.conf_file = config['file']#TODO rename as score config file
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

        # Initialize background scores and signature
        self._load_scores()

    def get_score_by_position(self, position: int) -> List[ScoreValue]:
        """
        Get all the posible scores at the given position

        :param position: Genomic position in the gene
        :return: A list of dicts like this {'ref': 'A', 'alt': 'T', 'value': 3.23, 'signature': {'sample1': 0.1, 'sample2': 0.4}}
        """
        return self.scores_by_pos.get(position, [])

    def get_all_positions(self) -> List[int]:
        """
        Get all possible score positions in the element

        :return: All the
        """
        return self.scores_by_pos.keys()

    def _read_score(self, row: list) -> float:
        """
        Parses one score line and returns the score value

        :param row: A row from the score file
        :return: The score parsed
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
        For each mutation in certain element:
        Looks in the scores file, for the scores associated with that element,
        chromosome, start position and end position.
        Fills scores_by_pos:
        { position : [ (ref, alt, score, { signature :  probability } ]
        compute_muts_statistics
        :return:
        """

        tb = tabix.open(self.conf_file)#conf_file is the file with the scores

        #Loop through a list of dictionaries from the elements dictionary
        for region in self.segments:
            try:
                #get all rows with certain chromosome and startcompute_muts_statistics and stop
                # between the element start -1 and stop
                for row in tb.query("{}{}".format(self.conf_chr_prefix, region['chrom']), region['start']-1, region['stop']):
                    value = self._read_score(row)

                    ref = row[self.conf_ref]
                    alt = row[self.conf_alt]
                    pos = int(row[self.conf_pos])

                    if self.conf_element is not None:
                        #Check that is the element we want
                        #TODO move up to reduce time
                        if row[self.conf_element] != self.element:
                            continue

                    if self.signature is not None:
                        ref_triplet = get_ref_triplet(row[self.conf_chr].replace(self.conf_chr_prefix, ''), int(row[self.conf_pos]) - 1)
                        ref = row[self.conf_ref]
                        alt = row[self.conf_alt]

                        if ref is not None and ref_triplet[1] != ref:
                            logging.warning("Background mismatch at position %d at '%s'", int(row[self.conf_pos]), self.element)

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
                logging.warning("Tabix error at {}='{}{}:{}-{}'".format(self.element, self.conf_chr_prefix, region['chrom'], region['start']-1, region['stop']))
                continue
