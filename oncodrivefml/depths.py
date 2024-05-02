"""
This module contains the methods associated with the
depths that are assigned to the mutations.

The depths are read from a file.
The file must be compressed using bgzip, and then indexed using tabix.
$ bgzip ..../all_samples.depth_per_site.tsv
$ tabix -b 2 -e 2 ..../all_samples.depth_per_site.tsv.gz

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


DepthValue = namedtuple('DepthValue', ['value'])
"""
Tuple that contains the reference, the alteration, the depth value

Parameters:
    value (float): depth value of that substitution
"""

depths_reader = None


class ReaderError(Exception):

    def __init__(self, msg):
        self.message = msg


class ReaderGetError(ReaderError):
    def __init__(self, chr, start, end):
        self.message = 'Error reading chr: {} start: {} end: {}'.format(chr, start, end)

class DepthTabixReader:

    def __init__(self, conf):
        self.file = conf['file']
        self.conf_chr_prefix = conf['chr_prefix']
        self.pos_pos = conf['pos']
        self.depth_pos = conf['depth']
        self.element_pos = conf['element']

    def __enter__(self):
        self.tb = tabix.open(self.file)
        self.index_errors = 0
        self.elements_errors = 0
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.index_errors > 0 or self.elements_errors > 0:
            raise ReaderError('{} index errors and {} discrepancies between the expected and retrieved element'.format(self.index_errors, self.elements_errors))
        return True

    def _read_row(self, row):
        depth = float(row[self.depth_pos])
        pos = None if self.pos_pos is None else int(row[self.pos_pos])
        element = None if self.element_pos is None else row[self.element_pos]
        return (depth, pos), element

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


def init_depths_module(conf):
    global depths_reader
    # TODO add an else case or fix this function
    depths_reader = DepthTabixReader(conf)


class Depths(object):
    """

    Args:
        element (str): element ID
        segments (list): list of the segments associated to the element
        config (dict): configuration

    Attributes:
        depths_by_pos (dict): for each positions get all possible changes

            .. code-block:: python

                    { position:
                        [
                            DepthValue(
                                value
                            )
                        ]
                    }
    """

    def __init__(self, element: str, segments: list, config: dict):

        self.element = element
        self.segments = segments
        
        # depth configuration
        self.conf_file = config['file']
        # self.conf_depth = config['depth']
        self.conf_chr = config['chr']
        self.conf_chr_prefix = config['chr_prefix']
        self.conf_pos = config['pos']
        self.conf_element = config['element']
        self.conf_extra = config['extra']

        # depths to load
        self.depths_by_pos = defaultdict(list)

        # Initialize background depths
        self._load_depths()

    def get_depth_by_position(self, position: int) -> List[DepthValue]:
        """
        Get all DepthValue objects that are associated with that position

        Args:
            position (int): position

        Returns:
            :obj:`list` of :obj:`DepthValue`: list of all DepthValue related to that position

        """
        return self.depths_by_pos.get(position, [])

    def get_all_positions(self) -> List[int]:
        """
        Get all positions in the element

        Returns:
            :obj:`list` of :obj:`int`: list of positions

        """
        return self.depths_by_pos.keys()

    def _load_depths(self):
        """
        For each position get all possible substitutions and for each
        obtains the assigned depth

        Returns:
            dict: for each positions get a list of DepthValue
            (see :attr:`depths_by_pos`)
        """

        try:
            with depths_reader as reader:
                for region in self.segments:
                    try:
                        for row in reader.get(region['CHROMOSOME'], region['START'], region['END'], self.element):
                            depth, pos = row
                            self.depths_by_pos[pos].append(DepthValue(depth))

                    except ReaderError as e:
                        logger.warning(e.message)
                        continue
        except ReaderError as e:
            logger.warning("Reader error: %s. Regions being analysed %s", e.message, self.segments)