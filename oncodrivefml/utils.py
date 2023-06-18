"""
This module contains some useful methods
"""

import logging
from collections import defaultdict
from datetime import datetime
from ago import human
from os.path import exists

from oncodrivefml import __logger_name__

logger = logging.getLogger(__logger_name__)


def executor_run(executor):
    """
    Method to call the run method

    Args:
        executor (:class:`~oncodrivefml.executors.bymutation.ElementExecutor`):

    Returns:
        :meth:`~oncodrivefml.executors.bymutation.ElementExecutor.run`

    """
    return executor.run()


def defaultdict_list():
    """
    Shortcut

    Returns:
        :class:`~collections.defaultdict` of :obj:`list`

    """
    return defaultdict(list)


def loop_logging(iterable, size=None, step=1):
    """
    Loop through an iterable object displaying messages
    using :func:`~logging.info`

    Args:
        iterable:
        size (int): Defaults to None.
        step (int): Defaults to 1.

    Yields:
        The iterable element

    """

    if size is None:
        size = len(iterable)

    i = 0
    start_time = datetime.now()
    for i, value in enumerate(iterable):
        if i % step == 0:
            logger.info("[%d of %d]", i+1, size)
        yield value
    logger.info("[%d of %d]", i+1, size)
    logger.debug("Time: %s", human(start_time))


def exists_path(path):
    return False if path is None else exists(path)




def load_depths(depths_filename, chr_prefix):
    """
    Function to load the depths file into a Python dictionary
    {
        1 : {POS1: DEPTH1, POS2: DEPTH2},
        2 : {POS1: DEPTH1, POS2: DEPTH2},
        3 : {POS1: DEPTH1, POS2: DEPTH2},
        ...
    }
    the chr_prefix is removed from the beginning of the chromosome string
    """
    depths_dict = {}
    with open(depths_filename, 'r') as f:
        for line in f:
            # print(line)
            chromosome, position, value = line.strip().split("\t")
            chromosome = chromosome.lstrip(chr_prefix)
            # print(chromosome, position, value)
            if chromosome not in depths_dict.keys():
                depths_dict[chromosome] = {}
            
            depths_dict[chromosome][position] = int(value)
    return depths_dict
