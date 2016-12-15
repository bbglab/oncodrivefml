"""
This module contains some useful methods
"""

import logging
from collections import defaultdict
from datetime import datetime
from ago import human

from oncodrivefml.executors.bymutation import ElementExecutor


def executor_run(executor: ElementExecutor):
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
            logging.info("[{} of {}]".format(i+1, size))
        yield value
    logging.info("[{} of {}]".format(i+1, size))
    logging.debug("Time: {}".format(human(start_time)))
