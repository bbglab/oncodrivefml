import logging
from collections import defaultdict
from datetime import datetime

from ago import human

from oncodrivefml.executors.bymutation import ElementExecutor


def executor_run(executor: ElementExecutor):
    return executor.run()


def defaultdict_list():
    return defaultdict(list)


def loop_logging(iterable, size=None, step=1):

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
