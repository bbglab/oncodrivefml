from collections import defaultdict

from oncodrivefml.executors import ElementExecutor


def executor_run(executor: ElementExecutor):
    return executor.run()


def defaultdict_list():
    return defaultdict(list)