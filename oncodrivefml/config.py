import logging
import os

import sys
from bgconfig import BGConfig, _file_name, _file_exists_or_die

file_exists_or_die = _file_exists_or_die
file_name = _file_name


def load_configuration(config_file):
    """
    Load the configuration file and check format.

    :param config_file: The path to the configuration file. If it's none, use default location.
    :return: The configuration dictionary
    """

    config_template = os.path.join(os.path.dirname(__file__), "oncodrivefml.conf.template")

    try:
        return BGConfig(config_template, config_file=config_file)
    except ValueError as e:
        logging.error(e)
        sys.exit(-1)



