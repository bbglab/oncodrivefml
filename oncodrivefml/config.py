import logging
import os

from shutil import copyfile

import sys
from validate import Validator
from os.path import expanduser, expandvars, join, sep, exists
from configobj import ConfigObj, interpolation_engines, flatten_errors
from bgdata.configobj import BgDataInterpolation
interpolation_engines['bgdata'] = BgDataInterpolation

CONFIG_SPECS = {
    'signature': {
        'method': "option('full', 'complement', 'bysample', 'file')",
        'path': "string(default=None)",
        'column_ref': "string(default=None)",
        'column_alt': "string(default=None)",
        'column_probability': "string(default=None)"
    },
    'score': {
        'file': 'string', 'chr': 'integer', 'chr_prefix': 'string', 'pos': 'integer', 'ref': 'integer',
        'alt': 'integer', 'score': 'integer', 'element': 'integer(default=None)', 'extra': 'integer(default=None)'
    },
    'background': {
        'sampling': 'integer', 'recurrence': 'boolean', 'range': 'integer(default=None)'
    },
    'statistic': {
        'method': "option('amean', 'gmean', 'maxmean')"
    },
    'settings': {
        'cores': 'integer(default=None)', 'drmaa': 'integer(default=None)', 'drmaa_maximum': 'integer(default=100)', 'queues': 'list(default=None)'
    }
}


def file_exists_or_die(path):
    path = expandvars(expanduser(path))
    if path is not None:
        if not exists(path):
            logging.error("File '{}' not found".format(path))
            sys.exit(-1)
    return path


def file_none_exists_or_die(path):
    if path is None:
        return None
    return file_exists_or_die(path)


def file_name(file):
    if file is None:
        return None
    return os.path.basename(file).split('.')[0]


def get_bbglab_home():

    bbglab_home = expandvars(expanduser(os.getenv('BBGLAB_HOME', '~/.bbglab/')))
    bbglab_home = bbglab_home.rstrip(os.path.sep)

    return join(sep, bbglab_home)


def get_oncodrivefml_config_file():

    file_path = 'oncodrivefml.conf'

    # Check if the configuration file exists in current folder
    if exists(file_path):
        return file_path

    # Check if exists in the default configuration folder
    bbglab_home = get_bbglab_home()
    file_path = join(bbglab_home, 'oncodrivefml.conf')
    if exists(file_path):
        return file_path

    # Otherwise, create the configuration file from the template
    if not exists(bbglab_home):
        os.makedirs(bbglab_home)
    config_template = os.path.join(os.path.dirname(__file__), "oncodrivefml.conf.template")
    copyfile(config_template, file_path)

    return file_path


def load_configuration(config_file):

    if config_file is None:
        config_file = get_oncodrivefml_config_file()

    config = ConfigObj(config_file, configspec=CONFIG_SPECS, interpolation='bgdata')
    res = config.validate(Validator(), preserve_errors=True)

    for section_list, key, error in flatten_errors(config, res):

        if key is not None:
            section_list.append(key)
        else:
            section_list.append('[missing section]')
        section_string = ' > '.join(section_list)

        if not error:
            error = 'Missing value or section.'

        logging.error("Config error at {} = {}".format(section_string, error))

    if res != True:
        sys.exit(-1)

    return config


