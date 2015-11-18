import sys
from setuptools import setup, find_packages
from oncodrivefml import __version__

if sys.hexversion < 0x03000000:
    raise RuntimeError('This package requires Python 3.0 or later.')

setup(
    name='oncodrivefml',
    version=__version__,
    packages=find_packages(),
    package_data={'oncodrivefml': ['*.txt.gz']},
    url="https://bitbucket.org/bbglab/oncodrivefml",
    download_url="https://bitbucket.org/bbglab/oncodrivefml/get/"+__version__+".tar.gz",
    license='UPF Free Source Code',
    author='Biomedical Genomics Group',
    author_email='nuria.lopez@upf.edu',
    description='',
    install_requires=[
        'configobj >= 5.0.6',
        'numpy >= 1.9.0',
        'scipy >= 0.14.0',
        'statsmodels >= 0.6.1',
        'pytabix >= 0.1',
        'bokeh >= 0.9.2',
        'pandas >= 0.15.2',
        'matplotlib >= 1.4.0',
        'intervaltree >= 2.1.0',
        'bgqmap >= 1.0.0',
        'bgdata >= 0.6.0',
        'itab >= 0.3.0'
    ],

    entry_points={
        'console_scripts': [
            'oncodrivefml = oncodrivefml.main:cmdline'
        ]
    }
)
