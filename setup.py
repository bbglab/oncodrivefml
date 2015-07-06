from distutils.core import setup
from setuptools import find_packages
from oncodrivefml import __version__

setup(
    name='oncodrivefml',
    version=__version__,
    packages=find_packages(),
    package_data={'oncodrivefml': ['*.txt.gz']},
    url='http://bg.upf.edu',
    license='UPF Free Source Code',
    author='Biomedical Genomics Group',
    author_email='nuria.lopez@upf.edu',
    description='',
    install_requires=[
        'configobj >= 5.0.6',
        'bgcore >= 0.6.0',
        'numpy >= 1.9.0',
        'scipy >= 0.14.0',
        'pandas >= 0.15.2', 'matplotlib', 'intervaltree'
    ],

    dependency_links=[
        'git+https://bitbucket.org/bbglab/bgcore.git@develop#egg=bgcore-0.6.0'
    ],

    entry_points={
        'console_scripts': [
            'oncodrivefml = oncodrivefml.main:cmdline'
        ]
    }
)
