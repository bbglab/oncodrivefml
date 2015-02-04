from distutils.core import setup
from setuptools import find_packages
from oncodrivefm2 import VERSION

setup(
    name='oncodrivefm2',
    version=VERSION,
    packages=find_packages(),
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
        'pandas >= 0.15.2'
    ],

    dependency_links=[
        'git+https://bitbucket.org/bbglab/bgcore.git@develop#egg=bgcore-0.6.0'
    ],

    entry_points={
        'console_scripts': [
            'oncodrive-fm2 = oncodrivefm2.oncodrivefm2:cmdline'
        ]
    }
)
