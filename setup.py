
from os import path

from setuptools import Extension
from setuptools import setup, find_packages
from oncodrivefml import __version__


directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()

setup(
    name='oncodrivefml',
    version=__version__,
    packages=find_packages(),
    package_data={'oncodrivefml': ['*.txt.gz', '*.conf.template', '*.conf.template.spec', '*.pyx', '*.json.gz']},
    url="https://bitbucket.org/bbglab/oncodrivefml",
    download_url="https://bitbucket.org/bbglab/oncodrivefml/get/"+__version__+".tar.gz",
    license='UPF Free Source Code',
    author='BBGLab (Barcelona Biomedical Genomics Lab)',
    author_email='bbglab@irbbarcelona.org',
    description='Identify signals of positive selection in somatic mutations',
    setup_requires=[
        'cython',
    ],
    install_requires=required,
    ext_modules=[Extension('oncodrivefml.walker_cython', ['oncodrivefml/walker_cython.pyx'])],
    entry_points={
        'console_scripts': [
            'oncodrivefml = oncodrivefml.main:cmdline'
        ]
    }
)
