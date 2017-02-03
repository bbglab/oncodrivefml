############
OncodriveFML
############

Recent years saw the development of methods to detect signals of positive selection in the pattern of somatic mutations in genes across cohorts of tumors, and the discovery of hundreds of driver genes. The next major challenge in tumor genomics is the identification of non-coding regions which may also drive tumorigenesis. We present OncodriveFML, a method that estimates the accumulated functional impact bias of somatic mutations in any genomic region of interest based on a local simulation of the mutational process affecting it. It may be applied to all genomic elements to detect likely drivers amongst them. OncodriveFML can discover signals of positive selection when only a small fraction of the genome, like a panel of genes, has been sequenced.

License
=======
OncodriveFML is made available to the general public subject to certain conditions described in its `LICENSE <LICENSE>`_. For the avoidance of doubt, you may use the software and any data accessed through UPF software for academic, non-commercial and personal use only, and you may not copy, distribute, transmit, duplicate, reduce or alter in any way for commercial purposes, or for the purpose of redistribution, without a license from the Universitat Pompeu Fabra (UPF). Requests for information regarding a license for commercial use or redistribution of OncodriveFML may be sent via e-mail to innovacio@upf.edu.

Installation
============

OncodriveFML depends on Python 3.4 and some external libraries. The easiest way to install all this software stack is using the well known `Anaconda Python distribution <http://continuum.io/downloads>`_.

Then to get OncodriveFML installed first clone the repository and then install it using ``pip``::

        $ git clone git@bitbucket.org:bbglab/oncodrivefml.git
        $ cd oncodrivefml
        $ pip install .

The first time that you run OncodriveFML it will download the genome reference from our servers. By default the downloaded datasets go to ``~/.bgdata`` if you want to move this datasets to another folder you have to define the system environment variable BGDATA_LOCAL with an export command. 

The following command will show you the command help::

	$ oncodrivefml --help 
    usage: oncodrivefml [-h] -i MUTATIONS_FILE -e ELEMENTS_FILE [-o OUTPUT_FOLDER]
                    [-c CONFIG_FILE] [--samples-blacklist SAMPLES_BLACKLIST]
                    [--debug]

    optional arguments:
    -h, --help            show this help message and exit
    -i MUTATIONS_FILE, --input MUTATIONS_FILE
                        Variants file
    -e ELEMENTS_FILE, --elements ELEMENTS_FILE
                        Genomic elements to analyse
    -o OUTPUT_FOLDER, --output OUTPUT_FOLDER
                        Output folder. Default to regions file name without
                        extensions.
    -c CONFIG_FILE, --configuration CONFIG_FILE
                        Configuration file. Default to 'oncodrivefml.conf' in
                        the current folder if exists or to
                        ~/.bbglab/oncodrivefml.conf if not.
    --samples-blacklist SAMPLES_BLACKLIST
                        Remove this samples when loading the input file
    --debug             Show more progress details

      
File formats
============

.. note::

   All the files can be compressed using GZIP (extension ".gz"), BZIP2 (ext. ".bz2") or LZMA (ext. ".xz")

Input file format
-----------------

The variants file is a text file with 5 columns separated by a tab character (without any header or comment):

* Column 1: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column 2: Mutation position. A positive integer.
* Column 3: Reference allelle. A single letter: A, C, G or T (upper case)
* Column 4: Alternate allelle. A single letter: A, C, G or T (upper case)
* Column 5: Sample identifier. Any alphanumeric string.
      
Regions file format
--------------------

The regions file is a text file with 4 columns separated by a tab character (without any header or comment):

* Column 1: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column 2: Start position. A positive integer.
* Column 3: End position. A positive integer.
* Column 4: Element identifier.
