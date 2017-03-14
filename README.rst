OncodriveFML
============

Recent years saw the development of methods to detect signals of positive selection in the pattern of somatic mutations in genes across cohorts of tumors, and the discovery of hundreds of driver genes. The next major challenge in tumor genomics is the identification of non-coding regions which may also drive tumorigenesis. We present OncodriveFML, a method that estimates the accumulated functional impact bias of somatic mutations in any genomic region of interest based on a local simulation of the mutational process affecting it. It may be applied to all genomic elements to detect likely drivers amongst them. OncodriveFML can discover signals of positive selection when only a small fraction of the genome, like a panel of genes, has been sequenced.

License
-------
OncodriveFML is made available to the general public subject to certain conditions described in its `LICENSE <LICENSE>`_. For the avoidance of doubt, you may use the software and any data accessed through UPF software for academic, non-commercial and personal use only, and you may not copy, distribute, transmit, duplicate, reduce or alter in any way for commercial purposes, or for the purpose of redistribution, without a license from the Universitat Pompeu Fabra (UPF). Requests for information regarding a license for commercial use or redistribution of OncodriveFML may be sent via e-mail to innovacio@upf.edu.

Installation
------------

OncodriveFML depends on Python 3.5 and some external libraries.
The easiest way to install all this software stack is using the well known `Anaconda Python distribution <http://continuum.io/downloads>`_::

    $ conda install -c bbglab oncodrivefml

OncodriveFML can also be installed using ``pip``::

    pip install oncodrivefml

Finally, you can get the latest code from the repository and install with ``pip``::

        $ git clone git@bitbucket.org:bbglab/oncodrivefml.git
        $ cd oncodrivefml
        $ pip install .

The first time that you run OncodriveFML it will download the genome reference from our servers.
By default the downloaded datasets go to ``~/.bgdata`` if you want to move this datasets to another folder you have to define the system environment variable BGDATA_LOCAL with an export command.

The following command will show you the command help::

	$ oncodrivefml --help 
