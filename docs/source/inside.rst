
Behind the scenes
=================

This section will point out some parts which
might be interesting if you are running
OncodriveFML yourself.

.. _inside cli:

Command line interface
----------------------

The command line interface of OncodriveFML overwrites some of the
parameters in the configuration file.

.. warning::

   This overwrite is performed regardless the parameter is set or not in the configuration file.

.. _inside cli noindels:

The flag :command:`--no-indels` also affects the
:ref:`indels configuration parameters <config indels>`.
Particularly, it has effect on the ``include`` option.
The use of this flag discards the analysis of indels
by setting ``include = False``.

.. _inside cli signature:

Using the :command:`--signature` of the command line,
set the signature configuration to
``method = "file"`` and ``path = "<provided path>"``

.. note:: Signatures provided as an external file are not normalized.


.. _inside cli signature correction:

The :ref:`table below <cli signature correction>` shows the effects of the
:command:`--signature-correction` flag in the :ref:`signature configuration <config signature>`:


.. table:: Effects of --signature-correction
   :name: cli signature correction

   ======================  ========================================
   Value                   Effect in signature
   ======================  ========================================
   wg                      ``normalize_by_sites = 'whole_genome'``
   wx                      ``normalize_by_sites = 'whole_exome'``
   ======================  ========================================

.. note:: This option does not have any impact if signatures
   are passed with the ``--signature`` option.

.. _inside bgdata:

BgData
------

OncodriveFML uses external data retrieved using the `BgData package <https://bitbucket.org/bgframework/bgdata>`_.
You can download and check this data yourself. If you want to
use different data, you can download the source code
and modify the code to use your own data.

Reference genome
^^^^^^^^^^^^^^^^

As March 2017 BgData includes three reference genomes: *HG18*, *HG19*
and *HG38*.

.. code-block:: bash

   bgdata datasets/genomereference/hg38


If you want to use a different genome, you need to
modify the code in the :mod:`oncodrivefml.signature` module.


Gene stops
^^^^^^^^^^

OncodriveFML also uses a tabix file that contains the
positions and the alterations of the gene stops.


.. code-block:: bash

   bgdata datasets/genestops/hg38
