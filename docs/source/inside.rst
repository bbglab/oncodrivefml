
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

.. _inside cli type:

The :ref:`following table <cli type>` shows the
modifications introduced
in the :ref:`indels configuration parameters <config indels>`
by the :command:`--type` flag:


.. table:: Effects of --type
   :name: cli type

   ======================  ========================================
   Value                   Effect in configurtion of indels
   ======================  ========================================
   coding                  ``method = 'stop'``
   noncoding               ``method = 'max'``
   ======================  ========================================

.. _inside cli noindels:

The flag :command:`--no-indels` also affects the
:ref:`indels configuration parameters <config indels>`.
Particularly, it has effect on the ``include`` option.
The use of this flag discards the analysis of indels
by setting ``include = False``, while not using it
includes indels (``include = True``).

.. _inside cli sequencing:

The :ref:`table below <cli sequencing>` shows the effects of the
:command:`--sequencing` flag in the :ref:`signature configuration <config signature>`:


.. table:: Effects of --sequencing
   :name: cli sequencing

   ======================  ========================================
   Value                   Effect in signature
   ======================  ========================================
   wgs                     ``normalize_by_sites = 'whole_genome'``
   wes                     ``normalize_by_sites = 'whole_exome'``
   targeted                ``normalize_by_sites = None``
   ======================  ========================================

.. _inside cli signature:

Using the :command:`--signature` of the command line,
set the signature configuration to
``method = "file"`` and ``path = "<provided path>"``

.. note:: Signatures provided as an external file are not normalized.

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

   bgdata datasets genomereference hg19


If you want to use a different genome, you need to
modify the code in the :mod:`oncodrivefml.signature` module.


Gene stops
^^^^^^^^^^

OncodriveFML also uses a tabix file that contains the
positions and the alterations of the gene stops.


.. code-block:: bash

   bgdata datasets genestops hg19
