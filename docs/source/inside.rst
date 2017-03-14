
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

   This overwrite is perform regardless the parameter is set or not in the configuration file.


In the :ref:`following table <cli indels>`, the effects of the :command:`--indels` flag
on the :ref:`indels configuration parameters <config indels>`
are shown:



.. table:: Effects of --indels
   :name: cli indels

   ======================  ========================================
   Value                   Effect in indels
   ======================  ========================================
   discard                 ``include = False``
   coding                  ``include = True``
                           ``method = 'stop'``
   noncoding               ``include = True``
                           ``method = 'max'``
   ======================  ========================================


The :ref:`table below <cli sequencing>` shows the effects of the
:command:`--sequencing` flag in the :ref:`signature configuration <config signature>`:


.. table:: Effects of --sequencing
   :name: cli sequencing

   ======================  ========================================
   Value                   Effect in signature
   ======================  ========================================
   genome                  ``normalize_by_sites = 'whole_genome'``
   exome                   ``normalize_by_sites = 'whole_exome'``
   other                   none
   ======================  ========================================

Finally, the use of the :command:`--debug` flag
sets the level of the *console* handler in the :ref:`logging section <config logging>`
to ``'DEBUG'``.

.. _inside pickles:

Pickle files
------------

OncodriveFML can create and use intermediate files
to speed up computations that use the same files.

The regions file is loaded using the BgParsers library,
so the cache of that file is out of the scope of
OncodriveFML. In short, the file will be cached
the first time you use it and rebuild
if you change its name or content.

There are 2 other items for which OncodriveFML
can create or use a cache-like files to speed up future executions.
Those files are saved in (or loaded from) the same folder
as the mutations file.
However, the systems is not as sophisticated as the BgParsers and may
lead to different issues.
To generate these cache-like files
you need to run OncodriveFML with the
:command:`--generate-pickle` option
(see :ref:`command line interface <oncodrive help cmd>`).

.. warning::

   Using this option can speed up computations as some steps
   can be replaced by a single file read. However, changes
   in the input files are not noticed by these pickle files
   unless you rename them.
   Thus we recommend it use only for advanced users understanding
   the process.

Mutations
^^^^^^^^^

One of the pickle files that can be created contains
a dictionary with the mutations mapped to the genomic
elements being analysed and some other useful metadata
(such as number of indels or SNP mutations).
This file is: ``<mutations file>+__mapping__+<elements file>``.
This file is helpful to skip the loading and mapping
mutations step.
If exists next to the mutations file, OncodriveFML loads it
as long as it does not receive any file with the blacklisted samples.

Signature
^^^^^^^^^

The other pickle file created is the
signature pickle.
It is only created for signature methods: ``full`` and ``complement``
It name is: ``<mutations file>+_signature_+<method>+_+<classifier>``.
See :ref:`signature configuration <config signature>` for more details
about the methods, classifiers... for the signature.

If exists next to the mutations file, OncodriveFML loads it
as long as it does not receive any file with the blacklisted samples
and the ``only_mapped_mutations`` option is not used
(see :ref:`signature configuration <config signature>`).

.. _inside bgdata:

BgData
------

OncodriveFML uses external data retrieved using the BgData package.
You can download and check this data yourself. If you want to
use different data, you can download the source code
and replace the code to use your own data.

Reference genome
^^^^^^^^^^^^^^^^

As March 2017 BgData includes three reference genomes: *HG18*, *HG19*
and *HG38*.

.. code-block:: bash

   bgdata datasets genomereference hg19


If you want to use a different genome, you need to
modify the code in the :mod:`oncodrivefml.signature` module.

Signature correction
^^^^^^^^^^^^^^^^^^^^

BgData includes the counts of the triplets
in whole exome and whole genome.

.. code-block:: bash

   bgdata datasets exomesignature hg19

   bgdata datasets genomesignature hg19


Those counts are used to compute trinucleotides
frequencies and perform signature correction
(find more details in the :ref:`signature <signature>` section
and in the :ref:`signature configuration <config signature>`).

Gene stops
^^^^^^^^^^

OncodriveFML also uses a tabix file that contains the
positions and the alteration of the gene stops.


.. code-block:: bash

   bgdata datasets genestops hg19 TODO
