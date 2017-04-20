
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


The flag :command:`--no-indels` also affects the
:ref:`indels configuration parameters <config indels>`.
Particularly, it has effect on the ``include`` option.
The use of this flag discards the analysis of indels
by setting ``include = False``, while not using it
includes indels (``include = True``).

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
lead to few issues.
To generate these cache-like files
you need to run OncodriveFML with the
:command:`--generate-pickle` option
(you can list all the options :ref:`using the help <help cmd>`).

.. warning::

   Using this option can speed up computations as some steps
   can be replaced by a single file read. However, changes
   in the input files are not noticed by these pickle files
   unless you rename them.
   Thus we recommend its use only to advanced users that understand
   the process.

Mutations
^^^^^^^^^

One of the pickle files that can be created contains
a dictionary with the mutations mapped to the genomic
elements being analysed and some other useful metadata
(such as the number of indels or SNP mutations).
This file, named ``<mutations file>+__mapping__+<elements file>``,
is helpful to skip the steps of loading and mapping
mutations.
If this file is in the same location as the mutations file, OncodriveFML loads it
as long as it does not receive any file with blacklisted samples.

Signature
^^^^^^^^^

The other pickle file created is the
signature pickle.
It is only created for signature methods: ``full`` and ``complement``
Its name is: ``<mutations file>+_signature_+<method>+_+<classifier>``.
See :ref:`signature configuration <config signature>` for more details
(methods, classifiers, etc.) about the signature.

If this file is located in the same directory as the mutations file, OncodriveFML loads it
as long as it does not receive any file with  blacklisted samples
and the ``only_mapped_mutations`` option is not used
(see :ref:`signature configuration <config signature>`).

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

Signature correction
^^^^^^^^^^^^^^^^^^^^

BgData includes the counts of the triplets
in whole exome and whole genome.

.. code-block:: bash

   bgdata datasets exomesignature hg19

   bgdata datasets genomesignature hg19


Those counts are used to compute the trinucleotides
frequencies and to perform signature correction
(find more details in the :ref:`signature <signature>` section
and in the :ref:`signature configuration <config signature>`).

Gene stops
^^^^^^^^^^

OncodriveFML also uses a tabix file that contains the
positions and the alterations of the gene stops.


.. code-block:: bash

   bgdata datasets genestops hg19
