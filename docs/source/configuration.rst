.. _project configuration:

Configuration
=============

The method behaviour can be modified through a configuration file.

.. warning::

   Using the command line interface overwrites some setting in the configuration file.
   Check how the command line interface changes the configuration in the :ref:`command line interface <inside cli>`
   section.

Check the :file:`oncodrivefml_v2.conf.template` that is included in the package
to find an example of the configuration file.

This section will explain each of the parameters in the configuration file:

.. _config genome:

Genome
------

.. literalinclude:: oncodrivefml_v2.conf.template
   :language: text
   :lines: 1-4

The genome section makes reference to the reference genome
used by OncodriveFML.

The reference genome has been obtained from `<http://hgdownload.cse.ucsc.edu/downloads.html>`_.

Currently, only ``HG19`` and ``HG38`` are fully supported.
Use :code:`build = 'hg19'` or :code:`build = 'hg38'` to use any of them.

There is a partial support for other genomes.
The support is only partial because the values for the position and alterations
of the stops in the these genomes have not been computed yet. If you want to
run OncodriveFML with any of these genomes, make sure you do not use
the ``stop`` method for the indels (:ref:`ref <config indels>`).

.. warning::

   If you decide to use a reference genome other than ``HG19``, make sure that the
   scores file you use is compatible with it.
   And you update all related parameters.

.. _config signature:

Signature
---------

.. literalinclude:: oncodrivefml_v2.conf.template
   :language: text
   :lines: 6-14

The signature represents the probability of a certain nucleotide
to mutate taking into account its context [#context]_.

The signature can be configured using the following parameters:

method
  Method used to compute the signature. Options are:

  - ``method = 'none'``: all changes
    have equal probability.
    This approach is recommended for small datasets.

  - ``method = 'complement'``:
    use a 96 matrix with the signatures complemented.

  - ``method = 'full'``:
    use a 192 matrix with all the possible signatures.

  - ``method = 'bysample'``:
    equivalent to ``method = 'complement'`` and
    ``classifier = 'SAMPLE'``

  - ``method = 'file'``:
    use a precomputed signature. This option requires
    to add the path to the file as ``path = '/path/to/file'``.
    Note that this option *can* be overwritten by the
    :ref:`--signature option <inside cli signature>`
    in the command line interface.


classifier
  The signature by default is computed using the whole dataset.
  However, you can group the mutations in categories that
  correspond to any of the values in
  ``SAMPLE``, ``CANCER_TYPE`` and ``SIGNATURE``
  columns (if provided).

  If a file with the signature is provided, and
  that signature has also been computed using groups,
  the same classifier must be specified.

normalize_by_sites
  Compute a normalization of the signature.
  This option appears commented because it *is* overridden
  by the :ref:`--sequencing option <inside cli sequencing>`
  of the command line interface.

  If you provide an external file with the signature,
  it will never be corrected, regardless of the value
  of this option.



The recommended approach is to compute your own signature
(e.g. using the `bgsignature package <https://bitbucket.org/bgframework/bgsignature>`_)
and pass it to OncodriveFML.


.. _config score:

Score
-----

The score section is used to know
which scores are going to be used.

.. literalinclude:: oncodrivefml_v2.conf.template
   :language: text
   :lines: 17-40

The scores should be a file that for a given position, in a given chromosome, 
gives a value to every possible alteration.

Some of the parameters in this section are optional,
while others are mandatory.

file
  It is a string and represents the path to the scores file.

format
  Indicates the format of the file. Options are:

  - ``format = 'tabix'`` indicates that the file is a
    tab separated file compressed with bgzip.
    This means that a .tbi index file should be present in the same location.
  - ``format = 'pack'`` is a
    binary format we have implemented to reduce the file size.
    It is only available for specific scores.
    Thus, if you want to use your own file, use the `tabix <http://www.htslib.org>`_ format.

chr
  Column in the file where the chromosome is indicated.
chr_prefix
  When querying the tabix file for a specif chromosome
  OncodriveFML only uses the number of the chromosome or 'X' or 'Y'. If the
  tabix file requires a prefix before the chromosome, use this option. For instance, if the chromosomes in the
  tabix file are labeled as ``chr1``, ``chr2``, .., ``chrY``, set this option to: ``chr_prefix = 'chr'``. 
  If this is not the case, use an empty string: ``chr_prefix = ''``.
pos
  Column that indicates the position of the scored alteration in the chromosome.
ref
  Column that contains the reference allele. It is optional.
alt
  Column that contains the alternate allele. It is optional.
  If is not specified, it is assumed that the 3 possible changes have the same score.
score
 Column that contains the score.
element
 Column that contains the element identifier. It is optional.
 If it is provided and the value does not match with the one from the regions,
 these scores are discarded.

OncodriveFML uses two additional parameters
(``mean_to_stop_function`` and ``minimum_number_of_stops``),
which are related only to the ``stop`` method
for :ref:`computing the indels <analysis indel>`.

mean_to_stop_function
  When analysing a certain gene, OncodriveFML might need
  to score an indel according to the value of the stops in the gene.
  It might happen that the number of stops is 0
  or is below a certain threshold.
  In such cases, OncodriveFML uses the function specified
  in ``mean_to_stop_function`` parameter to assign a score from the mean value
  of all the stops in the gene.

  .. only:: html

     Download the :download:`IPython notebook <_static/ScoresFunction.html>` that
     has been created with the functions computed for
     *CADD1.0* and *CADD1.3*, or :ref:`see it <nb scoresfunct>`.

minimum_number_of_stops
  When analysing a certain gene, OncodriveFML gets all the scores
  associated with the mutations that produce a stop in that gene.
  ``minimum_number_of_stops`` indicates the minimum number of stops
  that a gene is required to have in order to avoid using the function above.

.. attention:: These parameters must also be adjusted for
   each scores file.

.. _config statistic:

Statistic
---------

The statistic section is related to the configuration
of the analysis

.. literalinclude:: oncodrivefml_v2.conf.template
   :language: text
   :lines: 43-73

There a different parameters you can configure:

method
  Represents the type of operation that is applied to
  observed and simulated scores before comparing them.
  Options are:

  - ``method = 'amean'``: arithmetic mean
  - ``method = 'gmean'``: geometric mean


discard_mnp
  Indicates whether to include or not MNP mutations in the analysis.

  - ``discard_mnp = False``: include them
  - ``discard_mnp = True``: discard them


per_sample_analysis
  In some cases, you might be interested in performing the
  analysis per sample. This means that all the mutations that come
  from the same sample are reduced to a single score.
  This score can be computed as:

  - ``per_sample_analysis = 'max'``: maximum score
  - ``per_sample_analysis = 'amean'``: arithmetic mean
  - ``per_sample_analysis = 'gmean'``: geometric mean

  Comment this option if you are not interested in this type of analysis.


OncodriveFML includes a few more parameters
that are related to how many simulations are performed.

sampling
  Represents the minimum number of
  simulations to be performed.

sampling_max
  Represents the maximum number
  of simulations to be performed.

sampling_chunk
  Represents the maximum size (in millions)
  that a single process can handle. This value is
  used to keep the memory usage within certain limits.

  .. note::

     With a value of 100, each process takes less than 4 GB
     of RAM. We have not considered the memory taken by
     the main process.

sampling_min_obs
  Represents the minimum
  number of observations [#obs]_. When it is reached,
  no more simulations are performed.

.. _config indels:

Indels
^^^^^^

The indels subsection of statistic contains
the configuration for the analysis of indels.

.. literalinclude:: oncodrivefml_v2.conf.template
   :language: text
   :lines: 63-73

OncodriveFML accepts various parameters related to the indels:


include
  Indicates
  whether to include indels in the analysis
  or not.

  - ``include = True``: include indels in the analysis
  - ``include = False``: exclude indels from the analysis

  This option *is* overridden by the
  :ref:`--no-indels flag <inside cli noindels>`
  of the command line interface.

method
  Indicates how to simulate the indels.

  - ``method = 'max'``: simulates the
    indels as a set of
    substitutions.
    Indels that are simulated as substitutions [#indelsSubs]_
    follow the same signature patter as the mutatinal signature.
  - ``method = 'stop'``: simulates
    indels as stops. This option is recommended
    for simulating indels in coding regions.

  Check the :ref:`analysis of indels <analysis indel>`
  section to find more details.

  This option *is* overridden by
  the :ref:`--type opton <inside cli type>`
  in the command line interface.

max_consecutive
  OncodriveFML discards indels that fall in
  repetitive regions. OncodriveFML considers that
  an indel is in a repetitive region when the
  same sequence of the indel appears consecutively 
  in a genomic element a certain number of times 
  (or even more).
  The maximum number of consecutve repetitions can be 
  set with the ``max_consecutive`` option.
  OncodriveFML will not discard any indel
  due to repetitive regions if you set
  ``max_consecutive = 0``.


.. _exomic frameshift rate:

gene_exomic_frameshift_ratio
  Indicates which mutations influence the :ref:`probabilities <analysis probs>`
  for frameshift indels and substitutions.

  - ``gene_exomic_frameshift_ratio = False``: the probabilities are taken
    from the mapped mutations discarding those whose length is
    multiple of 3. Note that in order to work properly,
    this option should be set only when the regions file corresponds to
    coding regions.
  - ``gene_exomic_frameshift_ratio = True``:
    probabilities are taken from the observed mutations rate in each region.

  This option is harmless when ``method = 'max'``,
  as all indels are simulated as substitutions.

stops_function
  The *observed* score of an indel that is computed with the
  ``method = 'stop'`` option is related to the score of the stops
  in its gene. You can decide how this relation is by choosing
  a function that is applied to all stops scores in the gene.

  - ``stops_function = 'mean'``: associates the indel to a value that is
    equal to the *mean* of all stop scores in the gene
  - ``stops_function = 'median'``: associates the indel to a value that is
    equal to the *median* of all stop scores in the gene
  - ``stops_function = 'random'``: associates the indel to a value that is
    a *random value between the maximum and the minimum* of all stop scores in the gene
  - ``stops_function = 'random_choice'``: associates the indel to a value that is
    a *random value* between all the possible stop scores in the gene

.. _config settings:

Settings
--------

To configure the system where the analysis is performed
OncodriveFML includes the setting section:

.. literalinclude:: oncodrivefml_v2.conf.template
   :language: text
   :lines: 76-78

Use the ``cores`` option to indicate how many cores to
use. You can comment this option in order to use
all the available cores.

The command line :command:`--cores` option
*can* override this value.

.. note::

   OncodriveFML works on shared memory systems
   using the :mod:`multiprocessing` module.

----

.. [#context] Previous and posterior nucleotides

.. [#obs] An observation is counted when a simulated value,
   after applying the function in ``method`` to the simulated scores, 
   is higher than the result of applying the same function to the 
   observed scores.

.. [#indelsSubs] All indels are simulated as substitutions when
   ``method = 'max'``. Indels that are in-frame
   are also simulated as substitutions when ``method = 'stop'``.
