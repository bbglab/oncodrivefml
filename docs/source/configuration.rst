.. _project configuration:

Configuration
=============

The method behaviour can be modified through a configuration file.

.. warning::

   Using the command line interface overwrites some setting in the configuration file.
   Check how the command line interface changes the configuration in the :ref:`command line interface <inside cli>`
   section.

Check the :file:`oncodrivefml.conf.template` that is included in the package
to find an example of the configuration file.

This section will explain each of the parameters in the configuration file:

.. _config genome:

Genome
------

.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :lines: 1-4

The genome section makes reference to the reference genome
used by OncodriveFML.

The reference genome has been obtained from `<http://hgdownload.cse.ucsc.edu/downloads.html>`_.

Currently, only ``HG19`` is fully supported. Use :code:`build = 'hg19'` to use it.

There is partial support for ``HG18`` and ``HG38``.
The support is only partial because the values for the position and alterations
of the stops in the these genomes has not been computed yet. If you want to
run OncodriveFML with any of these genomes, make sure you do not use
the ``stop`` method for the indels(:ref:`ref <config indels>`).

.. warning::

   If you decide to use a refernce genome other than ``HG19``, make sure that the
   scores file you use in compatible with it.

.. _config signature:

Signature
---------

.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :lines: 8,12-13,32-33,32-38,40-42,45-48,51-52,60-61

The signature represents the probability of certain nucleotide
to mutate taking into account its context [#context]_.

You can choose one of the following options for the signature:

- Not use any signature, which is equivalent to assume that all changes
  have equal probability: ``method = 'none'``.
  This approach is recommended when the dataset is small.

- OncodriveFML can also compute a signature using the provided dataset.
  This option contains a set of parameters that you can use to decide how
  this computation is done.

  - Select one of the methods to compute the signature from the dataset:
    ``method = 'full'`` to count each mutation once and
    ``method = 'complement'`` to collapse complementary mutations.

    .. note::

       The option ``method = 'bysample'`` is equivalent to ``method = 'complement'``
       but forces the classifier (see below) to be ``SAMPLE``.

  - The classifier parameter indicates which column from the mutations file
    is used to group the mutations when computing the signature. E.g. grouping
    by ``SAMPLE`` generates one signature for each sample.
    Only ``SAMPLE``, ``CANCER_TYPE`` and ``SIGNATURE`` columns can be used.

  - You can decide to use only SNP (``include_mnp = False``)
    or also use MNP mutations (``include_mnp = True``).

  - You can choose between using only the mutation that are mapped
    to the regions under analysis (``only_mapped_mutations = True``)
    or use all the mutations (SNPs and optionally MNPs) in the dataset.

  - The signature can be corrected by the frequency of sites.
    If you do not specify anything, OncodriveFML will not correct the signature.
    Use ``normalize_by_sites = 'whole_genome'`` or ``'wgs'`` to correct
    by the frequencies in the whole genome.
    Use ``normalize_by_sites = 'whole_exome'`` or ``'wes'`` or ``'wxs'`` to correct
    by the frequencies in the exome.
    If you have specified ``only_mapped_mutations = True``, then the correction will
    be done by the frequencies of trinuceotides found in the regions under analysis,
    as long as you indicate one of the above mentioned values.

    .. note::

       The frequencies have been computed for genome build ``HG19``.
       If you want to check the values, use the :ref:`bgdata package <inside bgdata>`.

- The recommended approach is to use your own signatures.
  OncodriveFML has the option ``method = 'file'`` to
  load precomputed signatures from a file.
  This option requires a few additional parameters:

  - ``path``: path to the file containing the signature
  - ``colum_ref``: column that contains the reference triplet
  - ``column_alt``: column that contains the alternate triplet
  - ``column_probability``: column that contains the probability

    .. warning::

       Probabilities must sum to one.

.. _config score:

Score
-----

The score section is used to know
which scores are going to be used.

.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :lines: 65-67,75,77,79-80,82-83,85-86,88-89,91-92,94-95,101-103,105-106

The scores should be a file that by for a particular
SNP mutation in a given chromosome, at a given position,
for each alteration gives a value to that mutation.

Some of the parameters in this section are optional,
where other are mandatory.

- ``file`` is a string with the path to the scores file.
- ``format = 'tabix'`` indicates that the file is a
  tab separated file compress with bgzip and next to a Tabix .tbi index file
- ``chr`` column in the file where the chromosome is.
- ``chr_prefix``: when asking the tabix file for the specif chromosome
  OncodriveFML only passes the number of the chromosome or 'X' or 'Y'. If the
  tabix file requires anything in front, use this option. If is not the case,
  use the empty string: ``chr_prefix = ''``.
- ``pos`` column that indicates the position
- ``ref`` column that contains the reference allele. It is optional.
- ``alt`` column that contains the alternate allele. It is optional.
  If is not specified, it is assumed that the 3 possible changes have the same score.
- ``score`` column that contains the score.
- ``element`` column that contains the element identifier. It is optional.
  If it is provided and the value does not match with the one from the regions
  those scores are discarded.

OncodriveFML uses two additional parameters,
which are related only to the ``stop`` method
for :ref:`computing the indels <analysis indel>`.

- When analysing a certain gene, OncodriveFML might need
  to score an indel according to the value of the stops in the gene.
  It might happen that the number of stops is 0
  or is below a certain threshold.
  In such cases, OncodriveFML uses the function specified
  in this parameter to assign a value from the mean value
  of all the stops in the gene.

- When analysing a certain gene, OncodriveFML get all the scores
  associated with the mutations that produce a stop in that gene.
  ``minimum_number_of_stops`` indicates the minimum number of stops
  that a gene must have in order to avoid using the function above.

.. _config statistic:

Statistic
---------

The statistic section is related to the configuration
of the analysis

.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :lines: 110,112,114,119-121,130-132,134-135,137-138,140-141

There a different parameters you can configure:

- ``method`` represents the type of operation that is applied to
  observed and simulated scores before comparing them.
  The arithmetic mean (``method = 'amean'``) and
  the geometric mean (``method = 'gmean'``) are supported.
  The recommended one is the arithmetic mean.

- In some cases, you might be interested in performing they
  analysis per sample. It means that when analysing one region
  from all the mutations that come from the same sample,
  OncodriveFML keeps only 1. That one has a score
  which can be the maximum score of all the mutation that
  come from the sample sample (``per_sample_analysis = 'max'``),
  the arithmetic mean (``per_sample_analysis = 'amean'``) or
  the geometric mean (``per_sample_analysis = 'gmean'``).
  Comment this option if you do not want this type of analysis.

- MNP mutations can be included in the analysis optionally.
  ``discard_mnp = False`` to include them
  and ``discard_mnp = True`` to discard them.

OncodriveFML includes a few more parameters,
that are related to how many simulations are performed.
OnocdriveFML tries to keep memory usage within
certain limits.

- ``sampling`` represents the minimum number of
  simulations to be performed.

- ``sampling_max`` represents the maximum number
  of simulations to be performed.

- ``sampling_chunk`` represents the maximum size (in millions)
  that a single process can handle. This value is
  used to keep memory usage whitin certain limits.

- ``sampling_min_obs`` represents the minimum
  number of observations [#obs]_. When it is reached,
  no more simulations are performed.

.. _config indels:

Indels
^^^^^^

Thd indels subsection of statistic contains
the configuration for analysis indels.

.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :lines: 144-146,149-150,155-156,158-161,166-168,171-173,177-178,180-181

OncodriveFML accepts various parameters related to the indels

- The main option is ``include``, which indicates
  OncodriveFML whether to include indels in the anlysis
  or not.
  Use ``inlcude = True`` to include indels
  and ``inlcude = False`` to exclude them.

- OncodriveFML can simulate indels in two ways.
  ``method = 'max'`` simulates indels a set of
  substitutions. ``method = 'stop'`` simulates
  the indels as stop. This option is recommended
  to simulte indels in coding regions.
  Check the :ref:`analysis of indels <analysis indel>`
  section to find more details.

- OncodriveFML discards indels that fall in
  repetitive regions. OncodriveFML considers that
  an indel is in a repetitive region when the
  same sequence of the indel appears, at least,
  certain number of times consecutive in the gene
  following the direction of the strand.
  That number of times is set with the
  ``max_consecutive`` option.
  OncodriveFML will not discard any indel
  due to repetitive regions if you set
  ``max_consecutive = 0``.

- Indels that are simulated as substitutions [#indelsSubs]_
  can be simulated with equal probabilities in all the positions
  of the gene, or using the signature associated with them.
  E.g. if the signature is represented by the cancer type,
  and indel coming from a breast cancer will be simulated
  with the signature of that cancer.
  OncodriveFML has this option because as indels do not contribute to
  the signature, you might decide that they should be be
  simulated with none signature.
  Use ``simulate_with_signature = True`` to use the signature
  or ``simulate_with_signature = False`` to simulate
  equiprobable indels.

.. _exomic frameshift rate:

- ``gene_exomic_frameshift_ratio`` is a flag that indicates OncodriveFML
  which mutations influences the :ref:`probabilities <analysis probs>`
  for frameshift indels and substitutions.
  When ``gene_exomic_frameshift_ratio = False`` the probabilities are taken
  from the mapped mutations discarding those whose length is
  multiple of 3. Note that in order to work properly,
  this option should be set when the regions file corresponds to
  coding regions.
  If ``gene_exomic_frameshift_ratio = True``, the probabilities
  are taken from the observed mutations rate in each region.
  This option is harmless when ``method = 'max'``.

- The *observed* score of an indel that is computed with the
  ``method = 'stop'`` option is related to the score of the stops
  in its gene. You can decide how this relation is by choosing
  a function that is applied to all stops scores in the gene. E.g.
  ``stops_function = 'mean'`` associates the indel a value wich is
  equal to the mean of all stop scores in the gene.
  The options you can choose are:
  - ``'mean'`` for arithmetic mean
  - ``'median'`` for the median
  - ``'random'`` for a random value between the maximum and the minimum
  - ``'random_choice'`` for choosing a random value of all the possible ones


.. _config settings:

Settings
--------

To configure the system where the analysis is performed
OncodriveFML includes the setting section:

.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :lines: 194-197

Use the ``cores`` option to indicate how many cores to
use. You can comment this option in order to use
all available cores.

.. note::

   OncodriveFML works on shared memory systems
   using the :mod:`multiprocessing` module.

.. _config logging:

Logging
-------

The logging section is used to configure the
logging system of OncodriveFML.


.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :lines: 201-230

OncodriveFML uses the :mod:`logging` module.
In particular it loads the configuration file into a dictionary
and passes this section to :func:`~logging.config.dictConfig`.

You can change this section to other compatible configurations
to fit your needs.

All the logs are done using a logger named ``oncodrivefml``.
The logging system can be configured through the logging section of the
:ref:`configuration file <project configuration>`.

.. warning::

   If OncodriveFML detects that the run has already been calculated,
   the warning informing the user uses the root logger.





.. [#context] Previous and posterior nucleotides

.. [#obs] An observation is counted when the simulated scores,
   after applying the function in ``method``, is higher
   the result of applying that same function to the observed
   scores.

.. [#indelsSubs] All indels are simulated as substitutions when
   ``method = 'max'``. Indels that are in-frame
   are also simulated as substitutions when ``method = 'stop'``.
