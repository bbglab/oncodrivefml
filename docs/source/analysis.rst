
.. _analysis:

Analysis
========

This sections explains how OncodriveFML
compute the scores for the observed mutations
and how mutations are simulated.

The analysis is done for each element
independently.
The same number of observed mutations
is simulated within the element,
taking only the positions
indicated in the regions file.


Observed
--------

Single Nucleotide Polymorphism (SNP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SNP mutations are the simplest to compute.
To score them, OncodriveFML get the score
for the corresponding alteration in the position
of the mutation.

If there is not a score for that particular change,
the mutation is ignored [#obsIgnored]_.

Multi Nucleotide Polymorphism (MNP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MNP mutations are considered as set of SNPs.
The observed value is the maximum value of all
the changes produced by the MNP.

MNPs are ignored [#obsIgnored]_
when none of the changes it introduces
has a score.

.. _analysis indel:

Insertion or deletion (INDEL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Indels are scored in two different ways:
as substitutions or as stops.

As substitutions
   Indels that fall in non-coding regions or
   in-frame indels in coding regions are considered as
   a set of substitutions.
   Similarly to MNP mutations, the changes produced by
   the indel are computed as a set of SNPs mutation and OncodriveFML
   assigns the indel the maximum score of those changes.
   In an insertion, the reference genome is compared with the
   indel.
   In a deletion, the reference genome is compared with itself
   but shifted a number of position equal to the length of
   the indel.
   Only the changes produced in the length of the indel are considered.

   .. note::

      If none of the changes produced by the indel has
      a score, the indel is ignored [#obsIgnored]_.

.. _analysis indels as stop:

As stops
   Indels can be scored as stops in the analysis of coding regions
   and if their length is not a multiple of 3.
   In coding regions, a frameshift indel might cause,
   somewhere in the gene, a stop.
   This is why OncodriveFML can use this approach.
   The way OncodriveFML
   scores this type of indels is taking all the stop scores [#stopscores]_ in the gene
   under analysis and applying a user defined function to them.
   In some cases, OncodriveFML can infer a value for the scores of the stops using
   the mean score of all mutations in the gene. See the :ref:`configuration of indel <config indels>`
   section for further information.

   .. attention:: This feature is experimental.
      Thus, it is only available for ``hg19`` and ``hg38`` genomes,
      and it needs to be manually set up in the configuration
      using the configuration file.


Indels with a length higher than 20 nucleotides
are ignored [#obsIgnored]_.
This value can be configured in the
:ref:`configuration file <config indels>`.

Simulated
---------

The same number of mutations that are observed
and have a score are simulated.

To perform the simulation two arrays are computed:

- One contains the scores of all possible changes
  to be simulated.

- The other array contains the probabilities of each of those changes.

Using the probability array, a random sampling of the scores array is
done to obtain the simulated scores.

.. _analysis probs:

Probabilities
^^^^^^^^^^^^^

The probability array is computed taking into account different parameters.

If only substitutions are simulated,
either because the analysis excludes indels or
because they are simulated as substitutions,
the probabilities are:

.. math::

   p = p_{subs} * \frac{\sum_s {p_s*f_s}}{n_{substitutions}}

where ``s`` represents each of the signatures found in the gene in the observed mutations,
:math:`p_s` is the probability of a particular mutation to occur given the ``s`` signature,
:math:`n_{substitutions}` is the total number of substitutions,
and :math:`f_s` is the relative frequency of a particular signature ``s`` in the gene.


However, if you are not using any signature (see :ref:`singature configuration <config signature>`):

.. math::

   p = p_{subs} / {n_{substitutions}}


where :math:`{n_{substitutions}}` is the amount of substitutions in the gene.

----

However, if you configure indels to be analysed as *stops*
things are slightly more complex.
Substitution are simulated as explained above,
as well as in frame indels.
However, there is also a chance that a
the score of one stop is selected.

The probability associated to any of the stop scores is:

.. math::

   p = \frac{1}{n_{stops}} * p_{frameshift indel}

where :math:`p_{frameshift indel} + p_{subs} = 1`, and :math:`n_{stops}` is the number of
stop scores for that gene.

:math:`p_{frameshift indel}` represents the probability of simulating a frameshift indel in that gene,
and :math:`p_{subs}` represents the probability of simulating a substitution.


The probability of simulating a frameshift indel, also,
depends on whether you are analysing using the whole cohort percentages
or only the mutations observed in each gene.

- When using :ref:`exomic frameshift probabilities <exomic frameshift rate>`
  OncodriveFML computes how
  many indels you observe, and how many of those fall into the region
  you are analysing (which should be coding). Among the mapped indels
  OncodriveFML distinguishes between frameshift and in-frame indels.
  The ratio of frameshift indels against the total amount of mutations
  is used to compute :math:`p_{frameshift indel}`.

- When using the probabilities taken from the gene:

  .. math::

     p_{frameshift indel} = \frac{n_{observed frameshift indels}}{n_{observed mutations}}

  where :math:`n_{observed frameshift indels}` is the number of observed frameshift indels
  and :math:`n_{observed mutations}` is the number of observed mutations.

----



.. [#obsIgnored] When an observed mutation is ignored
   it means that it cannot be assigned a score, and thus
   it does not contribute to the observed scores and
   in the simulation the number of mutations simulated is
   one less for that region.

.. [#stopscores] The package BgData includes the precomputed
   position and alteration of the stops for the HG19 genome build.
   OncodriveFML makes use of it.
