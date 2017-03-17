
.. _signature:

Signature
=========

The signature is an array that assigns a probability to
a single nucleotide mutation taking into account its context [#context]_.
It represents the chance of certain mutation to occur within a context.

Check the different options for the signature in the
:ref:`configuration file <config signature>`.
In short, you can choose between not using any signature, using your own signature
or computing the signature from the mutations file.
Additionally, signatures can be grouped into different categories
such as the sample.

The signature array is computed by counting, for each Single Nucleotide Polymorphism,
the reference and alternated triplets.

.. note::

   OncodriveFML also uses the MNP mutations to compute the
   signature, by treating them as a set of separate SNPs.
   You can enable or disable this option with the ``include_mnp`` option in the
   :ref:`configuration file <config signature>`.

The counts are then divided by the total number of counts
to generate a frequency of triplets. For a mutation :math:`i`
the frequency is
:math:`f_i = \frac{m_i}{M}` where :math:`M = \sum_j m_j`, and
:math:`m_i` represent that number of times that the mutation
:math:`i` with its context [#context]_ has been observed.

Optionally, the signature can be corrected taking into
account the frequency of trinucleotides in the
reference genome.
OncodriveFML introduces this feature because the
distribution of triplets is not expected to be constant.
When using the command line interface, OncodriveFML
does this correction automatically according to
the value passed in the flag ``--sequencing``
(you can list all the options :ref:`using the help <help cmd>`).


Reasoning behind the correction
-------------------------------


Let's first take the conditional probability of a mutation (with contectx [#context]_)
to occur given the number of those triplets in the region:
:math:`p_i = p(m = i | T_i) = \frac{m_i}{T_i}`.

Then, the normalized frequency of the mutation :math:`i` is:
:math:`\overline{f_i} = \frac{m_i/T_i}{\sum_j m_j/T_j}`.

The results can be adapted in case our input are not absolute values but the relative frequencies.
:math:`f_i` is the frequency of mutations and :math:`t_i` the frequency of nucleotides:

.. math::

    f_i = \frac{m_i}{\sum_j m_j};      t_i = \frac{T_i}{\sum_j T_j} \simeq \frac{T_i}{N}

(:math:`N` is the number of nucleotides, :math:`\sum_j T_j = N - 2 \cdot s`, where :math:`s` is the number of segments)

Then:

.. math::

   \overline{f_i} = \frac{f_i/t_i}{\sum_j f_j/t_j}

Proof:

.. math::

   \frac{f_i/t_i}{\sum_j f_j/t_j} = \frac{\frac{\frac{m_i}{\sum_j m_j}}{\frac{T_i}{\sum_j T_j}}}{\sum_k \frac{\frac{m_k}{\sum_j m_j}}{\frac{T_k}{\sum_j T_j}}} = \frac{\frac{m_i}{T_i} \cdot \frac{\sum_j T_j}{\sum_j m_j}}{\sum_k (\frac{m_k}{T_k} \cdot \frac{\sum_j T_j}{\sum_j m_j})} = \frac{\frac{m_i}{T_i} \cdot \frac{\sum_j T_j}{\sum_j m_j}}{\frac{\sum_j T_j}{\sum_j m_j} \cdot \sum_k \frac{m_k}{T_k}} = \frac{m_i / T_i}{\sum_k m_k/T_k}

.. [#context] The context is formed by the previous and posterior nucleotides.
