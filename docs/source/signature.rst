
.. _signature:

Signature
=========

The signature is an array that assigns a probability to
a single nucleotide mutation taking into account its context [#context]_.
It represents the chance of a certain mutation to occur within a context.

Check the different options for the signature in the
:ref:`configuration file <config signature>`.
In short, you can choose between not using any signature, using your own signature
or computing the signature from the mutations file.
Additionally, signatures can be grouped into different categories
(such as the sample).

The signature is computed count all the Single Nucleotide Polymorphisms
in the input file, taking into account their context.
The counts are used to compute a frequency
:math:`f_i = \frac{m_i}{M}` where :math:`M = \sum_j m_j`, and
:math:`m_i` represent the number of times that the mutation
:math:`i` with its context [#context]_ has been observed.


Optionally, the signature can be corrected taking into
account the frequency of trinucleotides in the
reference genome.
OncodriveFML introduces this feature because the
distribution of triplets is not expected to be constant.
When using the command line interface, OncodriveFML
does this correction automatically according to
the value passed in the flag ``--signature-correction``
(you can list all the options :ref:`using the help <help cmd>`).

.. important:: Signature correction is done
   using precomputed counts of whole genome
   and whole exome of HG19 reference genome.

   This counts might be similar for other human genomes
   but ensure that correction is not done
   genomes of other species.
   Check the `command line <inside cli>`_
   and `configuration file <config signature>`_.


More complex signatures
(e.g. using only mutations that map to the regions
under analysis, or normalizing by the frequency
of trinucleotides in specific regions of the genome)
can be computed using the `bgsignature package <https://bitbucket.org/bgframework/bgsignature>`_
and passed to OncodriveFML via the `configuration file <config signature>`_.

.. c

	Reasoning behind the correction
	-------------------------------


	Let's first take the conditional probability of a mutation (with contectx [#context]_)
	to occur given the number of those triplets in the region:
	:math:`p_i = p(m = i | T_i) = \frac{m_i}{T_i}`.

	Then, the normalized frequency of the mutation :math:`i` is:
	:math:`\overline{f_i} = \frac{m_i/T_i}{\sum_j m_j/T_j}`.

	The results can be adapted in case our inputs are not absolute values but relative frequencies.
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


----

.. [#context] The context is formed by the previous and posterior nucleotides.
