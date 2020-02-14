Files
=====


.. _oncodrive file formats:

File formats
------------

.. note::

   All the files can be compressed using GZIP (extension ".gz"), BZIP2 (extension ".bz2") or LZMA (extension ".xz")

.. _files input format:

Input file format
^^^^^^^^^^^^^^^^^

The variants file is a text file with, at least, 5 columns separated by a tab character
(the header is required, but the order of the columns can change):

* Column CHROMOSOME: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column POSITION: Mutation position. A positive integer.
* Column REF: Reference allele [#refalt]_.
* Column ALT: Alternate allele [#refalt]_.
* Column SAMPLE: Sample identifier. Any alphanumeric string.
* Column CANCER_TYPE: Cancer type. Any alphanumeric string. Optional.
* Column SIGNATURE: User defined signature categories. Any alphanumeric string. Optional.

Mutations are expected to be in the positive strand.

.. note:: OncodriveFML, although reading the SAMPLE column, it does
   not perform a per-sample analysis.

   A **by-sample** option can be enabled in the configuration file,
   in which only one mutation per sample is included in the analysis.
   More details in the :ref:`configuration section <per sample analysis>`.


.. _files region format:

Regions file format
^^^^^^^^^^^^^^^^^^^

The regions file is a text file with, at least, 4 columns separated by a tab character
(header is required):

* Column CHROMOSOME]: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column START: Start position. A positive integer.
* Column END: End position. A positive integer.
* Column ELEMENT: Element identifier. Can appear multiple times if the
  element is divided in *segments*.


.. important:: Analysis is perform element-wise.
   One single element can have multiple *segments*
   (even if you do not provide an identifier for them).

   It is also important that different segments of
   the same element do not overlap.


Optional columns are:

* Column STRAND: Strand: ``+`` for positive, ``-`` for negative, ``.`` for unknown.
* Column SEGMENT: Segment identifier. Optional column.
* Column SYMBOL: Symbol, a different identifier for the element that will also be printed in the output file. Optional column.



Signature file format
^^^^^^^^^^^^^^^^^^^^^

The signature file is a JSON file, where
pairs of key-values represent the changes
and the probabilities of those changes.

Changes are represented as
``AAA>C`` (reference triplet, ``>`` and alternate).

See the `bgsignature package <https://bitbucket.org/bgframework/bgsignature>`_
for more information on how to create such signatures.


Output file format
^^^^^^^^^^^^^^^^^^

OncodriveFML generates a tabulated file with the results with the
extension ".tsv.gz". It is compressed with gzip.

Check the :ref:`output section <output files>` to find a detailed description
regarding the output.

----

.. [#refalt] The alleles consist on a single letter or a set of letters using A, C, G or T (upper case).
   Single Nucleotide Variants are indentified because both, REF and ALT contain only one letter.
   In Multi-Nucleotide Variants REF and ALT columns contain a set of letters of the same length.
   Insertions use ``-`` in the REF and a set of letters as ALT
   while deletions contain the set of deleted characters in the REF and ``-`` in the ALT columns.

