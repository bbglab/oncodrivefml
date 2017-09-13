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

The variants file is a text file with, at least, 5 columns separated by a tab character (the header is required, but the order of the columns can change):

* Column CHROMOSOME: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column POSITION: Mutation position. A positive integer.
* Column REF: Reference allele [#refalt]_.
* Column ALT: Alternate allele [#refalt]_.
* Column SAMPLE: Sample identifier. Any alphanumeric string.
* Column CANCER_TYPE: Cancer type. Any alphanumeric string. Optional.
* Column SIGNATURE: User defined signature categories. Any alphanumeric string. Optional.

.. _files region format:

Regions file format
^^^^^^^^^^^^^^^^^^^

The regions file is a text file with, at least, 4 columns separated by a tab character
(the column order must be preserved):

* Column 1 [CHROMOSOME]: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column 2 [START]: Start position. A positive integer.
* Column 3 [STOP]: End position. A positive integer.
* Column 4 [STRAND]: Strand: ``+`` for positive, ``-`` for negative, ``.`` for unknown.
* Column 5 [ELEMENT]: Element identifier.
* Column 6 [SEGMENT]: Segment identifier. Optional column.
* Column 7 [SYMBOL]: Symbol, a different identifier for the element that will also be printed in the output file. Optional column.


Output file format
^^^^^^^^^^^^^^^^^^

OncodriveFML generates a tabulated file with the results with the
extension ".tsv".

Check the :ref:`output section <output files>` to find a detailed description
regarding the output.


.. [#refalt] The alleles consist on a single letter or a set of letters using A, C, G or T (upper case).
   Single Nucleotide Variants are indentified because both, REF and ALT contain only one letter.
   In Multi-Nucleotide Variants REF and ALT columns contain a set of letters of the same length.
   Insertions use ``-`` in the REF and a set of letters as ALT
   while deletions contain the set of deleted characters in the REF and ``-`` in the ALT columns.

