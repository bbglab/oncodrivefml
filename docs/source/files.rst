Files
=====


.. _oncodrive file formats:

File formats
------------

.. note::

   All the files can be compressed using GZIP (extension ".gz"), BZIP2 (ext. ".bz2") or LZMA (ext. ".xz")


Input file format
^^^^^^^^^^^^^^^^^

The variants file is a text file with, at least, 5 columns separated by a tab character (header is required, but the order of the columns can change):

* Column CHROMOSOME: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column POSITION: Mutation position. A positive integer.
* Column REF: Reference allelle. A single letter: A, C, G or T (upper case)
* Column ALT: Alternate allelle. A single letter: A, C, G or T (upper case)
* Column SAMPLE: Sample identifier. Any alphanumeric string.

The variants file can contain more columns e.g. the cancer type. The more columns it contains, the more time it will take to read the file.


Regions file format
^^^^^^^^^^^^^^^^^^^

The regions file is a text file with, at least, 4 columns separated by a tab character
(header, between ``[]``, is optional, but the column order is fixed if header is not present):

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
