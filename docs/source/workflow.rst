How it works
============

This section will try to give an overview of
how OncodriveFML carries on the analysis.

The command line interface
--------------------------

.. _help cmd:

By typing ``oncodrivefml -h`` you will have a brief
description of how to use OncodriveFML:



Options:
  -i, --input MUTATIONS_FILE      Variants file  [required]
                                  (:ref:`see format <files input format>`)
  -e, --elements ELEMENTS_FILE    Genomic elements to analyse  [required]
                                  (:ref:`see format <files input format>`)

  -o, --output OUTPUT_FOLDER      Output folder. Default to regions file name
                                  without extensions.
  -c, --configuration CONFIG_FILE
                                  Configuration file. Default to
                                  'oncodrivefml_v2.conf' in the current folder if
                                  exists or to ~/.config/bbglab/oncodrivefml_v2.conf if
                                  not.
  --samples-blacklist SAMPLES_BLACKLIST
                                  Remove these samples when loading the input
                                  file.
  --signature SIGNATURE           File with the signatures to use

                                  See :ref:`details about the command line interface <inside cli>`
                                  to find more information about this option.

  --signature-correction [wg|wx]  Correct the computed signutares by genomic
                                  or exomic signtures. Only valid for human
                                  genomes (``hg19`` and ``hg38``)

                                  - *wg*: correction using whole genome counts
                                  - *wx*: correction using whole exome counts

                                  See :ref:`details about the command line interface <inside cli>`
                                  to find more information about this option.

  --no-indels                     Discard indels in your analysis
  --cores INTEGER                 Cores to use. Default: all
  --seed INTEGER                  Set up an initial random seed to have
                                  reproducible results

  --debug                         Show more progress details
  --version                       Show the version and exit.
  -h, --help                      Show this message and exit.





The files
---------

Input files
^^^^^^^^^^^

OncodriveFML makes use of three files:

Variants
   Also named as input.
   This file contains the observed mutations for the analysis.

Regions
   File containing the regions for the analysis.
   Only mutations that fall in these regions are analysed
   and only the genomic positions defined in this file are used
   for the simulation.

   You can define your own regions file
   based on your criteria. You can check
   an example of a regions file
   downloading `our example <https://bitbucket.org/bbglab/oncodrivefml/downloads/>`_.

   .. warning::

      It is not recommended to mix coding and
      non-coding regions in your regions file.
      In fact this will likely produce artifacts
      in the results as coding and non-coding regions
      of the genome have a very different functional
      impact scores. A good set of genomic regions should
      include elements that share biological functions
      (e.g. CDS, UTRs, promoters, enhancers, etc.).


Check the :ref:`formats for
the input files<oncodrive file formats>`.

Configuration
   The configuration file is also a key part of the run,
   and understanding how to adapt it to your needs is important.
   Check :ref:`this section <project configuration>`
   to find more details about it.

Output files
^^^^^^^^^^^^

Find information about the output :ref:`output files <output files>` section.

Workflow
--------

1. The first thing that is done by OncodriveFML is to load
   the configuration file.

2. The output is checked. The default behaviour is that
   OncodriveFML creates an output folder in the current directory
   with the same name as the elements file (without extension).

   If an output is provided and it exists and is a folder,
   OncodriveFML checks whether a
   file with the expected output name exits and, if so, it does not
   run. Otherwise, it assumes it is a path name an uses that as output.

   .. note::

      If the output does not exits, OncodriveFML only computes
      the tsv file with the results and skips the plots.


#. The regions file is loaded, and a tree with the intervals is created.
   This tree is used to find which mutations fall in the regions being
   analysed.

#. Loads the mutations file and keeps only the ones that fall into the regions
   being analysed.

#. Computes the signature (see the :ref:`signature <signature>` section),
   if not provided as an external file.

#. Analyses each region separately (only the ones that have mutations).
   In each region the analysis is as follow:

   1. Computes the score of each of the observed mutations.

   #. Simulates the same number of mutations in the segments of the region under analysis.
      Save the scores of each of the simulated mutations.
      The simulation is done several times.

   #. Applies a predefined function to the observed scores and to each of the simulated
      groups of scores.
      Counts how many times the simulated value is higher than, or equal to, the observed.

   #. From these counts, computes a P-value by dividing the counts by the number
      of simulations performed.

      .. warning::::

         As the statistical power is not infinite, the values carry an error.
         Due to this error, OncodriveFML does not provide P values of 0
         even if the counts are 0. OncodriveFML uses in those cases a count of 1.

   You can find more details in the :ref:`analysis section <analysis>`.

#. Joins the results and performs a multiple test correction.
   The multiple test correction is only done for regions with
   mutations from at least two samples.

   .. todo explain why

#. Creates the :ref:`output files <output files>`.

#. Checks that the output file does not contain
   missing or repeated genomic regions.
