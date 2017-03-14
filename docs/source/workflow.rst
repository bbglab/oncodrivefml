How it works
============

This section will try to give an overview of
how OncodriveFML performs the analysis.

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
  --indels
                                  Type of analysis performed with the indels [required]:

                                  - *discard* to exclude indels from the analysis
                                  - *coding* TODO
                                  - *noncoding* TODO

                                  See :ref:`details about the command line interface <inside cli>`
                                  to find more information about this option.
  --sequencing
                                  Type of sequencing  [required]:

                                  - *genome*: for datasets that come from whole genome analysis
                                  - *exome*: for datasets that come from whole exome analysis
                                  - *other*: for other types of datasets

                                  See :ref:`details about the command line interface <inside cli>`
                                  to find more information about this option.
  -o, --output OUTPUT_FOLDER      Output folder. Default to regions file name
                                  without extensions.
  -c, --configuration CONFIG_FILE
                                  Configuration file. Default to
                                  'oncodrivefml.conf' in the current folder if
                                  exists or to ~/.bbglab/oncodrivefml.conf if
                                  not.
  --samples-blacklist SAMPLES_BLACKLIST
                                  Remove these samples when loading the input
                                  file.
  --generate-pickle               Run OncodriveFML to generate pickle files
                                  that could speed up future executions and
                                  exit.
  --debug                         Show more progress details
  --version                       Show the version and exit.
  -h, --help                      Show this message and exit.


If you prefer to call OncodriveFML from a Python script,
you can download the source code, install it and call the
:func:`~oncodrivefml.main.main` function.

.. note::

   You might have notice that the :func:`~oncodrivefml.main.main`
   function accepts less parameters than the command line
   interface. This is because the command line interface
   modifies some parameters in the configuration, while
   calling direct Python code does not.
   Check :ref:`what is modified by the command line interface <inside cli>`.

   This implies that you should adapt the
   :ref:`configuration file <project configuration>`
   to your needs.


The files
---------

Check the :ref:`different formats for
the input files<oncodrive file formats>`.

The configuration file is also a key part of the run,
and understanding how to adapt it to your parameters is important.
Check :ref:`this section <project configuration>`
to find more details about it.

Output file
^^^^^^^^^^^

Find information about the output :ref:`output files <output files>` section.

Workflow
--------

1. The first thing that is done by OncodriveFML is loading
   the configuration and creating the output folder if it does not exist.

   .. note::

      If you have not provided any output folder, OncodriveFML
      will create one in the current directory with the same name
      as the elements file (without extension).

   If the output folder exits, OncodriveFML checks whether a
   file with the expected output name exits and, if so, it does not
   run.

#. The regions file is loaded, and a tree with the intervals is created.
   This tree is used to find which mutations fall in the regions being
   analysed.

#. Load the mutations file and keep only the ones that fall into the regions
   being analysed.

#. Compute the signature (see the :ref:`signature <signature>` section).

#. Analyse each region separately (only the ones that have mutations).
   In each region the analysis is as follow:

   1. Compute the score of each of the observed mutations.

   #. Simulate the same number of mutations in the segments of the region under analysis.
      Save the scores of each of the simulated mutations.
      The simulation is done several times.

   #. Apply a predefined function to the observed scores and to each of the simulated
      groups of scores.
      Count how many times the simulated value is higher or equal than the observed.

   #. From these counts, compute the P value by dividing the counts by the number
      of simulations performed.

      .. warning::::

         As the statistical power is not infinite, the values carry an error.
         Due to this error, OncodriveFML does not provide P values of 0
         even if the counts are 0. OncodriveFML uses in those cases a count of 1.

   You can find more details in the :ref:`analysis section <analysis>`.

#. Join the results and perform a multiple test correction.
   The multiple test correction is only done for regions with
   mutations from, at least, two samples.
   ## TODO explain why

#. Do some checks which include:

    #TODO

#. Create the :ref:`output files <output files>`.
