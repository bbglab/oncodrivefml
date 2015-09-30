# OncodriveFML #

Recent years saw the development of methods to detect signals of positive selection in the pattern of somatic mutations in genes across cohorts of tumors, and the discovery of hundreds of driver genes. The next major challenge in tumor genomics is the identification of non-coding regions which may also drive tumorigenesis. We present OncodriveFML, a method that estimates the accumulated functional impact bias of somatic mutations in any genomic region of interest based on a local simulation of the mutational process affecting it. It may be applied to all genomic elements to detect likely drivers amongst them. Here, we describe the results of its application to identify driver genes, promoters and UTRs in several malignancies, and compare it to current generation methods. Furthermore, we demonstrate that OncodriveFML can discover signals of positive selection when only a small fraction of the genome, like a panel of genes, has been sequenced.

## License ##
OncodriveFML is made available to the general public subject to certain conditions described in its [LICENSE](LICENSE). For the avoidance of doubt, you may use the software and any data accessed through UPF software for academic, non-commercial and personal use only, and you may not copy, distribute, transmit, duplicate, reduce or alter in any way for commercial purposes, or for the purpose of redistribution, without a license from the Universitat Pompeu Fabra (UPF). Requests for information regarding a license for commercial use or redistribution of OncodriveFML may be sent via e-mail to innovacio@upf.edu.

## Installation ##

OncodriveFML depends on Python 3.4 and some external libraries. The easiest way to install all this software stack is using the well known [Anaconda Python distribution](http://continuum.io/downloads#34).

Then to get OncodriveFML installed run the following command:

	$ pip install oncodrivefml

The first time that you run OncodriveFML it will download the genome reference from our servers. By default the downloaded datasets go to ~/.bgdata if you want to move this datasets to another folder you have to define the system environment variable BGDATA_LOCAL with an export command. 
The following command will show you the command help:

	$ oncodrivefml --help
    usage: oncodrivefml [-h] -i INPUT_FILE -r REGIONS_FILE -s SCORE_FILE
                        [-t SIGNATURE_FILE] [-o OUTPUT_FOLDER] [-n PROJECT_NAME]
                        [--geometric] [-mins MIN_SAMPLINGS] [-maxs MAX_SAMPLINGS]
                        [--samples-blacklist SAMPLES_BLACKLIST]
                        [--signature-ratio SIGNATURE_RATIO] [--no-figures]
                        [-D INDELS_FILE] [--indels-background INDELS_BACKGROUND]
                        [--cores CORES] [--debug] [--trace TRACE [TRACE ...]]
                        [--drmaa DRMAA] [--drmaa-max-jobs DRMAA_MAX_JOBS]
                        [--resume] [-q QUEUES]
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_FILE, --input INPUT_FILE
                            Variants file
      -r REGIONS_FILE, --regions REGIONS_FILE
                            Genomic regions to analyse
      -s SCORE_FILE, --score SCORE_FILE
                            Substitutions scores file
      -t SIGNATURE_FILE, --signature SIGNATURE_FILE
                            Trinucleotide signature file. Use 'compute' to compute
                            it from the whole file, use 'none' if you don't want
                            to use signature.
      -o OUTPUT_FOLDER, --output OUTPUT_FOLDER
                            Output folder. Default to 'output'
      -n PROJECT_NAME, --name PROJECT_NAME
                            Project name
      --geometric           Use geometric mean instead of arithmetic mean
      -mins MIN_SAMPLINGS, --min-samplings MIN_SAMPLINGS
                            Minimum number of randomizations (default is 10k).
      -maxs MAX_SAMPLINGS, --max-samplings MAX_SAMPLINGS
                            Maximum number of randomizations (default is 100k).
      --samples-blacklist SAMPLES_BLACKLIST
                            Remove this samples when loading the input file
      --signature-ratio SIGNATURE_RATIO
                            Folders with one fold change vector per element to
                            multiply to the signature probability
      --no-figures          Output only the tsv results file
      -D INDELS_FILE, --indels INDELS_FILE
                            Indels scores file
      --indels-background INDELS_BACKGROUND
                            Indels random background scores
      --cores CORES         Maximum CPU cores to use (default all available)
      --debug               Show more progress details
      --trace TRACE [TRACE ...]
                            Elements IDs to store files to trace and reproduce the
                            execution
      --drmaa DRMAA         Run in a DRMAA cluster using this value as the number
                            of elements to compute per job.
      --drmaa-max-jobs DRMAA_MAX_JOBS
                            Maximum parallell concurrent jobs
      --resume              Resume a DRMAA execution
      -q QUEUES             DRMAA cluster queues

      
## File formats ##

**TIP**:  All the files can be compressed using GZIP (extension ".gz"), BZIP2 (ext. ".bz2") or LZMA (ext. ".xz")

### Input file format ###

The variants file is a text file with 5 columns separated by a tab character (without any header or comment):

* Column 1: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column 2: Mutation position. A positive integer.
* Column 3: Reference allelle. A single letter: A, C, G or T (upper case)
* Column 4: Alternate allelle. A single letter: A, C, G or T (upper case)
* Column 5: Sample identifier. Any alphanumeric string.
      
### Regions file format ###

The regions file is a text file with 4 columns separated by a tab character (without any header or comment):

* Column 1: Chromosome. A number between 1 and 22 or the letter X or Y (upper case)
* Column 2: Start position. A positive integer.
* Column 3: End position. A positive integer.
* Column 4: Element identifier.

## Run an example ##

Download and extract example files:

    $ wget https://bitbucket.org/bbglab/oncodrivefml/downloads/oncodrivefml-examples.tar.gz
    $ tar xvzf oncodrivefml-examples.tar.gz
    
Run OncodriveFML like this:

    $ oncodrivefml -i gbm.txt.gz -r cds.regions.xz -t compute -s cadd.conf
    
Browse the results at the `output` folder.