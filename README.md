# OncodriveFML #

Recent years saw the development of methods to detect signals of positive selection in the pattern of somatic mutations in genes across cohorts of tumors, and the discovery of hundreds of driver genes. The next major challenge in tumor genomics is the identification of non-coding regions which may also drive tumorigenesis. We present OncodriveFML, a method that estimates the accumulated functional impact bias of somatic mutations in any genomic region of interest based on a local simulation of the mutational process affecting it. It may be applied to all genomic elements to detect likely drivers amongst them. OncodriveFML can discover signals of positive selection when only a small fraction of the genome, like a panel of genes, has been sequenced.

## License ##
OncodriveFML is made available to the general public subject to certain conditions described in its [LICENSE](LICENSE). For the avoidance of doubt, you may use the software and any data accessed through UPF software for academic, non-commercial and personal use only, and you may not copy, distribute, transmit, duplicate, reduce or alter in any way for commercial purposes, or for the purpose of redistribution, without a license from the Universitat Pompeu Fabra (UPF). Requests for information regarding a license for commercial use or redistribution of OncodriveFML may be sent via e-mail to innovacio@upf.edu.

## Installation ##

OncodriveFML depends on Python 3.4 and some external libraries. The easiest way to install all this software stack is using the well known [Anaconda Python distribution](http://continuum.io/downloads).

Then to get OncodriveFML installed first clone the repository and then install it using ``pip``:

        $ git clone git@bitbucket.org:bbglab/oncodrivefml.git
        $ cd oncodrivefml
        $ pip install .

The first time that you run OncodriveFML it will download the genome reference from our servers. By default the downloaded datasets go to ``~/.bgdata`` if you want to move this datasets to another folder you have to define the system environment variable BGDATA_LOCAL with an export command. 

The following command will show you the command help:

	$ oncodrivefml --help
    usage: oncodrivefml [-h] -i INPUT_FILE -r REGIONS_FILE -s SCORE_FILE
                        [-t SIGNATURE_FILE] [-o OUTPUT_FOLDER] [-n PROJECT_NAME]
                        [--statistic] [-mins MIN_SAMPLINGS] [-maxs MAX_SAMPLINGS]
                        [--samples-blacklist SAMPLES_BLACKLIST]
                        [--signature-ratio SIGNATURE_RATIO] [--no-figures]
                        [-D INDELS_FILE] [--indels-background INDELS_BACKGROUND]
                        [--cores CORES] [--debug] [--trace TRACE [TRACE ...]]
                        [--drmaa DRMAA] [--drmaa-max-jobs DRMAA_MAX_JOBS]
                        [--resume] [-q QUEUES]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    Mandatory options:
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
    
    General options:
      -o OUTPUT_FOLDER, --output OUTPUT_FOLDER
                            Output folder. Default to 'output'
      -n PROJECT_NAME, --name PROJECT_NAME
                            Project name
      --statistic           Statistic to use: amean, gmean, max
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
    
    Indels options:
      -D INDELS_FILE, --indels INDELS_FILE
                            Indels scores file
      --indels-background INDELS_BACKGROUND
                            Indels random background scores
    
    Execution options:
      --cores CORES         Maximum CPU cores to use (default all available)
      --debug               Show more progress details
      --trace TRACE [TRACE ...]
                            Elements IDs to store files to trace and reproduce the
                            execution
    
    Cluster options:
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
    
To run this example OncodriveFML needs all precomputed CADD scores that is a 80Gb file. It will be 
automatically download the first time that you run OncodriveFML, but if you want to speed up the download is better
if you first download it using our data package managment tool that is also installed when you install OncodriveFML.

Run this command to download the CADD scores file to the default bgdata folder `~/.bgdata`:
 
    $ bg-data -n 10 genomicscores cadd 1.3

Also if you want to speed up the download of the genome reference that is also needed, run this command:

    $ bg-data -n 10 datasets genomereference hg19

Now run OncodriveFML like this:

    $ oncodrivefml -i paad.txt.gz -r cds.regions.xz -t compute -s cadd.conf
    
    13:00:13 INFO: Loading regions
    13:00:30 INFO: Loading and mapping mutations
    13:00:30 INFO: [1 of 20763]
    13:00:44 INFO: [7333 of 20763]
    13:00:59 INFO: [14665 of 20763]
    13:01:12 INFO: [20763 of 20763]
    13:01:14 INFO: Computing signature
    13:01:17 INFO: Computing statistics
    13:01:19 INFO: [1 of 2560]
    13:02:17 INFO: [145 of 2560]
    13:02:41 INFO: [289 of 2560]
    13:03:04 INFO: [433 of 2560]
    13:04:11 INFO: [577 of 2560]
    13:05:38 INFO: [721 of 2560]
    13:05:38 INFO: [865 of 2560]
    13:05:38 INFO: [1009 of 2560]
    13:06:57 INFO: [1153 of 2560]
    13:06:57 INFO: [1297 of 2560]
    13:06:57 INFO: [1441 of 2560]
    13:06:57 INFO: [1585 of 2560]
    13:08:34 INFO: [1729 of 2560]
    13:08:34 INFO: [1873 of 2560]
    13:08:51 INFO: [2017 of 2560]
    13:08:51 INFO: [2161 of 2560]
    13:09:19 INFO: [2305 of 2560]
    13:09:42 INFO: [2449 of 2560]
    13:14:56 INFO: [2560 of 2560]
    13:14:56 WARNING: There are background positions without signature probability. We are using a probability of zero at these positions.
    13:14:56 WARNING: If you are computing the signature from the input file, most probable this means that you don't have enough mutations.
    13:14:56 WARNING: Try using a precomputed signature of a similar cancer type to improve the results.
    13:14:56 WARNING: The missing signatures are:
    13:14:56 WARNING:       ref: 'TAT' alt: 'TCT'
    13:14:56 WARNING:       ref: 'TAG' alt: 'TCG'
    13:14:56 WARNING:       ref: 'ATT' alt: 'AGT'
    13:14:56 WARNING:       ref: 'CAA' alt: 'CTA'
    13:14:56 WARNING:       ref: 'TTA' alt: 'TAA'
    13:14:56 INFO: Computing multiple test correction
    13:14:59 INFO: Creating figures
    13:15:08 INFO: Done
    
Ignore this warning that is due to the fact that we are using a small example dataset. You can browse the results in the `output` folder.