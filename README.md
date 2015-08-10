# OncodriveFML #

Recent years saw the development of methods to detect signals of positive selection in the pattern of somatic mutations in genes across cohorts of tumors, and the discovery of hundreds of driver genes. The next major challenge in tumor genomics is the identification of non-coding regions which may also drive tumorigenesis. We present OncodriveFML, a method that estimates the accumulated functional impact bias of somatic mutations in any genomic region of interest based on a local simulation of the mutational process affecting it. It may be applied to all genomic elements to detect likely drivers amongst them. Here, we describe the results of its application to identify driver genes, promoters and UTRs in several malignancies, and compare it to current generation methods. Furthermore, we demonstrate that OncodriveFML can discover signals of positive selection when only a small fraction of the genome, like a panel of genes, has been sequenced.

## Installation ##

OncodriveFML depends on Python 3 and some external libraries, [numpy](http://www.numpy.org/), [scipy](http://www.scipy.org/), [pandas](http://pandas.pydata.org/) and [statsmodels](http://statsmodels.sourceforge.net/).

The easiest way to install all this software stack is using the well known [Anaconda Python distribution](http://continuum.io/downloads#34).

Then to get OncodriveFML installed run the following command:

	$ pip install oncodrivefml

The first time that you run OncodriveFML it will download the genome reference from our servers. By default the downloaded datasets go to ~/.bgdata if you want to move this datasets to another folder you have to define the system environment variable BGDATA_LOCAL with an export command. 
The following command will show you the command help:

	$ oncodrivefml --help

	usage: oncodrivefml [-h] [-i INPUT_FILE] [-r REGIONS_FILE] [-t SIGNATURE_FILE]
                    [-s SCORE_FILE] [-o OUTPUT_FOLDER] [-n PROJECT_NAME]
                    [-mins MIN_SAMPLINGS] [-maxs MAX_SAMPLINGS]
                    [--cores CORES] [--debug] [--no-figures] [--drmaa DRMAA]
                    [--drmaa-max-jobs DRMAA_MAX_JOBS]
                    [--trace TRACE [TRACE ...]] [--geometric]

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_FILE, --input INPUT_FILE
                            Variants file (maf, vcf or tab formated)
      -r REGIONS_FILE, --regions REGIONS_FILE
                            Genomic regions to analyse
      -t SIGNATURE_FILE, --signature SIGNATURE_FILE
                            Trinucleotide signature file
      -s SCORE_FILE, --score SCORE_FILE
                            Tabix score file
      -o OUTPUT_FOLDER, --output OUTPUT_FOLDER
                            Output folder
      -n PROJECT_NAME, --name PROJECT_NAME
                            Project name
      -mins MIN_SAMPLINGS, --min_samplings MIN_SAMPLINGS
                            Minimum number of randomizations
      -maxs MAX_SAMPLINGS, --max_samplings MAX_SAMPLINGS
                            Maximum number of randomizations
      --cores CORES         Maximum CPU cores to use (default all available)
      --debug               Show more progress details
      --no-figures          Output only the tsv results file
      --drmaa DRMAA         Run in a DRMAA cluster using this value as the number
                            of elements to compute per job.
      --drmaa-max-jobs DRMAA_MAX_JOBS
                            Maximum parallell concurrent jobs
      --trace TRACE [TRACE ...]
                            Elements IDs to store files to trace and reproduce the
                            execution
      --geometric           Use geometric mean instead of arithmetic mean
      

