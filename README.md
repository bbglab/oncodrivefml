# OncodriveFML

Recent years saw the development of methods to detect signals of positive
selection in the pattern of somatic mutations in genes across cohorts of tumors,
and the discovery of hundreds of driver genes. The next major challenge in tumor
genomics is the identification of non-coding regions which may also drive
tumorigenesis. We present OncodriveFML, a method that estimates the accumulated
functional impact bias of somatic mutations in any genomic region of interest
based on a local simulation of the mutational process affecting it. It may be
applied to all genomic elements to detect likely drivers amongst them.
OncodriveFML can discover signals of positive selection when only a small
fraction of the genome, like a panel of genes, has been sequenced.

## License

OncodriveFML is made available to the general public subject to certain
conditions described in its [LICENSE](LICENSE). For the avoidance of doubt,
you may use the software and any data accessed through UPF software for
academic, non-commercial and personal use only, and you may not copy,
distribute, transmit, duplicate, reduce or alter in any way for commercial
purposes, or for the purpose of redistribution, without a license from the
Universitat Pompeu Fabra (UPF). Requests for information regarding a license for
commercial use or redistribution of OncodriveFML may be sent via e-mail to
innovacio@upf.edu.

## Usage

OncodriveFML is meant to be used through the command line.

By default, OncodriveFML is prepared to analyse mutations using HG19 reference
genome. For other genomes, update the [configuration] accordingly.

[configuration]: https://oncodrivefml.readthedocs.io/en/latest/configuration.html

### Running OncodriveFML without installation (using Docker)

You can run OncodriveFML without having to install anything in your machine if
you have [Docker installed](https://docs.docker.com/engine/install/).

This is how you would run the example included in this repository:

```bash
docker run --rm -i \
		-v ${BGDATA_LOCAL:-${HOME}/.bgdata}:/root/.bgdata \
      -v $(pwd)/example:/data \
      --workdir /data \
		bbglab/oncodrivefml:2.5.0 \
         -i paad.txt.gz -e cds.tsv.gz --signature-correction wx --seed 123 --force
```

`-v ${BGDATA_LOCAL:-${HOME}/.bgdata}:/root/.bgdata` will allow the docker
container to see the contents of your *bgdata* directory as defined by the
environment variable `BGDATA_LOCAL` (or if it is not defined, the default
`~/.bgdata`).

`-v $(pwd)/example:/data` will allow the docker container to see the example data
in `./example`. You would need to replace `$(pwd)/example` by the folder where
you have your own data.

`--workdir /data` will set the working directory to the data folders you
specified before.

The results will be saved in a folder named `cds` under the `./example` folder.

### Installation

OncodriveFML can work with the Python versions 3.8 up to 3.11 (included).

The easiest way to install all this software stack is using the well known
[Anaconda Python distribution](http://continuum.io/downloads>)

```bash
conda install -c bbglab oncodrivefml
```

OncodriveFML can also be installed using `pip`:

```bash
pip install oncodrivefml
```

Finally, you can get the latest code from the repository and install it in development mode:

```bash
git clone https://github.com/bbglab/oncodrivefml.git
cd oncodrivefml
make build-dev
source .venv/bin/activate
oncodrivefml --help
```

> [!NOTE]
> The first time that you run OncodriveFML it will download the genome reference
> from our servers. By default the downloaded datasets go to `~/.bgdata`, but if
> you want to move these datasets to another folder you have to define the
> system environment variable `BGDATA_LOCAL` with an export command.

#### Running the example

Download and extract the example files (if you cloned the repository skip this step):

```bash
wget https://github.com/bbglab/oncodrivefml/archive/refs/tags/2.5.0.tar.gz
tar xvzf 2.5.0.tar.gz
```

To run this example OncodriveFML needs all the precomputed *CADD* scores, that
is a **17Gb file**, that will be downloaded automatically, together with the
reference genome.

> [!WARNING]
> The CADD scores are originally from http://cadd.gs.washington.edu/ and are
> freely available for all non-commercial applications. If you are planning on
> using them in a commercial application, please contact them at
> http://cadd.gs.washington.edu/contact.

To run the example, we have included a bash script (`run.sh`)
that will execute OncodriveFML. The script should be executed in
the folder where the files have been extracted:

```bash
cd oncodrivefml-2.5.0/example
./run.sh
```

The results will be saved in a folder named `cds`.

> [!NOTE]
> It might fail to run in macOS. We recommend you to run it using the Docker image instead.
> See the section `Running OncodriveFML without installation (using Docker)` for details.

### Configuring OncodriveFML

Although OncodriveFML includes a predefined configuration file, it is highly
recommended to create one yourself. In fact, if you are interested in using a
reference genome other than HG19, or a score other than CADD 1.0, it is
mandatory. See the documentation for the [configuration] for more details.


### Documentation

Find OncodriveFML documentation in
[ReadTheDocs](http://oncodrivefml.readthedocs.io/en/latest/).

You can also compile the documentation yourself using
[Sphinx](http://www.sphinx-doc.org/en/stable/), if you have cloned the
repository. To do so, run the following command:

```bash
make docs
open docs/build/html/index.html 
```
