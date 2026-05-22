itBins is a superfast cli-tool for the automated curation of metagenome-assembled genomes (MAGs).

### Prerequesites:

* [python](https://www.python.org)
* [pandas](https://pandas.pydata.org)
* [numpy](https://www.numpy.org/)

### Install through bioconda (recommended in general):

itBins is ready to:

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/itbins/README.html)

```
mamba create -n itBinsEnv -c bioconda itbins
```

Activate the environment.

```
mamba activate itBinsEnv
```

Then run the following command in your project's directory.

```
itbins --example-task-file
```

This creates a basic config file for running itBins, called tasks.json. If you have the input files ready, you can now run itBins:

```
itbins -s -t tasks.json -b overview.txt -g SCGs.csv -o itBins_output.tsv > itBins.log
```

### Manual Install (recommended if the most up-to-date version is required):

itBins is mirrored to [Github](https://www.github.com/ProbstLab/itBins) from it's home on [Codeberg](https://www.codeberg.org/JMK/itBins). Usually the Codeberg repo will be the most recent and up to date. Set up an environment with mamba or a similar tool (using conda is no longer recommended).

```
mamba create -n itBinsEnv python=3.10.11 pandas=1.4.2 numpy=1.21.5
```

Activate the environment.

```
mamba activate itBinsEnv
```

Confirm that the Versions for python, pandas and numpy are correct.

```
mamba list
```

Then run the following command in your project's directory, substituting for wherever you downloaded itBins to.

```
python /[pathToItBins]/itBins.py --example-task-file
```

This creates a basic config file for running itBins, called tasks.json. If you have the input files ready, you can now run itBins:

```
python /[pathToItBins]/itBins.py -s -t tasks.json -b overview.txt -g SCGs.csv -o itBins_output.tsv > itBins.log
```

### Advanced Use:
coming soon ...

### Input Prep:
itBins uses the same input uBin uses. You can check out the [uBin-helperscripts](https://github.com/ProbstLab/uBin-helperscripts) for that. A python based version of these scripts will be made available in the near future. The preprint (see below) also includes details on the format of the input files.

### Configuration:
You can configure itBins by changing the tasks.json file, adjusting the order of tasks and their parameters. It is also possible to remove tasks, or force an early stop wit an inserted stop task while tinkering with your custom taskfile. Custom single-copy-gene sets can be configured in the taskfile aswell.

### Roadmap:
coming soon ...

### Citation:
A preprint of the paper is available on [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.30.715291v1).

