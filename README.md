itBins is a superfast cli-tool for the automated curation of metagenome-assembled genomes (MAGs).

### Prerequesites:

* [python](https://www.python.org)
* [pandas](https://pandas.pydata.org)
* [numpy](https://numpy.org/)

### Install:

Set up an environment with conda or a similar tool

```
conda create -n itBinsEnv python=3.10.11 pandas=1.4.2 numpy=1.21.5
```

Activate the environment

```
conda activate itBinsEnv
```

Then run the following command in your project's directory, substituting for wherever you keep itBins

```
python /[pathToItBins]/itBins_0_6_0.py --example-task-file
```

This creates a basic config file for running itBins, called tasks.json. If you have the input files ready, you can now run itBins:

```
python /[pathToItBins]/itBins_0_6_0.py -s -t tasks.json -b overview.txt -g SCGs.csv -o itBins_output.tsv > itBins.log
```

### Advanced Use:
coming soon ...

### Input Prep:
coming soon ...

### Configuration:
coming soon ...

### Roadmap:
coming soon ...

### Citation:
coming soon ...

