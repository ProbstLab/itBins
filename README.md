itBins is a superfast cli-tool for the automated curation of metagenome-assembled genomes (MAGs).

### Prerequesites:

* [python](https://www.python.org)
* [pandas](https://www.pandas.pydata.org)
* [numpy](https://www.numpy.org/)

### Install:

Itbins is mirrored to [Github](https://www.github.com/ProbstLab/itBins) from it's home on [Codeberg](https://www.codeberg.org/JMK/itBins). Get it from either. Set up an environment with mamba or a similar tool (using conda is no longer recommended).

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

Then run the following command in your project's directory, substituting for wherever you keep itBins.

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
coming soon ...

### Configuration:
You can configure itBins by changing the tasks.json file, adjusting the order of tasks and their parameters. It is also possible to remove tasks.

### Roadmap:
coming soon ...

### Citation:
A publication is coming soon, until then please cite the codeberg repository.

