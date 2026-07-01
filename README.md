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
python /[pathToItBins]/itBins.py -s -t tasks.json -b overview.txt -g SCGs.csv -o itBins_output.tsv &> itBins.log
```

The '-s' flag is used to create a summary file, which lists for example if a bin is suspected to be eukaryotic, and the scores before and after refinement. '-t' specifies the task file that should be used to specify parameters and task order. '-b' and '-g' specify the overview and single-copy-gene files respectively. '-o' specifies where to save output to. Saving output streams to logfile with '&>' is recommended.


### Input Prep:
itBins uses the same input uBin uses. You can check out the [uBin-helperscripts](https://github.com/ProbstLab/uBin-helperscripts) for that. A python based version of these scripts will be made available in the near future. The preprint (see below) also includes details on the format of the input files.

To generate the input files, the following data needs to be collected: Length & GC-content for each contig, both of which can be trivially extracted from the fasta files. The reads also need to be mapped against the contigs to obtain coverage data, we recommend BOWTIE2, but other mappers may also be used. Contig-to-bin assignments will be available from binning or bin aggregation, e.g., using DASTool. Taxonomic information is obtained by blasting aganinst a custom database (see the uBin-helperscripts for details). The data can then be put together as follows, in tab-separated format:


| scaffold  | length | gc   | coverage | taxonomy               | bin             |
|-----------|--------|------|----------|------------------------|-----------------|
| scaff\_01 | 12345  | 50.2 | 22.2     | bacteria;unclassified; | example\_bin\_1 |


In addition, single-copy-gene information will also be needed, it is similarily available by blasting aganst the SCGs sequences. The SCG data can be collected as follows, in comma separated format. Contigs without SCGs should still have an entry, just with just zeros. The gene names must match the custom set, if one is employed.


| scaffolds | b\_gene\_01 | b\_gene\_02 | b\_gene\_03 | b\_gene\_04 | b\_gene\_05 | ... |
|-----------|-------------|-------------|-------------|-------------|-------------|-----|
| scaff\_01 |           1 |           2 |           0 |           1 |           1 | ... |


### Advanced Use:
#### Modifying the task file & custom marker gene sets
The flags section allows setting two of the flags that define column names. It is followed by the parameters section, which lists minruns and maxruns, both of which should be set to 0, and the gene set parameters, which may be set to default or a comma separated list of gene names, the list must be enclosed in square brackets. Any set of genes may be used, but the names must match the column headers in the SCGs.csv file.

```
"archeal_set": "default",
"bacterial_set": ["b_gene_01", "b_gene_02", "b_gene_03", ...],
```

The taskfile lists all tasks to perform in the tasks section, with all their parameters in clear text format. NAME may be any name but must be unique in the task file.

```
"NAME": {"todo": "task_identifier",
         "parameter_1": value,
         "parameter_2": value,
         ...
        },
```

The tasks are performed in the order listed in the task file, by changing the order in the task file the order of operations can be modified. The stop task may be inserted to stop refinement with it without having to delete further tasks.


#### Binning success estimation
Rank abundance curves are calculated for three genes, bacterial rpS3, bacterial gyrA and archaeal rpS3Ae. The fractions of binned and unbinned genes in the whole dataset and in the top 70% of coverage are reported. Contigs with a coverage < 7 are excluded, and the whole assembly, including unbinned contigs musst be passed into itBins, otherwise the result will be skewed. Custom marker gene sets must include the three listed genes for the estimation to take place. If the fraction of binned gyrA is significantly different from the fraction of binned rpS3, the estimation may not be accurate.


### Citation:
A preprint of the paper is available on [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.30.715291v1).

