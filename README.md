#### itBins is a superfast cli-tool for the automated decontamination of metagenome-assembled genomes (MAGs)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/itbins/README.html)


## Reference

A preprint of the paper is available on [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.30.715291v1).

Please cite as:

Julian M. Künkel, Till L. V. Bornemann, Wei Xiu, Joern Starke, Tom L. Stach, André Rodrigues Soares, Jörg Schlötterer, Christin Seifert and Alexander J. Probst. Automated refinement of metagenomic bins and estimation of binning success using itBins. [https://doi.org/10.64898/2026.03.30.715291](https://doi.org/10.64898/2026.03.30.715291)


## Usage

```
usage: itBins [-b BIN_OVERVIEW_FILE] [-g SINGLE_COPY_GENE_FILE] [-d DASTOOL_BIN_NAME_COLUMN] [-u UBIN_BIN_NAME_COLUMN] [-p PREFIX] [-o [OUTPUT_FILE]] [--estimate [ESTIMATE_PATH]] [-s [SUMMARY_FILE]] [-t [TASK_FILE]]
              [--example-task-file [TASK_FILE]] [-i] [-q] [--log-to-stderr] [--no-output] [-h] [--version] [--manual]
              
-b --bins BIN_OVERVIEW_FILE					provide path to an overview file
-g --genes SINGLE_COPY_GENE_FILE			provide path to a singe copy gene file
-d --bin-DASTool DASTOOL_BIN_NAME_COLUMN	provide the column name for the uncurated bins
-u --bin-uBin UBIN_BIN_NAME_COLUMN			provide the column name for the curated bins
--scaff_prefix SCAFF_PREFIX					provide a prefix to attach to the names of scaffolds
--bin_prefix BIN_PREFIX						provide a prefix to attach to the names of bins
-p --prefix PREFIX							provide a prefix to attach to the names of bins and scaffolds
-o --output [OUTPUT_FILE]					provide file-path to save the modified table of scaffolds to
--estimate [ESTIMATE_PATH]					provide a path to save the binning success esitimation to
-s --summary [SUMMARY_FILE]					provide a file-path to save a summary file to. Default: "itBins_summary.tsv"
-t --task-file [TASK_FILE]					provide path to a task file (.json). Default: "./tasks.json".
--example-task-file [TASK_FILE]				provide path to create an example task file at, defaults to "./tasks.json"
-i --info --verbose							will print detailed info
-q --quiet									will not print progress messages
--log-to-stderr								redirect all printing to stderr
--no-output									supress all output 
-h --help									show help message and exit
--version									display version info
--manual									display manual
```

### Input Preparation

itBins requires two input files, the overview-file and the single-copy-gene-file.
File 1, the overview-file, is based on the following data:

* Length in bps for each contig (to be retrieved from fasta sequence information, we recommend [this script](https://github.com/ProbstLab/uBin-helperscripts/blob/master/bin/04_03fasta_length_individual.rb))
* Percent GC-content of each contig (to be retrieved from fasta sequence information, we recommend [this script](https://github.com/ProbstLab/uBin-helperscripts/blob/master/bin/04_02gc_count.rb))
* Coverage information for each contig, based on read mapping for quality controlled reads. This is average coverage of a contig per nucleotide. We recommend [BOWTIE2](https://github.com/BenLangmead/bowtie2), but other mappers may also be used. To caluclate the average coverage we recommend [this script](https://github.com/ProbstLab/uBin-helperscripts/blob/master/bin/04_01calc_coverage_v3.rb).
* Taxonomic information of each contig, which is obtained by searching aganinst a custom database (we recommend [this script](https://github.com/ProbstLab/uBin-helperscripts/blob/master/bin/07_00annotate_newdb.sh))
* Contig-to-bin assignments, i.e. which contig belongs to which bin,  usually available from binning tools or bin aggregation, e.g., using [DASTool](https://github.com/cmks/DAS_Tool)

The data must then be put together as follows, in tab-separated format (*e.g.*, using [this script](https://github.com/ProbstLab/uBin-helperscripts/blob/master/bin/08_00overview.sh)):



| scaffold  | length | gc   | coverage | taxonomy               | bin             |
|-----------|--------|------|----------|------------------------|-----------------|
| scaff\_01 | 12345  | 50.2 | 22.2     | bacteria;unclassified; | example\_bin\_1 |


File 2, the single-copy-gene-file, is based on the following data:

* Single-copy-gene (SCG) information (available by searching against the SCGs sequences, *e.g.*, using [this script](https://github.com/ProbstLab/uBin-helperscripts/blob/master/bin/08_02scg_metagenome.rb))

The SCG data must be collected as follows, in comma separated format. Contigs without SCGs should still have an entry, just with zeros. If you plan to use your own SCGs and not those provided in the uBin helperscripts, please see advanced use.


| scaffolds | b\_gene\_01 | b\_gene\_02 | b\_gene\_03 | b\_gene\_04 | b\_gene\_05 | ... |
|-----------|-------------|-------------|-------------|-------------|-------------|-----|
| scaff\_01 |           1 |           2 |           0 |           1 |           1 | ... |

itBins is meant to replace (for large datasets) and compliment the manual refinement tool [uBin](https://github.com/ProbstLab/uBin), which uses the same input files. Scripts for creating the input files are available here: [https://github.com/ProbstLab/uBin-helperscripts](https://github.com/ProbstLab/uBin-helperscripts). The preprint (see above) also includes details on the format of the input files.

### Output

* Contig-to-bin file with a new column for itBins' bin assignments
* Summary file with info on completeness and contamination statistics of each bin, flags for suspected eukaryotic bins
* Statistics on binning success (estimations)


### Usage example
Running itBins with default settings, and redirecting output to a log file.

```
itbins -b overview.txt -g SCGs.csv -o itbins_output.tsv -s itbins_summary.tsv &> itbins.log
```


## Advanced Usage

### Modifying the task file & custom marker gene sets
itBins performs a number of refinement tasks when processing each bin. By providing a task file, the order of these tasks and their parameters, and also general parameters of itBins, may be modified. A basic task file with default values can be produced as follows:

```
itbins --example-task-file
```

This creates a file called tasks.json. By modifying it it is possible to:

* overwrite flags
* set debugging parameters
* define custom single-copy-gene sets
* change what tasks are performed and their order
* change task parameters

The taskfile starts as follows:

```
"flags": {"d": "Bin",
          "u": "curated"},
"parameters": {"min_runs": 0,
               "max_runs": 0,
               "BSE_minimum_coverage": 7,
               "archeal_set": "default",
               "bacterial_set": "default"},
```

The "flags" section allows overriding two of the flags, `-d` and `-u`.
It is followed by the "parameters" section, which lists minruns (for debugging, skip that many bins) and maxruns (for debugging, stop processing after that many bins, 0: process all). The marker gene set parameters may be set to "default" or a comma separated list of gene names, the list must be enclosed in square brackets. Any set of genes may be used, but the names must match the column headers in the SCG-file.

```
"archeal_set": "default",
"bacterial_set": ["b_gene_01", "b_gene_02", "b_gene_03", ...],
```

In the tasks section all steps in the bin refinement process are listed. Depending on the exact task, differnt parameters may be adjusted. The 'NAME' of a task may be any name but must be unique in the task file, the 'task_identifier' determines what task is performed in the step.

```
"NAME": {"todo": "task_identifier",
         "parameter_1": value,
         "parameter_2": value,
         ...
        },
```

The tasks are performed in the order listed in the task file, by changing the order in the task file, the order of operations can be modified. The 'stop' task may be inserted to stop refinement at that point without having to delete tasks further down the list. This is useful for debugging a custom taskfile. Run itBins with your taskfile using the `-t` flag:

```
itbins -t taskfile.json -b overview.tsv -g SCGs.csv -o itbins_output.tsv -s itbins_summary.tsv &> itbins.log
```

### Binning success estimation
itBins estimates the success of binning based on coverage and marker genes. As a requirement, all contigs must be passed to itBins in the input files, not just binned contigs. For an overall estimation of binning success, itBins compares the number of binned contigs to the total number of contigs, *e.g.*, 345932/634933
and it’s fraction 0.54. At the same time, itBins also provides a value taking into account the average coverage and length of each contig resulting in a fraction, *e.g.*, 0.32.

For estimating the binned community members, itBins uses three marker-genes, bacterial *rpS3*, bacterial *gyrA* and archaeal *rpS3Ae*. Custom marker gene sets must include these three genes for the estimation to take place. *rpS3* and *rpS3Ae* represent bacteria and archaea repectively, *gyrA* is meant to be a canary value. If the results for *gyrA* strongly deviate from *rpS3*, the dataset may be skewed and the estimation inaccurate.

To calculate the binning success of the community, itBins firstly determines the number of marker genes that ended up in bins versus those that remained unbinned, *e.g.*, 312/503, as well as the fraction *e.g.*, 0.62. In a second estimate, itBins includes average coverage values of contigs on which the marker genes reside. This estimate reflects the different relative abundance of community members and provides a fraction only, *e.g.*, 0.81.

In total itBins calculates:

* The numbers of binned versus unbinned contigs and their fraction
* The fraction of binned contigs taking into account coverage and length
* The numbers of binned versus unbinned contigs with marker genes and their fraction
* The fraction of binned contigs with marker genes taking into account their average coverage


## Dependencies

* [python](https://www.python.org) 3.10.11
* [pandas](https://pandas.pydata.org) 1.4.2
* [numpy](https://www.numpy.org/) 1.21.5


## Installation through bioconda (recommended)

itBins is ready to:

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/itbins/README.html)

[Mamba](https://mamba.readthedocs.io/en/latest/index.html) is recommended over conda. First, create a new environment to contain your itBins installation:

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

For more details on running itBins see the 'Usage' and 'Advanced Usage' sections above.


## Manual Installation (if the most up-to-date version is required)

itBins is mirrored to [Github](https://www.github.com/ProbstLab/itBins) from it's home on [Codeberg](https://www.codeberg.org/JMK/itBins). Usually both will be up to date. Download the file `itBins.py`. Set up an environment with [Mamba](https://mamba.readthedocs.io/en/latest/index.html) or a similar tool (using conda is no longer recommended).

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

Then run the following command in your project's directory, substituting '[pathToItBins]' for wherever you downloaded `itBins.py` to.

```
python /[pathToItBins]/itBins.py --example-task-file
```

This creates a basic config file for running itBins, called tasks.json. If you have the input files ready, you can now run itBins:

```
python /[pathToItBins]/itBins.py -s -t tasks.json -b overview.txt -g SCGs.csv -o itBins_output.tsv &> itBins.log
```

The `-s` flag is used to create a summary file, which lists for example if a bin is suspected to be eukaryotic, and the scores before and after refinement. `-t` specifies the task file that should be used to specify parameters and task order. `-b` and `-g` specify the overview and single-copy-gene files respectively. `-o` specifies where to save output to. Saving output streams to logfile with `&>` is recommended.

For more details on running itBins see the 'Usage' and 'Advanced Usage' sections above.


## Troubleshooting

Most issues result from malformed input files. Make sure all contigs have one and only one entry in both the overview and the single-copy-gene file, and when concatenating files, avoid duplicate headers. If problems persist, feel free to [open an issue on Codeberg](https://codeberg.org/JMK/itBins/issues/new) or [reach out](https://www.uni-due.de/probst-lab/team.php).