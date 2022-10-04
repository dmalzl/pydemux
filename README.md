# pydemux
![pypi](https://img.shields.io/badge/pypi-v1.0.0-blue)
![python-version](https://img.shields.io/badge/Python->=3.9-blue)
![stable-version](https://img.shields.io/badge/version-1.0.0-blue)
An easy to use python tool for fast demultiplexing of SAM/BAM formated raw sequence data

# Installation
```
pip install pydemux
```
or clone the repository
```
git clone git@github.com:dmalzl/pydemux.git
```
and run
```
cd pydemux
pip install .
```
verify the installation by typing
```
pydemux -h
```

# Basic usage
pydemux can be used to demultiplex sequence data in SAM/BAM format for either single or paired-end reads (in case of paired-end the reads have to be interleaved). Basic command for this are listed below (you can test these using the files in the data directory)
```
pydemux single -b single_barcodes.tsv -o demux/ -s demux_stats.tsv single.bam
```
```
pydemux paired -b paired_barcodes.tsv -o demux/ -s demux_stats.tsv paired.bam
```
Note that for writing the files to a specific output directory the `/` needs to be included for the `-o/--output_prefix` commandline argument or they will be written to the current work directory. In order to compress the output one can also add the `-gz/--gzip` commandline argument to the command.

# Changing looked up SAMfile tags from which the barcodes are read
By default the algorithm looks for the barcode in the `BC` tag of each read in case of `pydemux single` and in the `BC` and `B2` tag of the paired-end reads in case of `pydemux paired`. This can be changed using either the `-t/--bctag` or `-t1/--bc1tag` and/or `-t2/--bc2tag` for `pydemux single` or `pydemux paired` respectively.

# Speeding up processing
In order to speed up demultiplexing the algorithm can be run concurrently using the `-p/--processes` argument.

# Optimizing yield
As with every sequencing-based data type, the barcodes are also prone to include sequencing errors. In order to optimize the read yield for each sample the algorithm allows for a given number of mismatches between the true and the sequenced barcodes which can be set with `-m/--mismatches`. By default only exact matches will be assigned. If you want to allow for one or more mismatches please make sure that the number of allowed mismatches does not exceed half of the minimum pairwise Hamming distance of all true barcodes minus 1 (i.e. `min(pairwise_hamming_distance(true_barcodes)) // 2 - 1`) since otherwise reads might be wrongly assigned
