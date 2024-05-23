# BWALIGN

This is a python tool designed for read alignment to a reference genome. It utilizes the seed and extend method to identify exact matches with a burrows-wheeler matching algorithm, then extends these exact matches across entire reads using a smith-waterman affine alignment algorithm.

## Installation Instructions

**TBC**

## Basic Usage

Basic usage is:

`bwalign <reference genome path> <fastq file>`

This will align reads to your reference genome using the default scoring parameters for affine alignment, which are: **TBD**

To run on an example dataset in this repo, do:

`bwalign foo.gz bar.fq`

## Complete Usage

Only the reference genome and single fastq file are required parameters. Otherwise, users can input the scoring parameters for the affine alignment portion of the tool.

* `-m MATCHREWARD`, `--match MATCHREWARD`: this defines the addition to total score upon finding a match in the alignment.
* `-i OPENPEN`, `--open OPENPEN`: defines the penalty to total score upon opening a gap in the alignment.
* `-e EXTPEN`, `--extension EXTPEN`: defines penalty to total score for extending a gap in the alignment.
* `-x MMPEN`, `--mismatch MMPEN`: defines penalty to total score for a mismatch in the alignment.
* `-o FILE`, `--output FILE`: writes the output to a file. By default, output writes to stdout.

## File Format

The output file will be a SAM file.

## Contributors

This repository was made in collaboration by Nabil Khoury, Yasmin Jaber, and Adrian Layer. Inspiration from [BWA-MEM](https://bio-bwa.sourceforge.net/bwa.shtml).