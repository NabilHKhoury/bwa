# BWALIGN

This is a python tool designed for read alignment to a reference genome. It utilizes the seed and extend method to identify exact matches with a burrows-wheeler matching algorithm, then extends these exact matches across entire reads using a smith-waterman affine alignment algorithm.

## Installation Instructions

Do:

`pip install bwalign`

If you do not have root access, you can include the user flag like so:

`pip install bwalign --user`

## Usage

Usage is as follows:

`bwalign <reference genome path> <fastq file>`

Both the reference genome and single fastq file are required parameters.

This will align reads to your reference genome using the default scoring parameters for affine alignment, which are: match reward (2), mismatch penalty (-1), gap opening penalty (-2), gap extension penalty (-1).

If you wish to run 

On a test dataset from this project, do:

`git clone git@github.com:NabilHKhoury/bwalign.git`

into your desired directory of choice, then `cd bwalign/data`

And run `bwalign` on the reference genome `random_reference.fasta` and `random_reads.fastq`.

## Output/File Format

The output file will be a SAM file with the name `output.sam` (will be changed in future releases). It will be outputted in the directory you ran `bwalign` from.

## Contributors

This repository was made in collaboration by Nabil Khoury, Yasmin A. Jaber, and Adrian Layer. Inspiration from [BWA-MEM](https://bio-bwa.sourceforge.net/bwa.shtml).