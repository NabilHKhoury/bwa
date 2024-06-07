# BWALIGN

This is a python tool designed for read alignment to a reference genome. It utilizes the seed and extend method to identify exact matches with a burrows-wheeler matching algorithm, then extends these exact matches across entire reads using a Needleman-Wunsch banded global alignment algorithm.

## Installation Instructions

Do:

`pip install bwalign`

If you do not have root access, you can include the user flag like so:

`pip install bwalign --user`

If on Linux/Unix, or any OS where bin is not on the PATH, you may need to add it manually. For `datahub.ucsd.edu` users, the command will be:

`export PATH=$PATH:/home/$USER/.local/bin` after running the above pip. Every time you close the terminal, you will need to rerun this.

## Usage

Usage is as follows:

`bwalign <reference genome path> <fastq file>`

Both the reference genome and single fastq file are required parameters.

This will align reads to your reference genome using our default scoring parameters for banded alignment, which are: match reward (2), mismatch penalty (-2), gap opening penalty (-2), band parameter (10).

If you wish to run on a test dataset from this project, do:

`git clone git@github.com:NabilHKhoury/bwalign.git`

into your desired directory of choice, then `cd bwalign/data`

And run `bwalign` on the reference genome `random_reference.fasta` and `random_reads.fastq`.

## Output/File Format

The output file will be a SAM file with the name `output.sam`. It will be outputted in the directory you ran `bwalign` from.

## Contributors

This repository was made in collaboration by Nabil Khoury, Yasmin A. Jaber, and Adrian Layer. Inspiration from [BWA-MEM](https://bio-bwa.sourceforge.net/bwa.shtml).