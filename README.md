# Major Project - Leftover Reads Inspector

This project was developed for my Major Project as part of my final year at Aberystwyth University in the academic year
2023-2024.

## Description

This tool allows users to investigate genome reads of any kind, although is intended focus on leftover
reads. These are the ones that are not used in an assembly. It provides functions to analyse their GC content,
k-mers, entropy, find similar sequences between files, and provide the ability to BLAST a file of sequences
against the SwissProt database. The tool can find the unmapped reads, by assembling with Megahit, and alignment
with Bowtie2. 

More information surrounding the format of the Data directory can be found within the README in the directory.

## Installation

The project files must first be downloaded, and can be done using the command below, once the files have been added
to GitHub:

`$ wget [link to this repository]`

This tool can be installed using the below command, when in the new downloaded folder:
`$ pip install .`

Or if in the directory above this one, use the below command:

`$ pip install [name of directory tool was downloaded into]`

## Usage

Once downloaded, the tool can be used across your system. This makes it useful to navigate into a directory containing
data to investigate. In order to assemble reads, the most basic command is:

`$ LeftoverReadsInspector -f [file to assemble] -assemble_and_find_unmapped True`

This functionality provides other options: `-assembly_file_name` which customises the name of the output file of contigs
by Megahit. `-mapped_reads_file_name` customises the name of the file of mapped reads outputted by Bowtie,
and `-unmapped_reads_file_name` customises the name of the file of unmapped reads outputted by Bowtie.

In order to find the GC content of a file of sequences, the command below should be used, and creates the output
with default options, which produces bar charts and histograms:

`$ LeftoverReadsInspector -f [file to inspect] -find_gc True`

Options can be used to change the output of this functionality, including the graphs.

You can use `$ LeftoverReadsInspector --help ` to see the full list of options for the tool.

### Acknowledgement

Created and authored by Thomas Collins. Suggestion of project and supervision provided by Dr Amanda Clare

