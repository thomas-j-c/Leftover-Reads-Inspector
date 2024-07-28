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

## Results

These are sample results taken from applying the utilities within the tool to datasets of the Staphylococcus aureus bacteria. 

![hist_hamming_distances_reads_in_assembly_vs_SRR022](https://github.com/user-attachments/assets/fdb5eb4b-74ae-42d5-bcee-7c25348301f5)
This shows the Hamming distance between reads that were used in the program's assembly, and the whole file of reads. 

![gc_content_hist_for_all](https://github.com/user-attachments/assets/e9eb6633-c7a4-4696-97bd-bca5568c4940) All reads.
![gc_cont_hist_mapped](https://github.com/user-attachments/assets/94be9e1e-9db6-4cb5-88a1-08889ed5e4eb) Mapped reads.
![gc_cont_hist_unmapped](https://github.com/user-attachments/assets/7f378866-18a1-4bda-89d6-d4c887eb163e) Unmapped reads.

This histogram shows the GC content of all of the reads contained in the genome file, and can be compared to the mapped and unmapped histograms.

![entropy_for_original_mapped_unmapped](https://github.com/user-attachments/assets/934a4175-8e42-49d9-b78f-a67e7201dbcf)

These strip plots show the entropy of the files of reads. 

### Other output. 

The program also has functionality to BLAST any file of reads vs the Swissprot database (files not included here due to the size of them), as well as CSV files containing all of the data used to create the above plots. 

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

