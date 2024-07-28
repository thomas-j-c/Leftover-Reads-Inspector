import subprocess
from utility.FileHandlingUtils import FileHandler
from utility.DataUtils import DataUtils


class Assembler:
    def __init__(self):
        dataUtils = DataUtils()
        if dataUtils.isPackageInstalled("megahit") and dataUtils.isPackageInstalled("bowtie2"):
            self.installed = True

        else:
            self.installed = False

    def callMegahit(self, inputFiles):
        """
        Calls Megahit to assemble the reads in the input files. These must be fastq files.
        :param inputFiles: List of input files
        """
        command = f"megahit -o ../Data/intermediary/megahit -r {inputFiles}"
        # Call Megahit to assemble the file(s) passed to it. This works because parsed.f is a list of files.
        subprocess.run(command.split(" "))

    def callBowtie(self, genomeFileName, inputFiles, unmappedReadsFileName, mappedReadsFileName):
        """
        Function to call bowtie2 to align to assembly created with callMegahit() function
        :param genomeFileName: Name of the binary files created by Bowtie2
        :param inputFiles: Input files
        :param unmappedReadsFileName: Name of output file of unmapped reads.
        :param mappedReadsFileName: Name of output file of mapped reads.
        """
        createGenome = f"bowtie2-build -q ../Data/intermediary/megahit/final.contigs.fa " \
                       f"../Data/intermediary/bowtie/{genomeFileName} "
        subprocess.run(createGenome.split(" "))  # Creates a custom bowtie genome, to allow it to find the mapped reads.

        findUnmapped = f"bowtie2 --quiet -t --no-hd -x ../Data/intermediary/bowtie/{genomeFileName} -U {inputFiles} " \
                       f"--un ../Data/output/bowtie/{unmappedReadsFileName}.fa " \
                       f"--al ../Data/output/bowtie/{mappedReadsFileName}.fa"
        # Searches the created genome against the input files. Again, supports list of .fasta files input so will work.

        subprocess.run(findUnmapped.split())

    def main(self, inputFiles, genomeFileName="bowtieGenome", unmappedFileName="unmappedReads",
             mappedFileName="mappedReads"):
        """
        The main function in this class. Calls all required fucntions to make an assembly,
        align the set of original reads, and find the mapped and unmapped reads.
        :param inputFiles: Input file(s)
        :param genomeFileName: Name of binary files produced by Bowtie2-build
        :param unmappedFileName: Name of output file of unmapped reads.
        :param mappedFileName: Name of output file of mapped reads.
        """
        if self.installed:
            filesForTools = ", ".join(inputFiles)
            # Creates assembly.
            self.callMegahit(filesForTools)
            # Finds mapped and unmapped reads.
            self.callBowtie(genomeFileName, filesForTools, unmappedFileName, mappedFileName)

        else:
            print("Please check installations of Bowtie2 and Megahit.")
