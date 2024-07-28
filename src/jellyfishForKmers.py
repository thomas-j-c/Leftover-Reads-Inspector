import subprocess
from utility.FileHandlingUtils import FileHandler
from utility.DataUtils import DataUtils


class JellyFish:
    def __init__(self):
        dataUtils = DataUtils()
        self.packageInstalled = dataUtils.isPackageInstalled("jellyfish")

    def runJellyfish(self, inputFile, k=7, jellyfishOutputFile=f"mer_counts", csvName="output_csv"):
        """
        Function that calls Jellyfish on input file. Jellyfish creates a .fa file, and a CSV file is made from that,
        in form of k-mer | Number of occurrences
        :param inputFile: Input file of data.
        :param k: Size of k for Jellyfish to search with
        :param jellyfishOutputFile: Name of output file of Jellyfish
        :param csvName: Name of output CSV created from the data found by Jellyfish.
        """
        # packageInstalled is a boolean value containing the output of isPackageInstalled("Jellyfish").
        if self.packageInstalled:
            cmd = "jellyfish"

            # Command to find the k-mers.
            subprocess.run([cmd, "count", f"-m {k}", "-s 100M", "-t 10", inputFile,
                            "-o../Data/intermediary/mer_counts.jf"], stdout=subprocess.PIPE)

            # Converts the binary file to a fasta file.
            subprocess.run([cmd, "dump", "../Data/intermediary/mer_counts.jf",
                            f"-o../Data/output/fa/{jellyfishOutputFile}"], stdout=subprocess.PIPE)

            # Takes the output fasta file and converts it to a CSV file, to allow further analysis.
            writer = FileHandler()
            writer.writeToCSVConvertData(["kmer", "occurrences"], inputData=f"../Data/output/fa/{jellyfishOutputFile}",
                                         outputFile=f"{csvName}")

        else:
            print("Package Jellyfish is not installed. Please install it using\n"
                  "sudo apt install jellyfish\n"
                  "on Ubuntu linux, before attempting command again.")
