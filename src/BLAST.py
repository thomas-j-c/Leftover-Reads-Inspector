import subprocess
import traceback

from utility.DataUtils import DataUtils

class BLAST:
    def __init__(self):
        check = DataUtils()
        self.packageInstalled = check.isPackageInstalled("ncbi-blast+")  # Ensures that the package is installed.

    def makeDB(self, filepath, dbName, dbType):
        """
        Function that makes the database, if the user requests this.
        :param filepath: The filepath to the file from which the database should be made. This needs to be a .gz file
        :param dbName: Name of the database to be created.
        :param dbType: Type of database to be created: Protein, Nucleotide
        """
        if self.packageInstalled:
            try:
                command = f"gunzip -c {filepath} | makeblastdb -dbtype {dbType} -title {dbName} -out {dbName}"
                subprocess.run(command.split(" "))

            except FileNotFoundError:
                print(f"File {filepath} not found. Please try again.")

        else:
            print("ncbi-blast+ not installed, or not found. "
                  "On Ubuntu, use 'sudo apt-get install ncbi-blast+' to install")

    def queryFile(self, program, filepath, outputName="blast_output", outputFormat=10, numThreads=6, makeDatabase=False,
                  fileForDB="", dbName="", dbType=""):
        """
        Queries either blastn or blastx on an input file. Allows user to customise the query sent to linux command line
        with parameters.
        :param program: Either "blastn" or "blastx", dependent on which command line option is chosen.
        :param filepath: File path passed to -f on command line. If multiple, then multiple iterations of this are ran.
        :param outputName: Name of the output file. Default blast_output.
        :param outputFormat: Format specifier for Blast. Passed to the function by the -outfmt option on command line
        :param numThreads: Number of threads for BLAST to use to compute. Default 6.
        :param makeDatabase: Bool value for whether a database should be created.
        :param fileForDB: The filepath to the file from which the database should be made. This needs to be a .gz file
        :param dbName: Name of the database to be created.
        :param dbType: Type of database to be created: Protein, Nucleotide
        """
        if self.packageInstalled:
            if makeDatabase:
                print("Beginning creation of BLAST database")
                self.makeDB(fileForDB, dbName, dbType)
                print("created BLAST database")
                print("Beginning BLAST search")

                # Command that calls blast on the custom database.
                # The f String is used to input the user's input to the command.
                command = f"{program} -db {dbName} -query {filepath} -out {outputName} -outfmt {outputFormat}"
                try:
                    subprocess.run(command.split(" "))

                except FileNotFoundError:
                    print(f"Query file {filepath} not found.")

            else:
                print("Beginning BLAST search")

                # Command that calls blast on the SwissProt database.
                command = f"{program} -db ../Data/swissprotdb/uniprot_sprot -query {filepath} -out " \
                          f"../Data/output/blast/{outputName}.csv -outfmt {outputFormat} -num_threads {numThreads}"

                subprocess.run(command.split())

        else:
            print("ncbi-blast+ not installed, or could not be found. "
                  "On Ubuntu, use 'sudo apt-get install ncbi-blast+' to install")

    def customQuery(self, query):
        """
        Allows user to enter a custom query to run on the command line.
        :param query: String query
        :return False if operation could not be performed
        """
        if self.packageInstalled:
            try:
                if "blast" not in query.split()[0]:  # Prevents users entering commands which are not BLAST commands
                    print("Query must begin with blastn or blastx")
                    return False

                subprocess.run(query.split(" "))

            except Exception:
                print("An error occurred with your query. Please try again. Stack trace will be printed.")
                print(traceback.format_exc())
                return False
