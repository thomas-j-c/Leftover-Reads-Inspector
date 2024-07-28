import csv
import subprocess
import sys
from itertools import islice
import pandas as pd


class FileHandler:

    """
    All getData... functions take the filepath as the filename parameter,
    and all return dictionaries returned by the getData... functions are in form: {sequence: [metadata]}
    """

    def getDataMAPFile(self, filename):
        maxValue = sys.maxsize
        reads = dict()

        with open(filename, "r") as file:
            csvReader = csv.reader(file, delimiter="\t")  # These files are tab delimited.
            for row in csvReader:
                try:
                    # Attempts to read the data using the max value. This allows large data to be read.
                    csv.field_size_limit(maxValue)
                except OverflowError:
                    # If the program encounters an OverflowError, then the max value is made smaller.
                    maxValue = int(maxValue / 10)

                metadata = [row[0], row[1], row[2], row[3], row[5]]
                reads.update({row[4]: metadata})

        return reads

    def getDataFASTQfile(self, filename):
        reads = dict()

        with open(filename, "r") as file:
            for row in file:
                if row[0] == '@':  # If the line begins with '@', we know it is the start of a new record.
                    lines = islice(file, 3)  # Gets the next 3 lines from the file.
                    metadata = [row.strip()]

                    for line in lines:
                        metadata.append(line.strip())

                    reads.update({metadata[1]: metadata})  # metadata[1] is the sequence.
        return reads

    def getDatasamFile(self, filename):
        data = dict()

        with open(filename, "r") as file:
            csvReader = csv.reader(file, delimiter="\t")  # These files are tab delimited.
            for row in csvReader:
                metaData = [row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[10]]
                data.update({row[9]: metaData})

        return data

    # Finds the sequences that are in the larger set but not the smaller set
    def getUnmappedReads(self, originalReadsFile, assembledReadsFile):
        """
        Function that finds the reads that are in a larger file and the smaller one. Purpose to find unmapped
        reads, hence the name.
        :param originalReadsFile: Larger input file.
        :param assembledReadsFile: Smaller input file.
        :return: A dictionary of the reads in larger file and are in the smaller file.
        """
        originalReads = self.getDataFromInputFile(originalReadsFile)
        mappedReads = self.getDataFromInputFile(assembledReadsFile)
        unmappedReads = dict()

        for read in originalReads:
            if read in mappedReads:
                continue
            else:
                unmappedReads.update({read: originalReads.get(read)})

        return unmappedReads

    def intersection(self, largerSet, smallerSet):
        """
        Function that finds the reads that are in a larger file, but not the smaller one. Purpose to find unmapped
        reads, hence the name.
        :param largerSet: Larger input file.
        :param smallerSet: Smaller input file.
        :return: A dictionary of the reads in larger file that are not in the smaller file.
        """
        leftover = dict()
        for read in largerSet:  # Looking through each item in the larger set of sequences
            if read not in smallerSet:  # If that read is not in the smaller set, ignore it
                continue

            else:
                print(largerSet)
                print(smallerSet)
                leftover.update({read: [int(largerSet.get(read)[0][1:]), int(smallerSet.get(read)[0][1:])]})

        return leftover

    def writeToFASTQ(self, data, outputName="output"):
        """
        Writes a dictionary of data to a file in the fastq format.
        :param data: Dictionary of data
        :param outputName: Name of output fastq file.
        """
        with open(f"../Data/output/fastq/{outputName}.fastq", "w+") as file:
            for read in data:
                readData = data.get(read)

                identification = readData[0].strip()
                positiveNegative = readData[2]
                quality = readData[3]
                sequence = read.strip()

            file.write(f"{identification}\n"
                       f"{sequence}\n"
                       f"{positiveNegative}\n"
                       f"{quality}\n")

    def convertSAMToFASTA(self, filepath, outputName="output_fasta"):
        """
        Uses the subprocess library and uses samtools to convert sam files to fasta. Made redundant by the --al and --un
        options for Bowtie2, however kept again for future proofing.
        :param filepath: Filepath of a sam file
        :param outputName: Name of output file
        """
        command = f"samtools -f 4 {filepath} > ../Data/output/fasta/{outputName}"
        subprocess.run(command.split())

    def getDataFAFile(self, fileName):
        data = dict()
        sequences = []
        ids = []
        currentContig = ""

        with open(fileName, "r") as file:
            for line in file:
                # If the line begins with a '>', we know until the next line beginning with this, it is a sequence
                if line[0] == ">":
                    ids.append(line.strip())
                    sequences.append(currentContig)  # Append the previous sequence to the list.
                    currentContig = ""

                else:
                    currentContig += line.strip()

        sequences.append(currentContig)

        sequences = [val for val in sequences if val != ""]

        for contig in range(0, len(ids)):
            data.update({sequences[contig]: [ids[contig]]})

        return data

    def convertCSVToDataFrame(self, data="None"):
        """
        Converts a CSV file to a pandas dataframe
        :param data: Filepath of CSV
        :return: The pandas dataframe.
        """
        return pd.read_csv(data)

    def getDataFromInputFile(self, filepath):
        """
        Returns the data from an input file, if it is of one of the accepted input types.
        :param filepath: The file to get data from.
        :return: The data, or an error message.
        """
        try:
            if filepath[-2:] == "fa" or filepath[-5:] == "fasta":
                return self.getDataFAFile(filepath)

            elif filepath[-5:] == "fastq" or filepath[-2:] == "fq":
                return self.getDataFASTQfile(filepath)

            elif filepath[-3:] == "map":
                return self.getDataMAPFile(filepath)

            elif filepath[-3:] == "csv":
                return self.convertCSVToDataFrame(filepath).to_dict()

            elif filepath[-3:] == "sam":
                return self.getDatasamFile(filepath)

            else:
                print("File passed must be fasta, fastq, map, sam, or csv")
                return None

        except FileNotFoundError:
            print(f"File path {filepath} not recognised as .fasta (fa), .fastq (fq), .map, .sam, or .csv/ "
                  f"no such file exists. Please try again.")
            return None

        except TypeError:
            print("None type error. Ensure the file passed is correct")
            return None

    def convertToCSVWriteable(self, fields, inputData):
        """
        Converts an input dictionary to form accepted by the csv.dictwriter() function, which is a 2d list of
        dictionaries, containing the column header as the key, and the value for that row as its value.
        :param fields: The column headers for the CSV.
        :param inputData: The dictionary to convert
        :return: A list of dictionaries.
        """
        if type(inputData) == str:
            data = self.getDataFAFile(inputData)  # If the input data is a .fa file

        elif type(inputData) == dict:
            data = inputData  # If the input data is a dictionary.

        else:
            return "Err: Input data must either be a filepath or dictionary."

        listOfDicts = []
        newFields = fields[1:]

        for row in data:
            dictOfData = dict()

            dictOfData.update({f"{fields[0]}": row})

            if type(data.get(row)) == list:
                for i in range(0, len(data.get(row))):
                    # Loop through the length of the list, starting from 1.
                    # Do the same fields[i+1]: but instead of value do list[i]
                    # Where list has come from data.get(row). So for i in range(1, len(data.get[row]):
                    dictValues = data.get(row)

                    if str(dictValues[i])[:1] == ">":
                        dictOfData.update({f"{newFields[i]}": f"{dictValues[i][1:]}"})

                    else:
                        dictOfData.update({f"{newFields[i]}": f"{dictValues[i]}"})

            else:
                dictOfData.update({f"{fields[1]}": data.get(row)})

            listOfDicts.append(dictOfData)

        return listOfDicts

    def writeToCSVConvertData(self, fieldNames, inputData, outputFile="output_csv.csv"):
        """
        Use if input data is not already in list of dicts form.
        :param fieldNames: The names of the columns for the CSV file.
        :param inputData: The input data. Can be a filepath or regular dictionary.
        :param outputFile: The name of the file to write to.
        :return: Nothing. Produces file in Data/output/ directory.
        """

        data = self.convertToCSVWriteable(fieldNames, inputData)
        with open(f"../Data/output/csv/{outputFile}.csv", "w") as csvFile:
            writer = csv.DictWriter(csvFile, fieldnames=fieldNames)
            # Writes the field headers.
            writer.writeheader()
            # Writes the dictionary of data to the csv file.
            # Data must be in the form of a list of dictionaries, with the keys being the same as the fields,
            # and values what should be written to those fields.
            writer.writerows(data)

    def writeToCSVNoConvert(self, fieldNames, inputData, outputFile="output_csv.csv"):
        """
        Use if data is in form of list of dictionaries.
        :param fieldNames: The names of the columns for the CSV file.
        :param inputData: The input data. Can be a filepath or regular dictionary.
        :param outputFile: The name of the file to write to.
        :return: Nothing. Produces file in Data/output/ directory.
        """
        with open(f"../Data/output/csv/{outputFile}.csv", "w") as csvFile:
            writer = csv.DictWriter(csvFile, fieldnames=fieldNames)
            writer.writeheader()
            writer.writerows(inputData)
