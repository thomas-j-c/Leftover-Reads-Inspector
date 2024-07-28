import subprocess
from .FileHandlingUtils import FileHandler
import pandas as pd

class DataUtils:
    # A function that checks whether a given package is installed
    def isPackageInstalled(self, package):
        """
        Function that checks if a given package is installed
        :param package: String name of package
        :return: False if package is not installed, otherwise True
        """
        try:
            # Makes the output of this call silent.
            null = open("/dev/null", "w")

            # Attempts to execute the package passed to it.
            subprocess.Popen(f"{package}", stdout=null, stderr=null)
            null.close()

            # If no error is raised, then the package is installed
            return True

        except OSError:
            print(f"{package} not installed.")
            return False

    # Returns the difference in the counts of Kmers, along with the number of occurrences in both sets.

    def compareRawKmerCountsFromJellyfish(self, dataOne, dataTwo):
        """
        Finds the difference of the raw counts of k-mers in the files
        :param dataOne: The first set of data
        :param dataTwo: The second set of data
        :return: A dictionary containing the k-mer, its number of occurrences in both sets, and the difference
        between them.
        """
        handler = FileHandler()
        setOne = handler.getDataFAFile(dataOne)
        setTwo = handler.getDataFAFile(dataTwo)

        if len(setOne) > len(setTwo):
            intersection = handler.intersection(setOne, setTwo)  # If the first set is larger than the second
        else:
            intersection = handler.intersection(setTwo, setOne)  # If the second set is larger than the first

        for kmer in intersection:
            occurrences = intersection.get(kmer)
            # Finds the difference between occurrences of k-mers, and appends it to the list.
            occurrences.append(abs(int(occurrences[0] - occurrences[1])))

            # Updates the dictionary to contain the new list.
            intersection.update({kmer: occurrences})

        return intersection

    def findTotalKmersFromJellyFishFA(self, data):
        """
        A function that finds the total number of k-mers in a file outputted by Jellyfish
        :param data: Dictionary containing k-mer: occurrences
        :return: int total number of k-mers
        """
        total = 0
        for i, kmer in enumerate(data):
            total += int(data.get(kmer)[0][1:])

        return total

    def normalise(self, occurrences, total):
        """
        Finds the normalised value of the occurrences of a k-mer.
        :param occurrences: Number of occurrences of individual k-mer.
        :param total: Total number of k-mers
        :return: A float value of te normalised value. 8dp
        """
        return float(("{:.8f}".format(occurrences / total)))

    def findDifferenceInNormalised(self, occurrences_setOne, total_setOne, occurrences_setTwo, total_set2):
        """
        Finds the difference between two normalised values.
        :param occurrences_setOne: Number of occurrences of a k-mer in set one
        :param total_setOne: Total number of k-mers in set one
        :param occurrences_setTwo: Number of occurrences of a k-mer in set two
        :param total_set2: Total number of k-mers in set two
        :return: Tuple containing the two normalised values, and the difference between them.
        """
        return (self.normalise(occurrences_setOne, total_setOne), self.normalise(occurrences_setTwo, total_set2),
                abs(self.normalise(occurrences_setOne, total_setOne) - self.normalise(occurrences_setTwo, total_set2)))

    def compareNormalisedKmerCounts(self, dataOne, dataTwo):
        """
        Updates the set of intersecting kmers with the normalised numbers, and the difference.

        This is required for the implementation of the writeToCSV function, as it must have the same number of
        column labels and instances of data.
        :param dataOne: The first input file.
        :param dataTwo: The second input file.
        """
        handler = FileHandler()
        setOne = handler.getDataFAFile(dataOne)
        setTwo = handler.getDataFAFile(dataTwo)

        if len(setOne) > len(setTwo):
            intersection = handler.intersection(setOne, setTwo)
        else:
            intersection = handler.intersection(setTwo, setOne)

        totalKmersSetOne = self.findTotalKmersFromJellyFishFA(setOne)
        totalKmersSetTwo = self.findTotalKmersFromJellyFishFA(setTwo)

        # The below functionality follows the same basic format as compareRawKmerCountsFromJellyfish.
        for kmer in intersection:
            occurrences = intersection.get(kmer)

            setOneNormalised, setTwoNormalised, difference = self.findDifferenceInNormalised(occurrences[0],
                                                                                             totalKmersSetOne,
                                                                                             occurrences[1],
                                                                                             totalKmersSetTwo)

            normalisedNums = [item for item in [setOneNormalised, setTwoNormalised, difference]]

            intersection.update({kmer: normalisedNums})

        return intersection

    def getRawAndNormalised(self, dataOne, dataTwo):
        """
        Uses helper functions to find all values for all k-mers.
        :param dataOne: The first input file.
        :param dataTwo: The second input file.
        :return:
        """
        rawData = self.compareRawKmerCountsFromJellyfish(dataOne, dataTwo)
        normalisedData = self.compareNormalisedKmerCounts(dataOne, dataTwo)

        allData = dict()
        for kmer in rawData:
            allData.update({kmer: rawData.get(kmer) + normalisedData.get(kmer)})

        return allData

    def convertToMPLFormat(self, data, dataHeader):
        """
        A function to convert different types of data into a forma that can be used by Matplotlib. Not widely used in
        the tool, however kept with expansion or future repurposing in mind.
        :param data: Input data to convert
        :param dataHeader: Part of the dat to extract.
        :return: The extracted data.
        """
        outputData = []

        if type(data) == dict:  # If the data passed is a dictionary.
            for row in [data]:
                outputData.append(row.get(dataHeader))

        elif type(data) == list:  # If the data is presented as a list of dictionaries.
            for row in data:
                outputData.append(row.get(dataHeader))

        elif type(data) == str:
            try:
                dataFrame = pd.read_csv(data)

                for i in range(0, len(dataFrame.index)):
                    outputData.append(dataFrame.iloc[i][dataHeader])

            except FileNotFoundError:  # If the data provided is not a csv file, then return nothing.
                print("File not found. Please try again")
                return

            except KeyError:
                print(f"Data header {dataHeader} passed is not contained in the passed file")
                return

        else:
            print("Data must be of type dict, list of dicts, or .csv filepath")
            return

        return outputData
