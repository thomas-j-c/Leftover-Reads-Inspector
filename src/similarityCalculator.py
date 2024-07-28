import difflib
import matplotlib.pyplot as plt

from utility.StatisticsUtils import Statistics
from utility.FileHandlingUtils import FileHandler
from utility.PlotsUtils import PlotsUtils

class FindMatches:
    def __init__(self, n=3, cutoff=0.6):
        """
        Initialises the class
        :param n: The number of closest matches for the program to find.
        :param cutoff: The similarity cut-off.
        """
        self.number = n
        self.cutoff = cutoff

    def findMatches(self, setOne, setTwo, sampleFirst=True, sampleSecond=True, percentage=0.001, seed="Random"):
        """
        Uses difflib.get_close_matches() function to find the n closest matches of a sequence.
        Loops through each of the sequences in the first file, looking for matches in the second.
        :param setOne: The first set, which will be looped through.
        :param setTwo: The second set, to be searched.
        :param sampleFirst: Whether to sample the first set or not. Useful for large sets.
        :param sampleSecond: Whether to sample the second set. Useful for large sets.
        :param percentage: The percentage of the set(s) to sample.
        :param seed: A seed to control the sample.
        :return: A dictionary containing all the sequences from the first file and its closest matches.
        """
        output = dict()
        stats = Statistics()

        if sampleFirst:
            setOne = stats.findRandomSampleDictionary(setOne, percentage, seed)

        if sampleSecond:
            setTwo = stats.findRandomSampleDictionary(setTwo, percentage, seed)

        setTwoSequences = setTwo.keys()
        count = 0
        for sequence in setOne:
            count += 1
            temp = difflib.get_close_matches(sequence, setTwoSequences, n=self.number, cutoff=self.cutoff)
            output.update({sequence: temp})

            print("{:.2f}% compared".format(count / len(setOne) * 100))
        return output

    def findHamming(self, data):
        """
        Finds the hamming distance of two passed sequences. Used in this context to find the distance between a
        sequence and its closest match.
        :param data: The data, output by findMatches()
        :return: The original data, with the distance added to the values.
        """
        stats = Statistics()

        for sequence in data:
            try:
                closest = data.get(sequence)[0]

            except IndexError:
                continue

            data.update({sequence: data.get(sequence) + [stats.findHamming(sequence, closest)]})

        print("Found hamming distances")
        return data

    def main(self, fileOne, fileTwo, sampleFirst=True, sampleComparison=True, percentage=0.001, seed="Random",
             outputCSVName="default_csv_name", histogramOfHammingDists=True, xLabel="Default", yLabel="Default",
             title="Default", histogramFileName="default_histogram"):
        """
        Main function for this class. By default, samples the files inputted to 0.001%, as this was designed for very
        large files.
        :param fileOne: First file for comparison. The sequences in this will be compared against the second file.
        :param fileTwo: Second file for comparison. These will be searched through to find matches.
        :param sampleFirst: Whether to sample the first input file.
        :param sampleComparison: Whether to sample the second input file.
        :param percentage: The percentage to sample the files by.
        :param seed: Seed for the random sample. Default value makes sample random.
        :param outputCSVName: Name of output CSV. It is in format: sequence | first match | second match | third match |
        Hamming distance
        :param histogramOfHammingDists: Whether to produce a histogram of the Hamming distances between closest matches
        :param xLabel: X label for histogram
        :param yLabel: Y label for histogram
        :param title: Title for histogram
        :param histogramFileName: File name of output Histogram
        """
        fileHandler = FileHandler()
        dataOne = fileHandler.getDataFromInputFile(fileOne)
        dataTwo = fileHandler.getDataFromInputFile(fileTwo)

        matches = self.findMatches(dataOne, dataTwo, sampleFirst, sampleComparison, percentage, seed)
        matches = self.findHamming(matches)

        fieldNames = ["Sequence"]
        fieldNames += [f"Match {i}" for i in range(1, self.number + 1)]
        fieldNames += ["Hamming Distance Of Closest"]
        print("Writing to CSV")
        fileHandler.writeToCSVConvertData(fieldNames, matches, outputFile=outputCSVName)

        if histogramOfHammingDists:
            print("Producing Histogram")
            df = fileHandler.convertCSVToDataFrame(f"../Data/output/csv/{outputCSVName}.csv")
            plots = PlotsUtils()

            fig, ax = plt.subplots()
            plots.makeHistogram(df["Hamming Distance Of Closest"], xLabel=xLabel, yLabel=yLabel,
                                title=title, ax=ax)

            plt.savefig(f"../Data/output/plots/{histogramFileName}.png")
