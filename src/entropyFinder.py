import math
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

from utility.PlotsUtils import PlotsUtils
from utility.FileHandlingUtils import FileHandler

class EntropyFinder:
    def findProbabilities(self, sequence):
        """
        Finds the probabilities of each character appearing in the sequence.
        Uses the Counter class to efficiently count the sequence.
        :param sequence: The sequence to search.
        :return: The probabilities of each character.
        """
        probabilities = dict()
        counts = Counter(sequence)  # Counts the occurrences of all sequences in the file.
        length = len(sequence)

        for letter in counts:
            # Finds the probability for each of A,C,G and T. Adds it to the dictionary.
            probabilities.update({letter: (counts.get(letter) / length)})

        return probabilities

    def findEntropy(self, sequence):
        """
        Uses the formula for Shannon Entropy to find the entropy for a sequence.
        :param sequence: A string for the entropy to be found of.
        :return: The entropy value, to 5dp.
        """
        probabilities = self.findProbabilities(sequence)
        entropy = 0

        for char in probabilities:
            # Finds the entropy using the equation in the Shannon Entropy paper
            entropy += probabilities.get(char) * (math.log2(probabilities.get(char)))

        return float("{:.5}".format(-entropy))

    def findAverageMinMaxOfDictionary(self, data):
        """
        A function that finds the average, minimum and maximum entropy values within the dictionary.
        :param data: The dictionary. Sequences must be the keys of it.
        :return: A tuple containing the average, minimum and maximum.
        """
        average = 0
        minimum = 2
        maximum = 0

        for sequence in data:
            # Finds the entropy of each sequence, and compares it to already found entropies.
            entropy = self.findEntropy(sequence)
            average += entropy

            if minimum > entropy:
                minimum = entropy

            if maximum < entropy:
                maximum = entropy

        average /= len(data)

        return average, abs(minimum), maximum

    def findEntropyOfEachSequence(self, data, outputFilename="all_entropies"):
        """
        Finds the entropy of each
        sequence in a dictionary.
        :param data: Filepath of the input dictionary.
        :param outputFilename: The name of the file outputted by the function.
        :return: The data that was written to the CSV, as a dictionary.
        """
        outputData = dict()

        for sequence in data:
            outputData.update({sequence: self.findEntropy(sequence)})

        writer = FileHandler()
        writer.writeToCSVConvertData(["Sequence", "Entropy"], outputData, outputFile=outputFilename)

        return outputData

    def findOutliersAndWriteToCSV(self, data, outputFilename="entropy_outliers"):
        """
        Finds the upper and lower quartile of the dataset. Outputs those reads with entropy that are not in
        that range.
        :param data: The input data, as a dictionary.
        :param outputFilename: The name of the output file.
        :return: Nothing, but produces a CSV file with name outputFilename.
        """
        upperQuartile, lowerQuartile = np.percentile(np.array(list(data.values())), [75, 25])
        print(f"Upper Quartile of Entropy: {upperQuartile}\n"
              f"Lower Quartile: {lowerQuartile}")
        outliers = dict()

        for sequence in data:
            value = data.get(sequence)
            if value > upperQuartile or value < lowerQuartile:
                outliers.update({sequence: value})

        handler = FileHandler()
        handler.writeToCSVConvertData(["Sequence", "Entropy"], outliers, outputFile=f"{outputFilename}.csv")

    def plotDictToLineChart(self, data,  yLabel="Default", xLabel="Default", graphTitle="Default",
                            filename="line_chart"):
        """
        Plots an input dictionary to a line chart. Uses helper function to do this.
        :param data: Input dictionary
        :param yLabel: Y label for the plot
        :param xLabel: X label for the plot
        :param graphTitle: Title for the graph
        :param filename: Name for the output figure
        """
        plots = PlotsUtils()
        plots.makeLineChart(data, yLabel, xLabel, graphTitle, filename)

    def getDataForStripPlot(self, files):
        """
        Gets data for strip plots from CSV files. Uses convertCSVToDataFrame() helper function from file handling class
        to do this.
        :param files: Files of input
        :return: The data from input files.
        """
        dataframes = []
        handler = FileHandler()

        for file in files:
            dataframes.append(handler.convertCSVToDataFrame(file))

        return dataframes

    def createStripPlots(self, files, title="Strip plot", filename="strip_plots"):
        """
        Creates strip plots. If multiple files are inputted, multiple plots will be created on the same figure.
        :param files: Input files
        :param title: Title for the figure
        :param filename: The name of the saved figure
        """
        dataframes = self.getDataForStripPlot(files)
        plots = PlotsUtils()

        if len(files) > 1:  # If there are multiple input files, then multiple plots are created on one figure.
            fig, axs = plt.subplots(len(files), figsize=(15, 12))
            for index, ax in enumerate(axs):
                plots.makeStripPlot(dataframes[index], ax=axs[index])

            fig.supylabel("Strip plots of entropy for input files")

        else:  # If not, then only one plot is created.
            fig, axs = plt.subplots(len(files), figsize=(15, 8))
            plots.makeStripPlot(dataframes[0], ax=axs)
            fig.supylabel("Strip plot for entropy of input file")

        fig.suptitle(title)
        plt.savefig(f"../Data/output/plots/{filename}.png", bbox_inches="tight", dpi=330)

    def main(self, filepath, yLabelForGraph="Entropy", xLabelForGraph="Position in Dictionary",
             title="Default", plotName="line_chart", outlierFileName="entropy_outliers",
             allEntropyFileName="all_entropies", lineChart=False):
        """
        Main function for this class
        :param lineChart: Switch to create a line chart.
        :param filepath: Input filepath
        :param yLabelForGraph: Y label for line graph
        :param xLabelForGraph: X label for line graph
        :param title: Title for line graph
        :param plotName: Name of output line plot
        :param outlierFileName: Name of CSV file of outliers. Determined by the Upper and Lower quartile of the
        entropies in input file
        :param allEntropyFileName: Name of CSV file of all entropies.
        """
        handler = FileHandler()
        data = handler.getDataFromInputFile(filepath)

        average, minEntropy, maxEntropy = self.findAverageMinMaxOfDictionary(data)

        print(f"Average Entropy of Dataset: {average}\n"
              f"Minimum Entropy of Dataset: {minEntropy}\n"
              f"Maximum Entropy of Dataset: {maxEntropy}")

        allEntropies = self.findEntropyOfEachSequence(data, allEntropyFileName)

        if lineChart:
            self.plotDictToLineChart(allEntropies, yLabel=yLabelForGraph, xLabel=xLabelForGraph, graphTitle=title,
                                     filename=plotName)

        self.findOutliersAndWriteToCSV(allEntropies, outputFilename=outlierFileName)
