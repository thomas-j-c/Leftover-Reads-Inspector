import random

import numpy as np
import statsmodels.stats.proportion as proportionStats
import statsmodels.stats.weightstats as weightStats

from .FileHandlingUtils import *


class Statistics:
    def __init__(self, significanceLevel=0.05):
        self.significance = significanceLevel

    def areProportionsSignificant(self, occurrences, totals):
        """
        Null hypothesis: The proportions of the occurrences are not significantly different
        Alternate Hypothesis: The proportions are significantly different
        :param occurrences: The number of occurrences of the same kmer in both sets
        :param totals: The total number of kmers in both sets.
        :return: Whether the null hypothesis can be rejected, and the p value.
        """

        _, pValue = proportionStats.proportions_ztest(count=np.array(occurrences), nobs=np.array(totals))
        if pValue <= self.significance:
            return "Result is significant", float("{:.6}".format(pValue))

        else:
            return "Result is not significant", float("{:.6}".format(pValue))

    def independentT_Test(self, datasetOne, datasetTwo):
        """
        Null Hypothesis: The means are not significantly different.
        Alt. Hypothesis: The means are significantly different.
        :param datasetOne: The first dataset to test.
        :param datasetTwo: The second dataset to test.
        :return: A tuple containing the result and the p value.
        """
        pValue, _, _ = weightStats.ttest_ind(datasetOne, datasetTwo)

        if pValue <= self.significance:
            # Reject the null hypothesis, means are significantly different.
            return "Result is significant", float("{:.6}".format(pValue))

        else:
            # Fail to reject the null hypothesis, means are not significantly different.
            return "Result is not significant", float("{:.6}".format(pValue))

    def findRandomSampleDataframe(self, data, seed="Random"):
        """
        Finds a sample of the dataframe passed to it. By default, the sample is random, however the value can be
        specified, for testing purposes.
        :param data: The data to find the random sample of
        :param seed: The sample seed.
        :return: The sampled data, which is a dataframe.
        """
        if seed == "Random":
            return data.sample(frac=0.1)
        else:
            return data.sample(frac=0.1, random_state=seed)

    def findRandomSampleDictionary(self, data, percentage, seed="Random"):
        """
        Finds a random sample of a dictionary. Functionality is the same as in the findRandomSampleDataframe() function
        :param data: Input dictionary
        :param percentage: Percentage to sample
        :param seed: The seed to sample with
        :return: The sampled subset of the input data
        """
        if seed == "Random":
            return {k: data[k] for k in random.sample(list(data), round((percentage * len(data))))}

        else:
            random.seed(seed)
            return {k: data[k] for k in random.sample(list(data), round((percentage * len(data))))}

    def getTotalKmersDataFrame(self, dataframe: pd.DataFrame):
        totalSetOne = 0
        totalSetTwo = 0

        for i in range(len(dataframe)):
            # The number of raw occurrences will always be in columns 1 and 2.
            totalSetOne += dataframe.iloc[i, 1]
            totalSetTwo += dataframe.iloc[i, 2]

        return totalSetOne, totalSetTwo

    def testProportions(self, filepath, outputFile="hypothesis_test_results_for_mapped_vs_unmapped.csv", seed="Random"):
        """
        This function takes the data from the file given when class initialised and converts to dataframe.
        If there are more than 100 samples, then it takes a sample of 10% of those. This can be controlled with the
        seed parameter.
        It loops through this dataframe, and performs a proportion test using areProportionsSignificant() in this class.
        Then writes the returned data to a new csv file, name specified with the outputFile parameter.
        :param outputFile: The name of the output file the user would like.
        :param seed: The pseudo-random seed the user would like.
        :return: Nothing. Produces a file in the Data/output directory.
        """
        handler = FileHandler()

        dataframe = handler.convertCSVToDataFrame(filepath)
        if len(dataframe.index) > 100:
            dataframe = self.findRandomSampleDataframe(dataframe, seed)

        totalSetOne, totalSetTwo = self.getTotalKmersDataFrame(dataframe)

        dataToWrite = []

        for i in range(len(dataframe)):
            data = dict()
            # The arrays have to be converted to numpy arrays to ensure that they can be used by the functions.
            result, pvalue = self.areProportionsSignificant(np.array([dataframe.iloc[i, 1], dataframe.iloc[i, 2]],),
                                                            np.array([totalSetOne, totalSetTwo]))

            # Converts data to format used by the write to CSV function.
            data.update({"Kmer": dataframe.iloc[i, 0], "Verdict": result, "P-Value": pvalue})
            dataToWrite.append(data)

        handler.writeToCSVNoConvert(["Kmer", "Verdict", "P-Value"], dataToWrite,
                                    outputFile=outputFile)

    def extractSignificantKmers(self, data="hypothesis_test_results_for_mapped_vs_unmapped.csv",
                                outputFile="significant_kmers.csv"):
        """
        Extracts the kmers that have significant difference in their proportions in the sets.
        :param data: The csv file created by testOnMultipleInstances function above.
        :param outputFile: The name of the output file the user would like.
        :return: Nothing, but creates file in Data/output directory with the kmer and its p-value.
        """
        handler = FileHandler()
        dataframe = handler.convertCSVToDataFrame(f"../../Data/output/csv{data}")
        dataToWrite = []

        for i in range(0, len(dataframe)):
            if dataframe.iloc[i, 1] == "Result is significant":
                data = dict()
                data.update({"Kmer": dataframe.iloc[i, 0], "P-Value": dataframe.iloc[i, 2]})

                dataToWrite.append(data)

        handler = FileHandler()
        handler.writeToCSVNoConvert(["Kmer", "P-Value"], dataToWrite, outputFile=outputFile)

    def findHamming(self, sequenceOne, sequenceTwo):
        """
        Finds the Hamming distance of two sequences. Sequences must be of the same length.
        :param sequenceOne: The first sequence.
        :param sequenceTwo: The second sequence.
        :return: The distance.
        """
        distance = 0

        for i, _ in enumerate(sequenceOne):
            try:
                # If the characters at the same position are the same, then they have a distacne of 0.
                if sequenceOne[i] == sequenceTwo[i]:
                    continue

                else:
                    # If not, then the distance is incremented
                    distance += 1

            except IndexError:
                # If the end of one of the sequences is found, then an Index Error is raised. This provides robustness
                # for the program.
                return distance

        return distance
