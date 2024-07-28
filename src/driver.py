import argparse
import traceback

from jellyfishForKmers import JellyFish
from entropyFinder import EntropyFinder
from GCContent import GCCalculator
from similarityCalculator import FindMatches
from utility.DataUtils import *
from utility.StatisticsUtils import Statistics as Stats
from assembleAndFindUnmapped import Assembler
from BLAST import BLAST


def searchFilesWithJellyfish():
    finder = JellyFish()
    count = 0
    for file in parsed.f:
        print(file)
        try:
            if parsed.k:
                finder.runJellyfish(file, parsed.k, jellyfishOutputFile=f"{parsed.k}mer_counts{count}.fa",
                                    csvName=f"{parsed.k}mer_output{count}")

            else:
                finder.runJellyfish(file, jellyfishOutputFile=f"7mer_counts{count}.fa",
                                    csvName=f"7mer_output{count}")
        except:
            print("An error has occurred with Jellyfish. The accepted filetypes are fasta and fastq.")
            print(traceback.format_exc())

        count += 1


def proportionsStatsTest(filePath, outputfile="hypothesis_test_results", significanceLevel=0.05):
    """
    Runs a hypothesis test on the proportions of the k-mers in two files.
    :param filePath:
    :param outputfile: Name of the output file.
    :param significanceLevel: The level of significance to test on.
    """
    if significanceLevel != 0.05:
        statsFinder = Stats(significanceLevel=significanceLevel)

    else:
        statsFinder = Stats()

    try:
        statsFinder.testProportions(filepath=f"../Data/output/csv/{filePath}.csv", outputFile=outputfile)

    except FileNotFoundError:
        print(f"File ../Data/output/csv/{filePath} not found.")


def findEntropy(yLabelForGraph="Entropy", xLabelForGraph="Position in Dictionary",
                title="Default", plotName="line_chart", csvFileName="entropy_outliers", stripPlot=True,
                lineChart=False):
    """
    Provides an interface between command line and the entropy finding class. Produces a strip plot and
    a csv file for further inspection.
    :param yLabelForGraph: Y label for the graph.
    :param xLabelForGraph: X label for the graph.
    :param title: Title for the graph.
    :param plotName: Name of the output plot file.
    :param csvFileName: Name of the output CSV file.
    :param stripPlot: Switch to make the program produce a strip plot.
    :param lineChart: Switch to make the program produce a line chart.
    """
    if parsed.yLabel:
        yLabelForGraph = parsed.yLabel

    if parsed.xLabel:
        xLabelForGraph = parsed.xLabel

    if parsed.title:
        title = parsed.title

    if parsed.output_line_plot_name:
        plotName = parsed.output_line_plot_name

    if parsed.csv_filename:
        csvFileName = parsed.csv_filename

    if parsed.strip_plot:
        stripPlot = parsed.strip_plot

    if parsed.line_chart:
        lineChart = parsed.line_chart

    calculator = EntropyFinder()

    count = 0
    for file in parsed.f:
        calculator.main(file, yLabelForGraph, xLabelForGraph, title,
                        plotName=f"{plotName}_for_file{count}",
                        outlierFileName=f"{csvFileName}_for_file_{count}",
                        allEntropyFileName=f"all_entropies_for_file_{count}",
                        lineChart=lineChart)

        count += 1

    if stripPlot:
        files = [f"../Data/output/csv/all_entropies_for_file_{count}.csv" for count in range(0, len(parsed.f))]
        if parsed.output_strip_plot_name:
            calculator.createStripPlots(files, title=title, filename=parsed.output_strip_plot_name)

        else:
            calculator.createStripPlots(files, title=title, filename=parsed.output_strip_plot_name)


def findGC(yLabelForGraph="Entropy", xLabelForGraph="Position in Dictionary",
           title="Default", histName="output_histogram", barName="output_bar_chart",
           csvFileName="gc_counts_total_and_bases",
           histogram=True, barChart=True, binSize=10):
    """
    Provides the interface between the command line and the GC counting class.
    :param yLabelForGraph: The y labels for the output graphs.
    :param xLabelForGraph: The x labels for the output graphs.
    :param title: Titles for graphs.
    :param histName: The name of the output histogram file.
    :param barName: The name of the output bar chart.
    :param csvFileName: Name of the CSV file outputted by the class.
    :param histogram: Whether to produce a histogram.
    :param barChart: Whether to produce a bar chart.
    :param binSize: The bin size for the histogram. Default 10.
    """
    if parsed.yLabel:
        if type(parsed.yLabel) == list:
            yLabelForGraph = " ".join(parsed.yLabel)

        else:
            yLabelForGraph = parsed.yLabel

    if parsed.xLabel:
        if type(parsed.xLabel) == list:
            xLabelForGraph = " ".join(parsed.xLabel)

        else:
            xLabelForGraph = parsed.xLabel

    if parsed.title:
        if type(parsed.title) == list:
            title = " ".join(parsed.title)

        else:
            title = parsed.title

    if parsed.csv_filename:
        csvFileName = parsed.csv_filename

    if parsed.bar_chart:
        barChart = parsed.bar_chart

    if parsed.histogram:
        histogram = parsed.histogram

    if parsed.bin_size:
        binSize = parsed.bin_size

    if parsed.hist_name:
        histName = parsed.hist_name

    if parsed.bar_name:
        barName = parsed.bar_name

    calculator = GCCalculator()

    try:
        calculator.main(parsed.f, barChart, histogram, csvFileName=f"{csvFileName}", binSize=binSize,
                        yLabel=yLabelForGraph, xLabel=xLabelForGraph, title=title, histName=histName,
                        barName=barName)

    except FileNotFoundError:
        print(f"One file in files parsed not found for gc content calculation. Please try again.")
        print(traceback.format_exc())


def findSimilar(number=3, cutoff=0.6, sampleFirst=True, sampleSecond=True, percentage=0.001, seed="Random",
                outputCSVName="default_csv_name", produceHistogram=True, outputHistogramName="default_histogram",
                xLabel="Default", yLabel="Default", title="Default"):
    """
    Function that provides an interface between the command line and the similarity calculator class.
    :param number: The number of the closest matches to find. Default 3.
    :param cutoff: The value of the cutoff, for which after no more matches are found. Default 0.6.
    :param sampleFirst: Whether the program should sample the first set or not.
    :param sampleSecond: As above.
    :param percentage: The percentage of the sets to sample.
    :param seed: Sets the seed for the sample. Default random.
    :param outputCSVName: The name for the output CSV.
    :param produceHistogram: Whether to produce a histogram of the hamming distances between a sequence and its closest
    match.
    :param outputHistogramName: Name for the output histogram, if requested.
    :param xLabel: X Label for histogram.
    :param yLabel: Y Label for histogram.
    :param title: Title for histogram.
    :return:
    """

    if parsed.number_of_closest:
        number = parsed.number_of_closest

    if parsed.cutoff:
        cutoff = int(parsed.cutoff)

    if parsed.sample_first:
        sampleFirst = parsed.sample_first

    if parsed.sample_second:
        sampleSecond = parsed.sample_second

    if parsed.sample_percentage:
        percentage = float(parsed.sample_percentage)

    if parsed.sample_seed:
        seed = parsed.sample_seed

    if parsed.csv_filename:
        outputCSVName = parsed.csv_filename

    if parsed.histogram:
        produceHistogram = parsed.histogram

    if parsed.xLabel:
        xLabel = parsed.xLabel

    if parsed.yLabel:
        yLabel = parsed.yLabel

    if parsed.title:
        title = parsed.title

    if parsed.hist_name:
        outputHistogramName = parsed.hist_name

    files = parsed.f
    if len(files) != 2:
        print("Error: There must be 2 for comparison")

    matches = FindMatches(number, cutoff)

    try:
        matches.main(files[0], files[1], sampleFirst, sampleSecond, percentage, seed, outputCSVName,
                     histogramOfHammingDists=produceHistogram, xLabel=xLabel, yLabel=yLabel, title=title,
                     histogramFileName=outputHistogramName)

    except FileNotFoundError:
        print(f"File {files[0]} or {files[1]} not found. Please try again.")


def getUnmappedReads(outputFilename="unmapped_reads"):
    """
    Function to provide the functionality to get the unmapped reads from an assembly and a set of reads.
    :param outputFilename: The name of the file to output.
    :return:
    """
    if parsed.fastq_name:
        outputFilename = parsed.fastq_name

    files = parsed.f

    try:
        handler = FileHandler()
        handler.writeToFASTQ(handler.getUnmappedReads(files[0], files[1]), outputName=outputFilename)
        print(f"Unmapped reads found and written to '../Data/output/fastq/{outputFilename}.fastq'")

    except FileNotFoundError:
        print(f"File {files[0]} or {files[1]} not recognised as a valid filepath. Please try again.")


def loadCommandFile(filepath):
    try:
        # Automatically closes the file after exiting the indent.
        with open(filepath, "r") as commandFile:
            # Reads in all the lines from the file
            commands = commandFile.readlines()

        # For each line in the file.
        for commandLine in commands:
            variable = ""
            values = []

            # Splits the command into its elements, and loops through it.
            for index, command in enumerate(commandLine.split(" ")):
                # If the command begins with a "-" character, then it is a flag or option for the program.
                if command[0] == "-":
                    if index != 0:
                        if len(values) > 1 or variable == "f":
                            # Sets the attribute given by variable to the values, if the length of values is greater
                            # than one
                            parsed.__setattr__(variable, values)
                        else:
                            # Sets the value of the attribute variable to the first, and only, item in the values list.
                            parsed.__setattr__(variable, values[0])
                        values = []

                    # Variable is set to the string beginning with "-", but excluding the "-".
                    variable = command[1:].strip()
                else:
                    values.append(command.strip())

            if len(values) > 1:
                parsed.__setattr__(variable, values)
            else:
                parsed.__setattr__(variable, values[0])

            run()

    except FileNotFoundError:
        print("File Path Not Found. Please try again.")


def assembleAndFindUnmapped(assemblyFileName="bowtieGenome", mappedReadsFileName="mappedReads",
                            unmappedReadsFileName="unmappedReads"):
    # TODO: Allow user to add extra options to bowtie2, build and megahit
    if parsed.assembly_file_name:
        assemblyFileName = parsed.assembly_file_name

    if parsed.mapped_reads_file_name:
        mappedReadsFileName = parsed.mapped_reads_file_name

    if parsed.unmapped_reads_file_name:
        unmappedReadsFileName = parsed.unmapped_reads_file_name

    assembler = Assembler()
    assembler.main(parsed.f, assemblyFileName, unmappedReadsFileName, mappedReadsFileName)


def queryBLAST(program, blastOutputName="blast_output.csv", outfmt=10, makeDatabase=False, fileForDatabase="",
               databaseName="", databaseType="", threads=6):
    if parsed.blast_output_name:
        blastOutputName = parsed.blast_output_name

    if parsed.outfmt:
        outfmt = parsed.outmt

    if parsed.make_database:
        makeDatabase = parsed.make_database

    if parsed.file_for_custom_db:
        fileForDatabase = parsed.file_for_custom_db

    if parsed.db_name:
        databaseName = parsed.db_name

    if parsed.db_type:
        databaseType = parsed.db_type

    if parsed.threads:
        threads = parsed.threads

    blast = BLAST()
    if len(parsed.f) > 1:
        for file in parsed.f:
            blast.queryFile(program, file, blastOutputName, outfmt, makeDatabase, fileForDatabase, databaseName,
                            databaseType, threads)

    else:
        blast.queryFile(program, parsed.f[0], blastOutputName, outfmt, makeDatabase, fileForDatabase, databaseName,
                        databaseType, threads)


def run():
    if parsed.assemble_and_find_unmapped:
        print("Assembling")
        assembleAndFindUnmapped()
        parsed.assemble_and_find_unmapped = False

    elif parsed.blastn:
        print("Blasting with blastn")
        queryBLAST("blastn")

    elif parsed.blastx:
        print("Blasting with blastx")
        queryBLAST("blastx")

    elif parsed.get_unmapped:
        print("Getting Unmapped")
        getUnmappedReads()
        parsed.get_unmapped = False

    elif parsed.kmers:  # If user wants just the kmers to be found using jellyfish
        print("Finding kmers")
        searchFilesWithJellyfish()
        parsed.kmers = False

    elif parsed.compare_kmers:  # If the user wants the kmers found by jellyfish to be compared
        print("Comparing kmers")
        searchFilesWithJellyfish()

        handler = FileHandler()
        if parsed.k:
            if parsed.csv_filename:
                handler.writeToCSVConvertData(
                    ["kmer", "Num in set one", "Num in set Two", "Raw Difference", "Normalised in set One",
                     "Normalised set Two", "Norm Difference"], dataUtils.getRawAndNormalised(
                        f"../Data/output/fa/{parsed.k}mer_counts0.fa", f"../Data/output/fa/{parsed.k}mer_counts1.fa"),
                    outputFile=parsed.csv_filename)
            else:
                handler.writeToCSVConvertData(
                    ["kmer", "Num in set one", "Num in set Two", "Raw Difference", "Normalised in set One",
                     "Normalised set Two", "Norm Difference"], dataUtils.getRawAndNormalised(
                        f"../Data/output/fa/{parsed.k}mer_counts0.fa", f"../Data/output/fa/{parsed.k}mer_counts1.fa"),
                    outputFile=f"{parsed.k}mer_counts_raw_and_normalised")

        else:
            if parsed.csv_filename:
                handler.writeToCSVConvertData(
                    ["kmer", "Num in set one", "Num in set Two", "Raw Difference", "Normalised in set One",
                     "Normalised set Two", "Norm Difference"], dataUtils.getRawAndNormalised(
                        f"../Data/output/fa/7mer_counts0.fa", f"../Data/output/fa/7mer_counts1.fa"),
                    outputFile=parsed.csv_filename)

            else:
                handler.writeToCSVConvertData(
                    ["kmer", "Num in set one", "Num in set Two", "Raw Difference", "Normalised in set One",
                     "Normalised set Two", "Norm Difference"], dataUtils.getRawAndNormalised(
                        f"../Data/output/fa/7mer_counts0.fa", f"../Data/output/fa/7mer_counts1.fa"),
                    outputFile="7mer_counts_raw_and_normalised")

        parsed.compare_kmers = False

    elif parsed.stats_for_kmers:
        print("Finding stats for kmers")
        searchFilesWithJellyfish()

        handler = FileHandler()
        if parsed.k:
            if parsed.csv_filename:
                handler.writeToCSVConvertData(
                    ["kmer", "Num in set one", "Num in set Two", "Raw Difference", "Normalised in set One",
                     "Normalised set Two", "Norm Difference"], dataUtils.getRawAndNormalised(
                        f"../Data/output/fa/{parsed.k}mer_counts0.fa", f"../Data/output/fa/{parsed.k}mer_counts1.fa"),
                    outputFile=f"{parsed.csv_filename}")

                if parsed.s_l:
                    proportionsStatsTest(f"{parsed.csv_filename}", significanceLevel=parsed.s_l)

                else:
                    proportionsStatsTest(f"{parsed.csv_filename}")

            else:
                handler.writeToCSVConvertData(
                    ["kmer", "Num in set one", "Num in set Two", "Raw Difference", "Normalised in set One",
                     "Normalised set Two", "Norm Difference"], dataUtils.getRawAndNormalised(
                        f"../Data/output/fa/{parsed.k}mer_counts0.fa", f"../Data/output/fa/{parsed.k}mer_counts1.fa"),
                    outputFile=f"{parsed.k}mer_counts_raw_and_normalised")

                if parsed.s_l:
                    proportionsStatsTest(f"{parsed.k}mer_counts_raw_and_normalised", significanceLevel=parsed.s_l)

                else:
                    proportionsStatsTest(f"{parsed.k}mer_counts_raw_and_normalised")

        else:
            if parsed.csv_filename:
                handler.writeToCSVConvertData(
                    ["kmer", "Num in set one", "Num in set Two", "Raw Difference", "Normalised in set One",
                     "Normalised set Two", "Norm Difference"], dataUtils.getRawAndNormalised(
                        f"../Data/output/fa/7mer_counts0.fa", f"../Data/output/fa/7mer_counts1.fa"),
                    outputFile=f"{parsed.csv_filename}")

                if parsed.s_l:
                    proportionsStatsTest(f"{parsed.csv_filename}", significanceLevel=parsed.s_l)

                else:
                    proportionsStatsTest(f"{parsed.csv_filename}")

            else:
                handler.writeToCSVConvertData(
                    ["kmer", "Num in set one", "Num in set Two", "Raw Difference", "Normalised in set One",
                     "Normalised set Two", "Norm Difference"], dataUtils.getRawAndNormalised(
                        f"../Data/output/fa/7mer_counts0.fa", f"../Data/output/fa/7mer_counts1.fa"),
                    outputFile="7mer_counts_raw_and_normalised")

                if parsed.s_l:
                    proportionsStatsTest("7mer_counts_raw_and_normalised", significanceLevel=parsed.s_l)

                else:
                    proportionsStatsTest("7mer_counts_raw_and_normalised")

        parsed.stats_for_kmers = False

    elif parsed.find_entropy:
        print("Finding entropy")
        findEntropy()
        parsed.find_entropy = False

    elif parsed.find_gc:
        print("Finding GC content")
        findGC()
        parsed.find_gc = False

    elif parsed.find_similar:
        print("Finding similar sequences")
        findSimilar()
        parsed.similarity_calculator = False


if __name__ == "__main__":
    arguments = argparse.ArgumentParser()
    dataUtils = DataUtils()

    arguments.add_argument("-command_file", help="Use this command if you have a pre-defined"
                                                 "command file. See the README for further information.")

    arguments.add_argument("-f", nargs="+", help="Enter filepath(s) of file(s) to analyse")

    arguments.add_argument("-get_unmapped", help="Set to True to find unmapped reads. To -f pass first the "
                                                 "file of original reads, and second the file of reads from assembly."
                                                 "Produces a fastq file as output, with name 'unmapped_reads'. Use "
                                                 "-fastq_name to change this. ")
    arguments.add_argument("-fastq_name", help="Name of the fastq file outputted by -get_unmapped.")

    arguments.add_argument("-kmers", help="Calls Jellyfish on files specified in -f, returns .fa files and csv files "
                                          "of them.", required=False)
    arguments.add_argument("-k", help="The value of k for jellyfish to search with: 7 as default.", required=False)

    arguments.add_argument("-compare_kmers", help="Compares sets of kmers in files provided.\n"
                                                  "Uses Jellyfish to get the kmers.\n"
                                                  "Returns a CSV file of the raw counts and difference,\n"
                                                  " as well as the normalised values.",
                           required=False)
    arguments.add_argument("-csv_filename", help="Name of the output csv file from call. Used by count_gc, "
                                                 "compare_kmers, and similarity_calculator.",
                           required=False)

    arguments.add_argument("-stats_for_kmers", help="Takes the CSV file from -compare_kmers. \n"
                                                    "Performs hypothesis test on the normalised values from random "
                                                    "subset.\n "
                                                    "Null Hypothesis: The proportions of the occurrences are not "
                                                    "significantly different. \n "
                                                    "Alt Hypothesis: The proportions are significantly different. \n"
                                                    "Produces CSV file in Data/output directory.",
                           required=False)
    arguments.add_argument("-s_l", help="Sets the significance level for statistical tests.\n"
                                        "Defualt: 0.05",
                           required=False)
    arguments.add_argument("-o", help="Specifes output filename for hypothesis tests csv. "
                                      "Default: hypothesis_test_results",
                           required=False)

    arguments.add_argument("-find_entropy", help="This finds the entropy for all of the files inputted. "
                                                 "Returns a line chart showing all of the entropies of the sequences "
                                                 "in the file, and strip plot showing the same. "
                                                 "Returns a CSV file of the sequences that lie outside of the IQR.")

    arguments.add_argument("-yLabel", help="Y label for graph produced by tool.")
    arguments.add_argument("-xLabel", help="X label for graph produced by tool.")
    arguments.add_argument("-title", help="Title for graph produced by tool.")
    arguments.add_argument("-line_chart", help="Set to True for a line chart to be produced.")
    arguments.add_argument("-output_line_plot_name", help="Sets the name for the plot outputted.")
    arguments.add_argument("-strip_plot", help="If set to false, the strip plot will not be made.")
    arguments.add_argument("-output_strip_plot_name", help="Name for the strip plot.")

    arguments.add_argument("-find_gc", help="Calls the main function of the GC counting class. "
                                            "Always produces a CSV file of the sequences and their overall GC content."
                                            "By default, produces a bar chart and histogram of the data.")
    arguments.add_argument("-bar_chart", help="Set to False to stop bar chart being created by GC content.")
    arguments.add_argument("-histogram", help="Set to False to stop histogram being created by GC content, "
                                              "or Hamming Distance of similarity.")
    arguments.add_argument("-bin_size", help="Set the bin size for the GC counts to be binned into.")
    arguments.add_argument("-hist_name", help="Sets the name for the histogram outputted by -find_gc and "
                                              "-similarity_calculator")
    arguments.add_argument("-bar_name", help="Sets the name for the bar chart outputted by -find_gc")

    arguments.add_argument("-find_similar", help="Compares sequences in first file to those in second. Returns"
                                                          " a CSV file of the sequence in the first file with"
                                                          " closest matches from the second.")
    arguments.add_argument("-number_of_closest", help="The number of closest matches to find. Default is 3.")
    arguments.add_argument("-sample_first", help="True to sample the first file inputted. Default True.")
    arguments.add_argument("-sample_second", help="True to sample the second file inputted. Default True.")
    arguments.add_argument("-sample_percentage", help="The percentage to sample from the files. Useful for large files."
                                                      " Default 1 percent.")
    arguments.add_argument("-sample_seed", help="Seed to control the samples. Default is completely random.")
    arguments.add_argument("-cutoff", help="The cutoff to stop searching for the closest match at. Default 0.6.")

    arguments.add_argument("-assemble_and_find_unmapped", help="This uses Megahit to assemble the file(s) specified in "
                                                               "-f. These files must be in .fasta format. Then uses "
                                                               "Bowtie to find the mapped reads, which are then "
                                                               "outputted. Finally, finds the unmapped reads using the "
                                                               "file(s) of input reads, and the mapped reads.")
    arguments.add_argument("-assembly_file_name", help="Name for Bowtie genome. Default bowtieGenome.")
    arguments.add_argument("-mapped_reads_file_name", help="Name for file outputted by Bowtie. Default bowtieoutput.")
    arguments.add_argument("-unmapped_reads_file_name", help="Name for file of unmapped reads. Default unmappedReads")

    arguments.add_argument("-blastn", help="Set to True to use the blastn program")
    arguments.add_argument("-blastx", help="Set to True to use the blastx program")
    arguments.add_argument("-blast_output_name",
                           help="Specifies the name of the blast output. Default blast_output.csv")
    arguments.add_argument("-outfmt", help="The output format specifier to pass to blast. Default 10")
    arguments.add_argument("-make_database",
                           help="Set to True to make custom database. Default, blast uses the swissprot database")
    arguments.add_argument("-file_for_custom_db", help="Pass filepath to this to use for custom blast database")
    arguments.add_argument("-db_name", help="Name for custom database")
    arguments.add_argument("-db_type", help="Type of custom database")
    arguments.add_argument("-threads", help="Number of threads for blast to search with. Default 6")

    parsed = arguments.parse_args()

    if parsed.command_file:
        loadCommandFile(parsed.command_file)

    else:
        run()
