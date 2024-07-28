import matplotlib.pyplot as plt
import matplotlib
import seaborn

matplotlib.use('TkAgg')


class PlotsUtils:
    def __init__(self):
        self.convertedData = []
        self.header = ""

    def makeHistogram(self, data, yLabel, xLabel, title, ax, barColour="blue"):
        """
        Plots a histogram on the axes object passed to it, using the data passed.
        :param data: The data to plot.
        :param yLabel: Y label for the plot
        :param xLabel: X label for the plot
        :param title: Title for the plot
        :param ax: The axes to plot the histogram on
        :param barColour: The colour of the bars
        """
        ax.hist(data.dropna(), bins="auto", log=False)
        ax.set_title(f"{title}")
        ax.set_xlabel(f"{xLabel}")
        ax.set_ylabel(f"{yLabel}")
        ax.plot()

    def makeBarChart(self, xData, yData, ax, yLabel="Default Y Label", xLabel="Default X Label", title="Default title"):
        """
        Plots a histogram on the axes object passed to it, using the data passed.
        :param xData: The data to plot on the x-axis.
        :param yData: The data to plot on the y-axis. Passed to the "height" parameter
        :param yLabel: Y label for the plot
        :param xLabel: X label for the plot
        :param title: Title for the plot
        :param ax: The axes to plot the bar chart on
        """
        plt.setp(ax.get_xticklabels(), fontsize=8)
        ax.bar(x=xData, height=yData)
        ax.set_ylabel(f"{yLabel}")
        ax.set_xlabel(f"{xLabel}")
        ax.set_title(f"{title}")
        ax.plot()

    def makeLineChart(self, data, yLabel="Default", xLabel="Default", title="Default", filename="line_chart"):
        """
        A function that creates a line chart with the data passed to it.
        :param data: The data for use on the y axis
        :param yLabel: The Y label for the plot
        :param xLabel: The X label for the plot
        :param title: The title for the plot
        :param filename: The name to save the plot
        """
        fig, ax = plt.subplots(figsize=(30, 12))
        yData = data.values()
        xData = [i for i in range(0, len(data))]

        ax.plot(xData, yData, linewidth=0.5)
        ax.set_title(f"{title}")
        ax.set_ylabel(f"{yLabel}")
        ax.set_xlabel(f"{xLabel}")
        ax.plot()
        plt.savefig(f"../Data/output/plots/{filename}.png")

    def makeStripPlot(self, dataframe, ax):
        """
        Makes a strip plot using the dataframe passed to it, upon the axis provided.
        :param dataframe: The dataframe containing the data
        :param ax: The axes object to plot the strip plot to.
        :return: A strip plot object.
        """
        return seaborn.stripplot(dataframe, ax=ax, orient="h", size=0.5, native_scale=False)

