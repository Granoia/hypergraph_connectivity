#plots a histogram according to the data obtained by connectivity_survey.py
#much of this code was repurposed from /home/bavent/cbb_tests/signaling-hypergraphs/plot_unregulated-signaling.py


import sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sys.path.append("/home/franzese/halp/halp/halp")
from directed_hypergraph import DirectedHypergraph
import utilities.directed_statistics as stats

import connectivity_survey

HGRAPH_FILENAME = connectivity_survey.INFILE
INFILE = "cdata_" + HGRAPH_FILENAME + ".csv"
INFILE_PATH = "/home/franzese/hypergraph_connectivity/cdata_outfiles/"
OUTFILE_PATH = "/home/franzese/hypergraph_connectivity/survey_plots/"



def main():
    with open(INFILE_PATH + INFILE, 'r') as f:
        v = np.loadtxt(f, delimiter=",", dtype='int', comments="#", skiprows=1, usecols=range(1,3))

        data_columns = [("B-connected nodes", 0), ("F-connected nodes", 1)]

        for data_column in data_columns:
            title = data_column[0]
            column_number = data_column[1]

            print(title + " distribution")
            print("-" * len(title + " distribution"))

            v_hist = np.ravel([row[column_number] for row in v])   # 'flatten' v
            max_val = max(v_hist)
            min_val = min(v_hist)
            fig = plt.figure()
            plt.title(HGRAPH_FILENAME+"\n" + title + " distribution")
            plt.xlabel("number of {}".format(title))
            plt.ylabel("frequency")
            ax1 = fig.add_subplot(111)
            hist, bins = np.histogram(v_hist, bins=(max_val-min_val+1), range=(min_val-.5, max_val+.5))
            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            plt.bar(center, hist, align='center', width=width)
            plt.xlim([min_val -.5, max_val+.5])

            print("Min: {}".format(min(v_hist)))
            print("Max: {}".format(max(v_hist)))
            print("Mean: {}".format(np.mean(v_hist)))
            print("Median: {}\n".format(np.median(v_hist)))

            plt.savefig(OUTFILE_PATH + HGRAPH_FILENAME +"_"+ title.replace(' ','_') + "_distribution.png")
            plt.close()


        
        #the code below is from bavent's original file. I think makes a log transformed version of the above plot, but I haven't repurposed it yet since I'm not planning on making a log transformed plot. But it may be useful in the near future so I'm leaving it here. 
        """
        # plt.show()
        fig = plt.figure()
        plt.title("Unregulated Signaling Pathway\n" + title + " distribution")
        plt.xlabel("number of {}".format(title))
        plt.ylabel("frequency")
        ax1 = fig.add_subplot(111)
        hist, bins = np.histogram(v_hist, bins=(max_val-min_val+1), range=(min_val-.5, max_val+.5))
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width, log=True)
        plt.xlim([min_val -.5, max_val+.5])
        plt.yscale("symlog")
        plt.savefig("log-png_plots/log-" + title + " distribution.png")
        plt.close()
        """


if __name__ == "__main__":
    main()

