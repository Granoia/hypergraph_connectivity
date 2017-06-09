#Given the name of a hypergraph infile located in /home/franzese/hypergraph_datafiles, this code creates a survey of the connectivity of that hypergraph
#For each node in the hypergraph, take a B-visit and F-visit (note: F-visit is currently an incorrect algorithm, I have a proof of this in my notebook), and record how many nodes are found by each.

#Saves the summary as a .csv file in the cdata_outfiles subdirectory
#each line of the csv file looks like this:
#node, number of b-visited nodes, number of f-visited nodes

#Much of this code was repurposed from /home/bavent/cbb_tests/signaling-hypergraphs/test_unregulated_signaling.py

INFILE = "NCIPID_p_regulators"
INFILE_PATH = "/home/franzese/hypergraph_datafiles/"
OUTFILE_PATH = "/home/franzese/hypergraph_connectivity/cdata_outfiles/"


import sys

sys.path.append("/home/franzese/halp/halp/halp")
from directed_hypergraph import DirectedHypergraph

import algorithms.directed_paths as directed_paths


def main():
        H = DirectedHypergraph()
        H.read(INFILE_PATH + INFILE + ".txt", ";", "\t")

        out_file = open(OUTFILE_PATH + "cdata_" + INFILE + ".csv", 'w')
        out_file.write("source node,b-connected nodes,f-connected nodes\n")

        for source_node in H.get_node_set():
	        line = str(source_node) + ","
        
                # B-conn execution
                b_visited_nodes, Pv, Pe, v = directed_paths.b_visit(H, source_node)
                line += str(len(b_visited_nodes)) + ","

                # F-conn execution
                f_visited_nodes, Pv, Pe, v = directed_paths.f_visit(H, source_node)
                line += str(len(f_visited_nodes)) + "\n"

                out_file.write(line)

if __name__ == "__main__":
        main()
