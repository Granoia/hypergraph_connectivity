#This file finds the connectivity outliers from the survey and looks them up in the data mapping. 

import connectivity_survey

HGRAPH_NAME = connectivity_survey.INFILE
INFILE = "cdata_" + HGRAPH_NAME + ".csv"
INFILE_PATH = connectivity_survey.OUTFILE_PATH

REF_PATH = "/home/franzese/hypergraph_datafiles/NCIPID_data_mapping/"

pathify = lambda s: REF_PATH+s
depathify = lambda s: s[len(REF_PATH):]

node_file = "signaling-hypergraph-hypernodes.txt"
hedge_file = "signaling-hypergraph-hyperedges.txt"
elements_file = "ncipid-elements.txt"
reactions_file = "ncipid-reactions.txt"
subpathways_file = "ncipid-subpathways.txt"
controls_file = "ncipid-controls.txt"
complexes_file = "ncipid-complexes.txt"

ls = [elements_file, complexes_file, subpathways_file, controls_file, reactions_file]

ref_ls = [pathify(s) for s in ls]


def find_entry(id, reference_files=ref_ls):
    #given a pid_number, searches through ncipid-elements.txt to find the entry associated with the id
    #returns that entry
    for filename in reference_files:
        with open(filename,'r') as f:
            for line in f:
                if line.split('\t')[0] == id:
                    return line.strip() + '\t' + depathify(filename) + '\n'
    s = id + " not found.\n"
    return s

def find_colmax_node(col):
    with open(INFILE_PATH + INFILE, 'r') as f:
        f.readline()
        l = f.readline().strip().split(',')
        max_node = l[0]
        max_val = float(l[col])
        for line in f:
            l = line.strip().split(',')
            if float(l[col]) > max_val:
                max_node = l[0]
                max_val = float(l[col])

        return max_node, max_val


class Row:
    def __init__(self, name, bcon, fcon):
        self.name = name
        self.bcon = int(bcon)
        self.fcon = int(fcon)

    def __str__(self):
        return self.name

def catalog_rows(file = INFILE_PATH + INFILE):
    row_ls = []
    with open(file, 'r') as f:
        f.readline()
        for line in f:
            l = line.strip().split(',')
            row_ls.append(Row(l[0],l[1],l[2]))
    return row_ls
            

foo = lambda x: x.bcon


def main():
    row_ls = catalog_rows()
    row_ls.sort(key=lambda row: row.bcon, reverse=True)
    print("##########################")
    print("#TOP 5 B-CONNECTED NODES:#")
    print("##########################")
    for i in range(5):
        print('-----------------------------------------------')
        print(find_entry(str(row_ls[i])))
        print("B-connected nodes: ",row_ls[i].bcon)
        print('-----------------------------------------------')

    print
    print

    row_ls.sort(key=lambda row: row.fcon, reverse=True)
    print("##########################")
    print("#TOP 5 F-CONNECTED NODES:#")
    print("##########################")
    for i in range(5):
        print('-----------------------------------------------')
        print(find_entry(str(row_ls[i])))
        print("F-connected nodes: ",row_ls[i].fcon)
        print('-----------------------------------------------')
        
    #max_node, max_val = find_colmax_node(1)
    #entry = find_entry(max_node, ref_ls)
    #print(entry.split('\t')[1],max_val)

if __name__ == "__main__":
    main()
