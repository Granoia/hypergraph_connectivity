import survey_summary

NODEFILE = "signaling-hypergraph-hypernodes.txt"
HEDGEFILE = "signaling-hypergraph-hyperedges.txt"

INFILE_PATH = "/home/franzese/hypergraph_datafiles/"


def parse_hedges(file=INFILE_PATH + HEDGEFILE):
    #extracts data out of a text file and organizes it into a list of hedge objects which is returned
    hedges_ls = []
    with open(file, 'r') as f:
        f.readline()
        for line in f.readlines():
            h_ls = line.split("\t")
            tails = h_ls[0].split(";")
            heads = h_ls[1].split(";")
            posReg = h_ls[2].split(";")
            negReg = h_ls[3].split(";")
            ID = h_ls[4]
            #pathways = h_ls[5].split(";")      #pathways field is in Reactome but not NCIPID
            #hedges_ls.append(hedge(heads,tails,posReg,negReg,ID,pathways))

            hedges_ls.append(hedge(heads,tails,posReg,negReg,ID))
            
    return hedges_ls


class hedge:
    def __init__(self,head,tail,posReg,negReg,ID,pathways=None):
        self.head_ls = head
        self.tail_ls = tail
        self.posReg = posReg
        self.negReg = negReg
        self.ID = ID
        self.pathways = pathways


def parse_nodes(filename=INFILE_PATH + NODEFILE):
    #extracts information from the list of hypernodes and puts them into node objects.
    node_ls = []
    with open(filename, 'r') as file:
        file.readline()
        for line in file.readlines():
            curr_ls = line.split("\t")
            node_ls.append(node(curr_ls[0]))
    node_ls.sort(key = lambda n: n.name)
    return node_ls

class node:
    #node objects for the graph G. they know who they are adjacent to and whether they are in a fragment
    def __init__(self, name):
        self.name = name
        self.entrances = []
        self.exits = []
        self.pRegulates = []
        self.pRegulatedBy = []
        self.nRegulates = []
        self.nRegulatedBy = []
        self.adj_nodes = []
        self.hedges = []
        self.distance = -1
        self.parent = None
        self.in_frag = None
        self.index = None
        self.lowlink = None    #index and lowlink are for Tarjan's algorithm, implemented in tarjan.py

    def __str__(self):
        n = lambda node: node.name
        s = ""
        s += "Name: "+str(self.name)
        s += "\tEntrances: "+str([n(x) for x in self.entrances])
        s += "\tExits: "+str([n(x) for x in self.exits])
        s += "\tpRegulates: "+str([n(x) for x in self.pRegulates])+"\tpRegulatedBy: "+str([n(x) for x in self.pRegulatedBy])
        s += "\tnRegulates: "+str([n(x) for x in self.nRegulates])+'\tnRegulatedBy: '+str([n(x) for x in self.nRegulatedBy])
        s += '\n' +'Adjacent Nodes: ' + str([n(x) for x in self.adj_nodes])
        return s
    

def populate_help(node_ls, hedge_ls,curr_hedge, part):
    #helper function for populate_nodes()
    if part == "tail_ls":
        lookup = curr_hedge.tail_ls
    elif part == "posReg":
        lookup = curr_hedge.posReg
    elif part == "negReg":
        lookup = curr_hedge.negReg
    
    for item in lookup:
        if item != "None":
            tail_node = binary_search_names(node_ls, item)
            if tail_node == None:
                print("Error in first step of populate_help! Failed to look up " + item)
            else:    
                for h in curr_hedge.head_ls:
                    if h != "None":
                        head_node = binary_search_names(node_ls, h)
                        if head_node == None:
                            print("Error in second step of populate_help! Failed to look up " + h)
                        tail_node.hedges.append(curr_hedge)
                        head_node.hedges.append(curr_hedge)
                        
                        tail_node.adj_nodes.append(head_node)
                        head_node.adj_nodes.append(tail_node)
                        

                        if part == "tail_ls":
                            tail_node.exits.append(head_node)
                            head_node.entrances.append(tail_node)
                        elif part == "posReg":
                            tail_node.pRegulates.append(head_node)
                            head_node.pRegulatedBy.append(tail_node)
                        elif part == "negReg":
                            tail_node.nRegulates.append(head_node)
                            head_node.nRegulatedBy.append(tail_node)
        
def populate_nodes(node_ls, hedge_ls, regulators = False):
    #attaches the node objects together using a list of hedges (which can be obtained using a function in parseCount.py). In doing so converts the hypergraph back to a standard graph.
    for hedge in hedge_ls:
        populate_help(node_ls, hedge_ls, hedge, "tail_ls")
        if regulators == True:
            populate_help(node_ls, hedge_ls, hedge, "posReg")
            #populate_help(node_ls, hedge_ls, hedge, "negReg")

def binary_search_names(ls, target, start=0, end=None):
    #binary searches a list of nodes organized by name
    if end == None:
        end = len(ls)

    mid = (end+start)//2
    if ls[mid].name == target:
        return ls[mid]
    elif start == end:
        return None
    elif ls[mid].name < target:
        return binary_search_names(ls, target, mid+1, end)
    elif ls[mid].name > target:
        return binary_search_names(ls, target, start, mid)



###################################
#compound graph conversion section#
###################################

#the original code for parse_nodes() results in a compound graph
#here I will write code which when given the compound graph will return the broken down version of the graph
#since data mappings tend to vary between databases, this code will likely be specific to the NCIPID system


def scroll_through_nodes(nodes):
    #this is just a diagnostic function that determines which files all the lookups are coming from in the nodes
    #every id from the ncipid nodes comes from the mapping file for elements, subpathways, or complexes
    file_set = set()
    for n in nodes:
        entry = survey_summary.find_entry(n.name)
        file_set.add(entry.strip().split('\t')[-1])
    print(file_set)


def scroll_through_complexes():
    #another diagnostic function
    #found that recursive division of complexes is not necessary.
    file_set = set()
    suspicious = []
    recursive = []
    not_suspicious_but_recursive = []
    with open(survey_summary.ref_ls[1],'r') as f:
        f.readline()
        for line in f:
            compound_elements = line.strip().split('\t')[4].split(';')
            for e in compound_elements:
                entry = survey_summary.find_entry(e)
                entry_file = entry.strip().split('\t')[-1]
                file_set.add(entry_file)
                if entry_file == 'ncipid-complexes.txt':
                    recursive.append(entry)
                    if entry.strip().split('\t')[4] != 'None':
                        not_suspicious_but_recursive.append(entry)
    print(file_set)
    for r in recursive:
        print r
    print(len(recursive))

    print
    print
    print('##############################')
    print('#not suspicious but recursive#')
    print('##############################')

    for e in not_suspicious_but_recursive:
        print e
    print(len(not_suspicious_but_recursive))




def is_compound(n):
    #uses the data mapping to determine whether a given node is a compound or an atomic element
    #if atomic element, returns False
    #if compound, returns a list of individual atomic elements in the compound
    #some entries in the ncipid-complexes file list no component elements. I wasn't exactly sure what should be done with these but I decided to treat them as atomic elements
    entry = survey_summary.find_entry(n.name)
    entry_file = entry.strip().split('\t')[-1]
    if entry_file == 'ncipid-complexes.txt':
        component_str = entry.strip().split('\t')[4]
        if component_str == 'None':
            return False
        component_ls = component_str.split(';')
        return component_ls
    else:
        return False


def unify(a, b):
    #returns the set union of two lists as a list
    a_set = set(a)
    b_set = set(b)
    return list(a_set | b_set)


def find_by_name(name, nodes):
    for n in nodes:
        if n.name == name:
            return n
    return None



def unify_attributes(n, c):
    #updates missing connections for several node attributes
    n.entrances = unify(n.entrances, c.entrances)
    n.exits = unify(n.exits, c.exits)
    n.pRegulates = unify(n.pRegulates, c.pRegulates)
    n.pRegulatedBy = unify(n.pRegulatedBy, c.pRegulatedBy)
    n.nRegulates = unify(n.nRegulates, c.nRegulates)
    n.nRegulatedBy = unify(n.nRegulatedBy, c.nRegulatedBy)
    n.adj_nodes = unify(n.adj_nodes, c.adj_nodes)
    n.hedges = unify(n.hedges, c.hedges)

def apply_to_attributes(f, n):
    n.entrances = f(n.entrances)
    n.exits = f(n.exits)
    n.pRegulates = f(n.pRegulates)
    n.pRegulatedBy = f(n.pRegulatedBy)
    n.nRegulates = f(n.nRegulates)
    n.nRegulatedBy = f(n.nRegulatedBy)
    n.adj_nodes = f(n.adj_nodes)
    n.hedges = f(n.hedges)

def takeaway(ls, element):
    new = [x for x in ls if x != element]
    return new

def convert_compound(c, component_ls, nodes):
    obj_components = []
    for n_str in component_ls:
        n_obj = find_by_name(n_str, nodes)
        if n_obj: 
            obj_components.append(n_obj)
        else:
            #if one of the nodes in the complex doesn't yet exist on the graph, make a new node object and add it to the graph
            n_obj = node(n_str)
            obj_components.append(n_obj)
            nodes.append(n_obj)
            
            
    #this adds the complex's pointers to the component elements of the broken up complex
    for n in obj_components:
        unify_attributes(n,c)

    #these for loops alter the objects that the complex was pointing at so that they reflexively point at the component elements of the broken up complex
    for n in c.entrances:
        n.exits = unify(n.exits, obj_components)
    for n in c.exits:
        n.entrances = unify(n.entrances, obj_components)

    for n in c.pRegulates:
        n.pRegulatedBy = unify(n.pRegulatedBy, obj_components)
    for n in c.pRegulatedBy:
        n.pRegulates = unify(n.pRegulates, obj_components)

    for n in c.nRegulates:
        n.nRegulatedBy = unify(n.nRegulatedBy, obj_components)
    for n in c.nRegulatedBy:
        n.nRegulates = unify(n.nRegulates, obj_components)

    for n in c.adj_nodes:
        n.adj_nodes = unify(n.adj_nodes, obj_components)
        

    

    #removes all pointers that point at the broken up complex
    def remover(ls):
        return takeaway(ls, c)
    for n in nodes:
        apply_to_attributes(remover, n)

    return takeaway(nodes, c) #gives back the node list without the complex

    


def non_compound_graph(nodes):
    for n in nodes:
        compound_check = is_compound(n)
        print(compound_check)
        if compound_check:
            nodes = convert_compound(n, compound_check, nodes)
    return nodes
            


    



def find_frags(G):
    #uses BFS to find the weakly connected components in G.
    #returned as a list of fragment objects
    frags=[]
    for node in G:
        if node.in_frag == None:
            frags.append(frag_BFS(G,node))
    return frags

def BFS_help(n,curr,Q,frag_ls):
    if n.distance == -1:
        n.distance = curr.distance + 1
        n.parent = curr
        n.in_frag = True
        Q.enqueue(n)
        frag_ls.append(n)

def frag_BFS(G,root,weak_component=True,reg_connect=True):
    #runs BFS starting from root and ignoring edge directionality to find a weakly connected component
    frag_ls=[root]            #starting from root and compile each node reached this way into a returned fragment object.
    Q = queue()
    root.distance = 0
    Q.enqueue(root)
    while not Q.is_empty():
        curr = Q.dequeue()
        if not weak_component and reg_connect:
            for n in unify(curr.exits, curr.pRegulates):
                BFS_help(n,curr,Q,frag_ls)
        elif not weak_component and not reg_connect:
            for n in curr.exits:
                BFS_help(n,curr,Q,frag_ls)
        else:
            for n in curr.adj_nodes:
                BFS_help(n,curr,Q,frag_ls)
                
        """
        else:
            for n in curr.adj_nodes:
                if n.distance == -1:
                    n.distance = curr.distance + 1
                    n.parent = curr
                    n.in_frag = True
                    Q.enqueue(n)
                    frag_ls.append(n)
        """
    new_frag = fragment(frag_ls)
    return new_frag

def reset_nodes(nodes):
    for n in nodes:
        n.distance = -1
        n.parent = None
        n.in_frag = None

def BFS_each_node(nodes, outfile,reg_connect=True):
    with open(outfile, 'w') as out:
        out.write('#name, reachable nodes\n')
        for n in nodes:
            current_frag = frag_BFS(nodes,n, False, reg_connect)
            out.write(str(n.name)+','+str(current_frag.size)+'\n')
            reset_nodes(nodes)

class fragment:
    #a fragment object catalogs a single weakly connected component on G
    def __init__(self,nodes):
        self.node_ls = nodes
        self.size = len(nodes)

    def find_node_by_name(self, target_name):
        for node in self.node_ls:
            if node.name == target_name:
                return node
        print("No node with given name in this fragment")
    
    def find_hedges(self):         #gets a list of all hedges included in a fragment
        all_hedges = []      #has duplicates
        frag_hedges = []     #will have duplicates filtered out
        for node in self.node_ls:
            all_hedges += node.hedges
        all_hedges.sort(key = lambda h: h.ID)
        i = 0
        all_len = len(all_hedges)
        if len(all_hedges) != 0:
            singleton = all_hedges[0]
            frag_hedges.append(singleton)
            while i < all_len:
                if all_hedges[i].ID != singleton.ID:
                    singleton = all_hedges[i]
                    frag_hedges.append(singleton)
                i += 1
        return frag_hedges

    def hedge_size(self):         #finds the number of hedges in the fragment
        frag_hedges = self.find_hedges()
        return len(frag_hedges)


class queue:
    def __init__(self):
        self.ls = []

    def is_empty(self):
        if len(self.ls) == 0:
            return True
        else:
            return False
    
    def enqueue(self, item):
        self.ls.append(item)

    def dequeue(self):
        if self.is_empty():
            return None
        else:
            return self.ls.pop(0)


def get_frag_sizes(frag_ls,min=1):
    sizes = []
    for frag in frag_ls:
        if frag.size >= min:
            sizes.append(frag.size)
    return sizes





#aggregated executive functions

def construct_compound_graph(nodefile, hedgefile, preg=False):
    hedges = parse_hedges(hedgefile)
    nodes = parse_nodes(nodefile)
    populate_nodes(nodes, hedges, preg)
    return nodes

def construct_noncompound_graph(nodefile, hedgefile, preg=False):
    nodes = construct_compound_graph(nodefile, hedgefile, preg)
    return non_compound_graph(nodes)







def main():
    hedges = parse_hedges()
    nodes = parse_nodes()
    populate_nodes(nodes, hedges, True)
    #nodes = non_compound_graph(nodes)
    
    #frags = find_frags(nodes)
    #frags.sort(key = lambda n: n.size)
    #frag_sizes = get_frag_sizes(frags)
    #frag_sizes.sort()
    #print(frag_sizes)
    #print(len(nodes))
    #print
    #print

    BFS_each_node(nodes, 'test_cdata_should_be_noreg_cgraph.txt', reg_connect=False)
    #scroll_through_nodes(nodes)
    #scroll_through_complexes()
    

if __name__ == "__main__":
    main()
