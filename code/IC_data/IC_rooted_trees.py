import networkx as nx 
import math
import numpy as np 
from collections import Counter

def make_digraph_from_edges(G_edges) :
    G = nx.DiGraph()
    for edge in G_edges :
        G.add_edge(edge[0],edge[1])
        
    return G




# Rooted trees...
cascade_dictionary = {}

# size 1
cascade_dictionary[1] = {}

G1 = nx.DiGraph()
G1.add_node(0)
cascade_dictionary[1][1] = G1


# size 2
cascade_dictionary[2] = {}

G2edges = [[0,1]]
G2 = make_digraph_from_edges(G2edges)
cascade_dictionary[2][2] = G2

# size 3
cascade_dictionary[3] = {}

G3edges = [[0,1],[0,2]]
G3 = make_digraph_from_edges(G3edges)
cascade_dictionary[3][3] = G3

G4edges = [[0,1],[1,2]]
G4 = make_digraph_from_edges(G4edges)
cascade_dictionary[3][4] = G4

# size 4
cascade_dictionary[4] = {}

G5edges = [[0,1],[0,2],[0,3]]
G5 = make_digraph_from_edges(G5edges)
cascade_dictionary[4][5] = G5

G6edges = [[0,1],[0,2],[1,3]]
G6 = make_digraph_from_edges(G6edges)
cascade_dictionary[4][6] = G6

G7edges = [[0,1],[1,2],[2,3]]
G7 = make_digraph_from_edges(G7edges)
cascade_dictionary[4][7] = G7

G8edges = [[0,1],[1,2],[1,3]]
G8 = make_digraph_from_edges(G8edges)
cascade_dictionary[4][8] = G8

# size 5
cascade_dictionary[5] = {}

G9edges = [[0,1],[0,2],[0,3],[0,4]]
G9 = make_digraph_from_edges(G9edges)
cascade_dictionary[5][9] = G9

G10edges = [[0,1],[0,2],[0,3],[1,4]]
G10 = make_digraph_from_edges(G10edges)
cascade_dictionary[5][10] = G10

G11edges = [[0,1],[0,2],[2,3],[1,4]]
G11 = make_digraph_from_edges(G11edges)
cascade_dictionary[5][11] = G11

G12edges = [[0,1],[0,2],[1,3],[3,4]]
G12 = make_digraph_from_edges(G12edges)
cascade_dictionary[5][12] = G12

G13edges = [[0,1],[0,2],[1,3],[1,4]]
G13 = make_digraph_from_edges(G13edges)
cascade_dictionary[5][13] = G13

G14edges = [[0,1],[1,2],[1,3],[2,4]]
G14 = make_digraph_from_edges(G14edges)
cascade_dictionary[5][14] = G14

G15edges = [[0,1],[1,2],[2,3],[2,4]]
G15 = make_digraph_from_edges(G15edges)
cascade_dictionary[5][15] = G15

G16edges = [[0,1],[1,2],[2,3],[3,4]]
G16 = make_digraph_from_edges(G16edges)
cascade_dictionary[5][16] = G16

G17edges = [[0,1],[1,2],[1,3],[1,4]]
G17 = make_digraph_from_edges(G17edges)
cascade_dictionary[5][17] = G17


# size 6
cascade_dictionary[6] = {}

G18edges = [[0,1],[0,2],[0,3],[0,4],[0,5]]
G18 = make_digraph_from_edges(G18edges)
cascade_dictionary[6][18] = G18

G19edges = [[0,1],[0,2],[0,3],[0,4],[1,5]]
G19 = make_digraph_from_edges(G19edges)
cascade_dictionary[6][19] = G19

G20edges = [[0,1],[0,2],[0,3],[1,4],[2,5]]
G20 = make_digraph_from_edges(G20edges)
cascade_dictionary[6][20] = G20

G21edges = [[0,1],[0,2],[0,3],[1,4],[1,5]]
G21 = make_digraph_from_edges(G21edges)
cascade_dictionary[6][21] = G21

G22edges = [[0,1],[0,2],[0,3],[1,4],[4,5]]
G22 = make_digraph_from_edges(G22edges)
cascade_dictionary[6][22] = G22

G23edges = [[0,1],[0,2],[1,3],[2,4],[2,5]]
G23 = make_digraph_from_edges(G23edges)
cascade_dictionary[6][23] = G23

G24edges = [[0,1],[0,2],[1,3],[3,4],[2,5]]
G24 = make_digraph_from_edges(G24edges)
cascade_dictionary[6][24] = G24

G25edges = [[0,1],[0,2],[1,3],[1,4],[3,5]]
G25 = make_digraph_from_edges(G25edges)
cascade_dictionary[6][25] = G25

G26edges = [[0,1],[0,2],[1,3],[3,4],[3,5]]
G26 = make_digraph_from_edges(G26edges)
cascade_dictionary[6][26] = G26

G27edges = [[0,1],[0,2],[1,3],[3,4],[4,5]]
G27 = make_digraph_from_edges(G27edges)
cascade_dictionary[6][27] = G27

G28edges = [[0,1],[0,2],[1,3],[1,4],[1,5]]
G28 = make_digraph_from_edges(G28edges)
cascade_dictionary[6][28] = G28

G29edges = [[0,1],[1,2],[1,3],[1,4],[2,5]]
G29 = make_digraph_from_edges(G29edges)
cascade_dictionary[6][29] = G29

G30edges = [[0,1],[1,2],[1,3],[2,4],[2,5]]
G30 = make_digraph_from_edges(G30edges)
cascade_dictionary[6][30] = G30

G31edges = [[0,1],[1,2],[1,3],[2,4],[3,5]]
G31 = make_digraph_from_edges(G31edges)
cascade_dictionary[6][31] = G31

G32edges = [[0,1],[1,2],[1,3],[2,4],[4,5]]
G32 = make_digraph_from_edges(G32edges)
cascade_dictionary[6][32] = G32

G33edges = [[0,1],[1,2],[2,3],[2,4],[2,5]]
G33 = make_digraph_from_edges(G33edges)
cascade_dictionary[6][33] = G33

G34edges = [[0,1],[1,2],[2,3],[2,4],[3,5]]
G34 = make_digraph_from_edges(G34edges)
cascade_dictionary[6][34] = G34

G35edges = [[0,1],[1,2],[2,3],[3,4],[3,5]]
G35 = make_digraph_from_edges(G35edges)
cascade_dictionary[6][35] = G35

G36edges = [[0,1],[1,2],[2,3],[3,4],[4,5]]
G36 = make_digraph_from_edges(G36edges)
cascade_dictionary[6][36] = G36

G37edges = [[0,1],[1,2],[1,3],[1,4],[1,5]]
G37 = make_digraph_from_edges(G37edges)
cascade_dictionary[6][37] = G37

# Make dictionary to keep results...
result_dic = {}

for size in [1,2,3,4,5,6] :
    result_dic[size] = {}
    for graph_num in cascade_dictionary[size].keys() :
        result_dic[size][graph_num] = 0    



Rspread=0.800000
cascade_size = 1
# Make file for results


# Open file and make analysis
f = open('ICCascades_Rspread%.6f_MinimumCascadeSize%i.txt'%(Rspread,cascade_size),'r')
f.readline()
line_num = -1
for line in f : 
    #if (line_num == 30) :
    #    break
    line_num +=1
    #if (line_num/100 == line_num//100) :
    print("Doing line number",line_num)
    # Import only nodes that get neigbors
    line.strip()
    columns = line.split(' ')
    columns1 = [x for x in columns if not '\n' in x]
    columns2 = [x for x in columns1 if not 'k' in x]
    growth_array = [node.replace('s','') for node in columns2]

    if (len(growth_array)<=5) :
        # Create graph
        G = nx.DiGraph()
        G.add_node(0)
        node_num = 0
        for growth_node in growth_array :

            node_num+=1
            G.add_edge(int(growth_node),node_num)


        # Now check which rooted tree....
        size_graph = len(G.nodes())

        if (size_graph <= 6) :

            if (len(cascade_dictionary[size_graph])>1) :
                for graph_num in cascade_dictionary[size_graph].keys() :
                    if ( nx.is_isomorphic(cascade_dictionary[size_graph][graph_num],G) == True ) :
                        #print(size_graph)
                        result_dic[size_graph][graph_num] +=1
                        break
            else : 
                only_key = int(list(cascade_dictionary[size_graph].keys())[0])
                result_dic[size_graph][only_key] += 1            



    # Calculate results
    #results = depth_breadth_SV(G,0)

# Save results
destination='../Outputs/'
np.save(destination+'ICRootedTrees_Rspread%s_MinimumCascadeSize%i.npy'%(Rspread,cascade_size),result_dic)
