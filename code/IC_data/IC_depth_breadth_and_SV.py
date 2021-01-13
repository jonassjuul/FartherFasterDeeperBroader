import networkx as nx 
import math
import numpy as np 
from collections import Counter

def max_depth(G,root_node) :
    p=nx.shortest_path_length(G,source=root_node)
    return max(p.values())

def max_breadth(G,root_node) :
    p=nx.shortest_path_length(G,source=root_node)
    count = Counter(list(p.values()))
    most_common = count.most_common(1)
    return most_common[0][1]

'''def Structural_virality(G) :
    WI = nx.wiener_index(G)
    number_of_nodes = len(G.nodes())
    if (number_of_nodes>1) :
        return WI/((number_of_nodes)*(number_of_nodes-1))
    else :
        return WI'''

def Structural_virality(G):
    size = len(G.nodes())
    if size==1:
        return 0 ##virality is not defined for cascades of size 1,
    sv=nx.average_shortest_path_length(G)  #Note: this is very time-consuming for larger cascades
    return sv


def depth_breadth_SV(G,root_node) :
    # Depth
    p = nx.shortest_path_length(G,source=root_node)
    depth = max(p.values())

    # Breadth
    count = Counter(list(p.values()))
    most_common = count.most_common(1)
    breadth = most_common[0][1]

    # Structural Virality
    SV = Structural_virality(G)

    return [depth, breadth, SV]


Rspread=0.800000
cascade_size = 1
# Make file for results
f = open('Measures_ICCascades_Rspread%s_MinimumCascadeSize%i.txt'%(Rspread,cascade_size),'w')
f.write('Depth\tbreadth\Structural Virality')
f.close()


# Open file and make analysis
f = open('ICCascades_Rspread%.6f_MinimumCascadeSize%i.txt'%(Rspread,cascade_size),'r')
f.readline()
line_num = -1
for line in f : 
    line_num +=1
    print("Doing line number",line_num)
    # Import only nodes that get neigbors
    line.strip()
    columns = line.split(' ')
    columns1 = [x for x in columns if not '\n' in x]
    columns2 = [x for x in columns1 if not 'k' in x]
    growth_array = [node.replace('s','') for node in columns2]


    # Create graph
    G = nx.Graph()
    G.add_node(0)
    node_num = 0
    for growth_node in growth_array :

        node_num+=1
        G.add_edge(int(growth_node),node_num)

    # Calculate results
    results = depth_breadth_SV(G,0)

    # Save results
    f = open('Measures_ICCascades_Rspread%s_MinimumCascadeSize%i.txt'%(Rspread,cascade_size),'a')
    f.write('\n%i\t%i\t%.4f'%(results[0],results[1],results[2]))
    f.close()
