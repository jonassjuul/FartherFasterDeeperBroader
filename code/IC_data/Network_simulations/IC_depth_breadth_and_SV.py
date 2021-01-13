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

#mechanism = ['uniform','preferential attachment','copying'][1]

Rspread=0.550000
cascade_size = 1
# Make file for results

filenames = [
'Cascades_on_network_ICMODELathletes_edges_Rspread0.350000_MinimumCascadeSize1.txt',
'Cascades_on_network_ICMODELathletes_edges_Rspread0.400000_MinimumCascadeSize1.txt',
'Cascades_on_network_ICMODELCornell5_Rspread0.350000_MinimumCascadeSize1.txt',
'Cascades_on_network_ICMODELCornell5_Rspread0.400000_MinimumCascadeSize1.txt'
#'Cascades_on_network_ICMODELMSU24_Rspread0.350000_MinimumCascadeSize1.txt',
#'Cascades_on_network_ICMODELMSU24_Rspread0.400000_MinimumCascadeSize1.txt',
#'Cascades_on_network_ICMODELTexas84_Rspread0.350000_MinimumCascadeSize1.txt',
#'Cascades_on_network_ICMODELTexas84_Rspread0.400000_MinimumCascadeSize1.txt'
]

for filename in filenames :
    f = open('Measures_'+filename,'w')
    #f.write('\n%i\t%i\t%.4f'%(results[0],results[1],results[2]))
    f.write('Depth\tbreadth\Structural Virality')
    f.close()


    # Open file and make analysis
    f = open(filename,'r')
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
        growth_array = [[],[]]
        entry_num = -1


        for entry in columns2 :
            entry_num +=1
            new_columns = entry.split('s')

            if (entry_num == 0) :
                seed = int(new_columns[1])
                G = nx.Graph()

                G.add_node(seed)
            else : 
                growth_array[0].append(int(new_columns[0]))
                growth_array[1].append(int(new_columns[1]))

        #growth_array = [node.replace('s','') for node in columns2]



        # Create graph
        node_num = 0
        for turn in range (len(growth_array[0])) :

            #node_num+=1
            G.add_edge(growth_array[0][turn],growth_array[1][turn])

        # Calculate results
        results = depth_breadth_SV(G,seed)

        # Save results
        f = open('Measures_'+filename,'a')
        f.write('\n%i\t%i\t%.4f'%(results[0],results[1],results[2]))
        f.close()

