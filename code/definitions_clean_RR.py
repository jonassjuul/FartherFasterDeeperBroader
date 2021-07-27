import numpy as np 
import matplotlib.pyplot as plt 
import math 
import random
from scipy import stats
import datetime
import scipy
import copy
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import pandas as pd
import networkx as nx 
import matplotlib.patches as mpatches
from matplotlib.offsetbox import (DrawingArea, OffsetImage,AnnotationBbox)
# ------------------
# Functions
# ------------------

def match_to(treatment,control,match_to='treatment') :
    
    cascades_treatment = {}
    cascades_control = {}
    
    if (match_to == 'treatment') :
        for size in treatment.keys() :
            if (size in control.keys()) :
                cascades_treatment[size] = list(treatment[size].keys())
                cascades_control[size] = random.choices(list(control[size].keys()),k=len(cascades_treatment[size]))
                
    return cascades_treatment,cascades_control
    
def get_plot_arrays(res_arr) :
    values = np.arange(0,max(res_arr)+1,1)
    frequencies = np.zeros(len(values))
    
    for value in res_arr : 
        frequencies[value] += 1
    return values,frequencies

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

# Match to treatment (True)
def sample_topology(N_experiments,cascade_data,treatment_name = 'TRUE',control_name = 'FALSE',variables=['size','depth','max_breadth','virality'],also_below_3 = True) :
    go_ahead = False        
    if (also_below_3 == True) :
        go_ahead = True


    pvalues = {}#{'size':[],'depth':[],'max_breadth':[],'virality':[]}
    for variable in variables :
        pvalues[variable] = []
    
    for experiment in range (N_experiments) :
        if (experiment/100 == experiment//100) :
            print("Doing experiment",experiment,"out of",N_experiments)
        
        # Draw random subset of cascades
        random_treatment, random_control = match_to(treatment=cascade_data[treatment_name],control=cascade_data[control_name])


        # Make result dic
        result_dic = {treatment_name:{},control_name:{}}
        for key in result_dic.keys() :
            result_dic[key] = {}
            for variable in variables :
                result_dic[key][variable] = []
            #result_dic[key] = {'size':[],'depth':[],'max_breadth':[],'virality':[]}
        
        # Now extract interesting features from randomly selected cascades..
        for size in random_treatment.keys():
            if (go_ahead == True or size > 2) :
                for cascade_name in random_treatment[size] :
                    for topological_feature in result_dic[treatment_name].keys():
                        result_dic[treatment_name][topological_feature].append(cascade_data[treatment_name][size][cascade_name][topological_feature])


                if (treatment_name != control_name) :
                    for cascade_name in random_control[size] :
                        for topological_feature in result_dic[control_name].keys():
                            result_dic[control_name][topological_feature].append(cascade_data[control_name][size][cascade_name][topological_feature])
                    
        for feature in result_dic[treatment_name].keys() :
            if (feature != 'prufer_sequence') :
                # Calculate KS statistics
                KSstat,KSpvalue = scipy.stats.ks_2samp(result_dic[treatment_name][feature],result_dic[control_name][feature])

                # Save KS pvalue
                pvalues[feature].append(KSpvalue)

    return pvalues,result_dic

# SORT FROM SMALLEST TO LARGEST
def sort_cascades(result_dic) :
    keys = list(result_dic.keys())
    for topological_feature in result_dic[keys[0]].keys() :
        print(topological_feature)
        for key in keys :
            result_dic[key][topological_feature] = sorted(result_dic[key][topological_feature])
        
    return result_dic


# Plot the CCDF of arr
def plot_ccdf(arr,color) :
    plt.plot(arr,1-np.arange(0,len(arr),1)/len(arr),color)

# Time for spread.
# Match to treatment (True)

def calc_gmean(arr) :
    return scipy.stats.gmean(np.array(arr)/60)

def sample_time(N_experiments,cascade_data,treatment_name = 'TRUE',control_name = 'FALSE') :
        
    var_name = 'uu2time'
    
    for experiment in range (N_experiments) :
        if (experiment/100 == experiment//100) :
            print("Doing experiment",experiment,"out of",N_experiments)
        
        random_treatment, random_control = match_to(treatment=cascade_data[treatment_name],control=cascade_data[control_name])

        if (experiment == 0) :
            final_result_dic = {treatment_name:{},control_name:{}}
            for key in final_result_dic.keys() :
                for size in range (max(random_treatment.keys())):
                    final_result_dic[key][size] = []          

        # Make result dic
        result_dic = {treatment_name:{},control_name:{}}
        for key in result_dic.keys() :
            for size in range (max(random_treatment.keys())+1):
                result_dic[key][size] = []
        
        # Now extract interesting features from randomly selected cascades..
        for size in random_treatment.keys():
            for cascade_name in random_treatment[size] :
                for subsize in (cascade_data[treatment_name][size][cascade_name][var_name].keys()) :
                    #print(size,subsize)
                    result_dic[treatment_name][subsize].append(cascade_data[treatment_name][size][cascade_name][var_name][subsize])

            if (treatment_name != control_name) :

                for cascade_name in random_control[size] :
                    for subsize in (cascade_data[control_name][size][cascade_name][var_name].keys()) :

                        result_dic[control_name][subsize].append(cascade_data[control_name][size][cascade_name][var_name][subsize])
                    
        # Calculate geometric mean
        for size in final_result_dic[treatment_name].keys() :
            final_result_dic[treatment_name][size].append(calc_gmean(result_dic[treatment_name][size])+0)
            if (treatment_name != control_name):
                final_result_dic[control_name][size].append(calc_gmean(result_dic[control_name][size])+0)
        
    
    return final_result_dic,result_dic



# Make temporal analysis
def temporal_analysis(gmean_dic) :
    median_time = {'TRUE':[],'FALSE':[],'sizes':[]}
    time95 = {'TRUE':[],'FALSE':[],'sizes':[]}
    time05 = {'TRUE':[],'FALSE':[],'sizes':[]}


    for size in range (1,max(gmean_dic['TRUE'].keys())+1) :
        if (size in gmean_dic['TRUE'].keys()) :
            median_time['TRUE'].append(np.median(gmean_dic['TRUE'][size]))

            
            
            
            median_time['FALSE'].append(np.median(gmean_dic['FALSE'][size]))
            time95['FALSE'].append(np.percentile(gmean_dic['FALSE'][size],95))
            time05['FALSE'].append(np.percentile(gmean_dic['FALSE'][size],5))        
            median_time['sizes'].append(size)

    for key in median_time : 
        median_time[key] = np.array(median_time[key])
        time95[key] = np.array(time95[key])
        time05[key] = np.array(time05[key])

    return median_time,time05,time95

def temporal_analysis_single (cascade_data) :

    gmean_diconlyfalse,example_time=sample_time(1,cascade_data,treatment_name='FALSE',control_name='FALSE')
    gmean_diconlytrue,example_time=sample_time(1,cascade_data,treatment_name='TRUE',control_name='TRUE')


    median_timesingle = {'TRUE':[],'FALSE':[],'TRUEsizes':[],'FALSEsizes':[]}
    for size in range (1,max(gmean_diconlytrue['TRUE'].keys())+1) :
        if (size in gmean_diconlytrue['TRUE'].keys()) :
            median_timesingle['TRUE'].append(np.median(gmean_diconlytrue['TRUE'][size]))
            median_timesingle['TRUEsizes'].append(size)

    print(max(gmean_diconlyfalse['FALSE'].keys()))        
    for size in range (1,max(gmean_diconlyfalse['FALSE'].keys())+1) :

        if (size in gmean_diconlyfalse['FALSE'].keys()) :        
            
            median_timesingle['FALSE'].append(np.median(gmean_diconlyfalse['FALSE'][size]))
            median_timesingle['FALSEsizes'].append(size)

    return median_timesingle

def temporal_analysis_categories(gmean_dic_categories) :

    median_time_categories = {'Other':[],'Politics':[],'sizes':[]}
    time95_categories = {'Other':[],'Politics':[],'sizes':[]}
    time05_categories = {'Other':[],'Politics':[],'sizes':[]}


    for size in range (1,max(gmean_dic_categories['Other'].keys())+1) :
        if (size in gmean_dic_categories['Other'].keys()) :
            median_time_categories['Other'].append(np.median(gmean_dic_categories['Other'][size]))

            
            
            
            median_time_categories['Politics'].append(np.median(gmean_dic_categories['Politics'][size]))
            time95_categories['Politics'].append(np.percentile(gmean_dic_categories['Politics'][size],95))
            time05_categories['Politics'].append(np.percentile(gmean_dic_categories['Politics'][size],5))        
            median_time_categories['sizes'].append(size)

    for key in median_time_categories : 
        median_time_categories[key] = np.array(median_time_categories[key])
        time95_categories[key] = np.array(time95_categories[key])
        time05_categories[key] = np.array(time05_categories[key])

    return median_time_categories,time05_categories,time95_categories


def temporal_analysis_single_categories (cascade_data_categories) :

    gmean_diconlypolitics,example_time=sample_time(1,cascade_data_categories,treatment_name='Politics',control_name='Politics')
    gmean_diconlyother,example_time=sample_time(1,cascade_data_categories,treatment_name='Other',control_name='Other')



    median_timesingle_categories = {'Other':[],'Politics':[],'Othersizes':[],'Politicssizes':[]}
    for size in range (1,max(gmean_diconlyother['Other'].keys())+1) :
        if (size in gmean_diconlyother['Other'].keys()) :
            median_timesingle_categories['Other'].append(np.median(gmean_diconlyother['Other'][size]))
            median_timesingle_categories['Othersizes'].append(size)

    print(max(gmean_diconlypolitics['Politics'].keys()))        
    for size in range (1,max(gmean_diconlypolitics['Politics'].keys())+1) :

        if (size in gmean_diconlypolitics['Politics'].keys()) :        
            
            median_timesingle_categories['Politics'].append(np.median(gmean_diconlypolitics['Politics'][size]))
            median_timesingle_categories['Politicssizes'].append(size)

    return median_timesingle_categories

def temporal_analysis_FirstFinal(gmean_dic) :
    median_time = {'First':[],'Final':[],'sizes':[]}
    time95 = {'First':[],'Final':[],'sizes':[]}
    time05 = {'First':[],'Final':[],'sizes':[]}


    for size in range (1,max(gmean_dic['First'].keys())+1) :
        if (size in gmean_dic['First'].keys()) :
            median_time['First'].append(np.median(gmean_dic['First'][size]))

            
            
            
            median_time['Final'].append(np.median(gmean_dic['Final'][size]))
            time95['Final'].append(np.percentile(gmean_dic['Final'][size],95))
            time05['Final'].append(np.percentile(gmean_dic['Final'][size],5))        
            median_time['sizes'].append(size)

    for key in median_time : 
        median_time[key] = np.array(median_time[key])
        time95[key] = np.array(time95[key])
        time05[key] = np.array(time05[key])

    return median_time,time05,time95


def temporal_analysis_single_FirstFinal (cascade_data) :

    gmean_diconlyfalse,example_time=sample_time(1,cascade_data,treatment_name='Final',control_name='Final')
    gmean_diconlytrue,example_time=sample_time(1,cascade_data,treatment_name='First',control_name='First')


    median_timesingle = {'First':[],'Final':[],'Firstsizes':[],'Finalsizes':[]}
    for size in range (1,max(gmean_diconlytrue['First'].keys())+1) :
        if (size in gmean_diconlytrue['First'].keys()) :
            median_timesingle['First'].append(np.median(gmean_diconlytrue['First'][size]))
            median_timesingle['Firstsizes'].append(size)

    print(max(gmean_diconlyfalse['Final'].keys()))        
    for size in range (1,max(gmean_diconlyfalse['Final'].keys())+1) :

        if (size in gmean_diconlyfalse['Final'].keys()) :        
            
            median_timesingle['Final'].append(np.median(gmean_diconlyfalse['Final'][size]))
            median_timesingle['Finalsizes'].append(size)

    return median_timesingle


'''
# ------------------------
#  ROOTED TREES
# ------------------------
'''
def make_digraph_from_edges(G_edges) :
    G = nx.DiGraph()
    for edge in G_edges :
        G.add_edge(edge[0],edge[1])
        
    return G

def import_rootedtrees():
    # Create dictionary for cascades. 
    #     - First key: size
    #     - Second key: name
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

    return cascade_dictionary

# statistical tests
# MLE COMPARISON OF MULTINOMIALS 


# Calculate Pis.....
def MLE_multi_prob(A,B) :    
    # A is array with count of event i on entry i.
    # B likewise but for data set B
    
    # Example from test [Link: http://ugrad.stat.ubc.ca/~stat305/Mar26.pdf ]
    # See also separate notebook    
    
    nA = sum(A)
    nB = sum(B)
    Probs = []
    for i in range (len(A)) :
        Probs.append((A[i]+B[i])/(nA+nB))
    return Probs

def Chi2_multi(A,B) :
    # Example from test [Link: http://ugrad.stat.ubc.ca/~stat305/Mar26.pdf ]
    # See also separate notebook
    Pi = MLE_multi_prob(A,B)
    nA = sum(A)
    nB = sum(B)
    Chi2 = 0
    for i in range (len(A)) :
        
        if (Pi[i] >0) :
            Chi2 += (A[i]-Pi[i]*nA)*(A[i]-Pi[i]*nA)/(Pi[i]*nA)
            Chi2 += (B[i]-Pi[i]*nB)*(B[i]-Pi[i]*nB)/(Pi[i]*nB)
    return Chi2

def get_Chi2_pvalue(A,B,dof) :
    Chi2 = Chi2_multi(A,B)
    p = 1-stats.chi2.cdf(Chi2,dof)
    return p
    
def combine_p_values(p_vec) :
    print(p_vec)
    Y = -2*np.log(np.array(p_vec))
    print(Y)
    
    T = sum(Y)
    dof = 2*len(p_vec)
    
    p_combined = 1-stats.chi2.cdf(T,dof)
    return p_combined
    
def count_rooted_trees() :
    
    # import cascade data
    
    raw_data_anon_path = '' # Insert path to Vosoughi et al.'s file 'raw_data_anon.csv'
    if (raw_data_anon_path == '') :
        input('Path to raw_anon_data.csv must be specified. The data can be obtained using the link specified in the Acknowledgement section of Vosoughi et al. (DOI: 10.1126/science.aap9559)')
    
    metadata_file = raw_data_anon_path
    df= pd.read_csv(metadata_file)


    cascade_dictionary = import_rootedtrees()

    
    
    
    # STARS 
    # -------
    # Make dic for results...
    # ---------------------------

    stars_result_dic = {}
    result_dic = {}

    # Now go through all cascades.
    # ---------------------------

    # First is just a single node...
    current_cascade_id = df['cascade_id'][0]
    current_graph = nx.DiGraph()
    current_graph.add_node(df['tid'][0])
    veracity = 'FALSE'
    size_graph = 1
    current_root_node = 0
    is_star = True




    for line in range (1,len(df['cascade_id'])) :
        if (veracity not in stars_result_dic.keys()) :
            stars_result_dic[veracity] = {}   

        if (veracity not in result_dic.keys()) :
            result_dic[veracity] = {}
            
            for size in [1,2,3,4,5,6] :
                result_dic[veracity][size] = {}
                for graph_num in cascade_dictionary[size].keys() :
                    result_dic[veracity][size][graph_num] = 0           
            
    
            
        if (df['cascade_id'][line] == current_cascade_id) :
    
            current_graph.add_edge(df['parent_tid'][line],df['tid'][line])
            
            if (df['parent_tid'][line]!=current_root_node and df['parent_tid'][line] != -1) :
                is_star = False           

        else : 
            size_graph = len(current_graph.nodes())
            #if (size_graph <= 6) :
            if (size_graph not in stars_result_dic[veracity].keys()) :
                stars_result_dic[veracity][size_graph] = {'star':0,'total':0}

            stars_result_dic[veracity][size_graph]['total'] +=1
            if (is_star == True) :
                stars_result_dic[veracity][size_graph]['star'] +=1


            # Try finding motifs w size <= 6 here...
            # -----
            if (size_graph <= 6) :

                if (len(cascade_dictionary[size_graph])>1) :
                    for graph_num in cascade_dictionary[size_graph].keys() :
                        if ( nx.is_isomorphic(cascade_dictionary[size_graph][graph_num],current_graph) == True ) :
                            #print(size_graph)
                            result_dic[veracity][size_graph][graph_num] +=1
                            break
                else : 
                    only_key = int(list(cascade_dictionary[size_graph].keys())[0])
                    result_dic[veracity][size_graph][only_key] += 1        
            # -----
                
                
                
            current_graph = nx.DiGraph()
            current_graph.add_node(df['tid'][line])
            
            current_root_node = df['tid'][line] +0
            is_star = True

            current_cascade_id = df['cascade_id'][line]

            veracity = df['veracity'][line] 
        

    return result_dic,stars_result_dic

def rooted_trees_statistics(result_dic,keys) : 
    plot_arr = {keys[0]:{},keys[1]:{}}
    for veracity in [keys[0],keys[1]] :
        for size in result_dic[veracity].keys() :
            if(size not in plot_arr[veracity].keys()) :
                plot_arr[veracity][size] = []
            
            for key in result_dic[veracity][size].keys() :
                plot_arr[veracity][size] += [key]*result_dic[veracity][size][key]#/sum(result_dic[veracity][3].values())
                
    plot_vecs = {keys[0]:{},keys[1]:{}}
    for veracity in [keys[0],keys[1]] :
        for size in result_dic[veracity].keys() :
            if(size not in plot_vecs[veracity].keys()) :
                plot_vecs[veracity][size] = [[],[],[],[],[]]# name, val, std
            
            for key in result_dic[veracity][size].keys() :
                plot_vecs[veracity][size][0].append(key)# += [key]*result_dic[veracity][size][key]#/sum(result_dic[veracity][3].values())
                plot_vecs[veracity][size][1].append(result_dic[veracity][size][key]/sum(list(result_dic[veracity][size].values())))#/sum(list(result_dic[veracity][size].values())))
                plot_vecs[veracity][size][2].append(result_dic[veracity][size][key])#/sum(list(result_dic[veracity][size].values())))            
                plot_vecs[veracity][size][3].append(sum(list(result_dic[veracity][size].values()))*plot_vecs[veracity][size][1][-1]*(1-plot_vecs[veracity][size][1][-1]))

    # Get right errors...
    for veracity in [keys[0],keys[1]] :
        for size in result_dic[veracity].keys() :
            if( size not in plot_vecs[veracity].keys()) :
                plot_vecs[veracity][size] = [[],[],[],[]]# name, val, std
            
            for cascade_num in range (len(result_dic[veracity][size])) :
                var_val = 0
                for cascade_sum in range (len(result_dic[veracity][size])) :
                    if (cascade_num != cascade_sum ) :
                        var_val += abs(-plot_vecs[veracity][size][2][cascade_sum]/(sum(plot_vecs[veracity][size][2][:]))**2)**2*plot_vecs[veracity][size][3][cascade_sum]
                    
                    else : 
                        var_val += abs((sum(plot_vecs[veracity][size][2][:])-plot_vecs[veracity][size][2][cascade_sum])/(sum(plot_vecs[veracity][size][2][:]))**2)**2*plot_vecs[veracity][size][3][cascade_sum]
                        
                    
                
                plot_vecs[veracity][size][4].append(2*math.sqrt(var_val))


    p_TrueFalse_arr = []

    for size in [3,4,5,6] :
        TRUE_data = plot_vecs[keys[0]][size][2]
        FALSE_data = plot_vecs[keys[1]][size][2]    
        dof_data = len(plot_vecs[keys[0]][size][2]) -1
        p_TrueFalse = get_Chi2_pvalue(TRUE_data,FALSE_data,dof=dof_data)
        p_TrueFalse_arr.append(p_TrueFalse)

    p_TrueFalse_combined = combine_p_values(p_TrueFalse_arr)

    return plot_vecs,p_TrueFalse_combined,p_TrueFalse_arr

   
    
    
    
def count_rooted_trees_general(result_dic_SIR,keys) :
    
    # import cascade data
    raw_data_anon_path = '' # Insert path to Vosoughi et al.'s file 'raw_data_anon.csv'
    if (raw_data_anon_path == '') :
        input('Path to raw_anon_data.csv must be specified. The data can be obtained using the link specified in the Acknowledgement section of Vosoughi et al. (DOI: 10.1126/science.aap9559)')
    
    metadata_file = raw_data_anon_path
    df= pd.read_csv(metadata_file)


    cascade_dictionary = import_rootedtrees()

    
    
    
    # STARS 
    # -------
    # Make dic for results...
    # ---------------------------

    stars_result_dic = {}
    result_dic = {}

    # Now go through all cascades.
    # ---------------------------

    # First is just a single node...
    current_cascade_id = df['cascade_id'][0]
    current_graph = nx.DiGraph()
    current_graph.add_node(df['tid'][0])
    veracity = 'FALSE'
    size_graph = 1
    current_root_node = 0
    is_star = True




    for line in range (1,len(df['cascade_id'])) :
        if (veracity not in stars_result_dic.keys()) :
            stars_result_dic[veracity] = {}   

        if (veracity not in result_dic.keys()) :
            result_dic[veracity] = {}
            
            for size in [1,2,3,4,5,6] :
                result_dic[veracity][size] = {}
                for graph_num in cascade_dictionary[size].keys() :
                    result_dic[veracity][size][graph_num] = 0           
            
    
            
        if (df['cascade_id'][line] == current_cascade_id) :
    
            current_graph.add_edge(df['parent_tid'][line],df['tid'][line])
            
            if (df['parent_tid'][line]!=current_root_node and df['parent_tid'][line] != -1) :
                is_star = False           

        else : 
            size_graph = len(current_graph.nodes())
            #if (size_graph <= 6) :
            if (size_graph not in stars_result_dic[veracity].keys()) :
                stars_result_dic[veracity][size_graph] = {'star':0,'total':0}

            stars_result_dic[veracity][size_graph]['total'] +=1
            if (is_star == True) :
                stars_result_dic[veracity][size_graph]['star'] +=1


            # Try finding motifs w size <= 6 here...
            # -----
            if (size_graph <= 6) :

                if (len(cascade_dictionary[size_graph])>1) :
                    for graph_num in cascade_dictionary[size_graph].keys() :
                        if ( nx.is_isomorphic(cascade_dictionary[size_graph][graph_num],current_graph) == True ) :
                            #print(size_graph)
                            result_dic[veracity][size_graph][graph_num] +=1
                            break
                else : 
                    only_key = int(list(cascade_dictionary[size_graph].keys())[0])
                    result_dic[veracity][size_graph][only_key] += 1        
            # -----
                
                
                
            current_graph = nx.DiGraph()
            current_graph.add_node(df['tid'][line])
            
            current_root_node = df['tid'][line] +0
            is_star = True

            current_cascade_id = df['cascade_id'][line]

            veracity = df['veracity'][line] 
        

    return result_dic,stars_result_dic



        
'''
# ------------------------
#  TRUE / FALSE DATA
# ------------------------
'''

# IMPORTING DATA ..

#def import_truefalsedata_new() :
def import_truefalsedata() :
    raw_data_anon_path = '' # Insert path to Vosoughi et al.'s file 'raw_data_anon.csv'
    if (raw_data_anon_path == '') :
        input('Path to raw_anon_data.csv must be specified. The data can be obtained using the link specified in the Acknowledgement section of Vosoughi et al. (DOI: 10.1126/science.aap9559)')

    metadata_file = raw_data_anon_path
    print("Reading metadata...")
    print("Make sure that the Prufer Sequence has been saved for each cascade. Key should be 'prufer_sequence'. Networkx can help get the prufer sequence for each cascade when build with Vosoughi et al.'s code.")

    fin=open(metadata_file,'r')
    lines=fin.readlines()
    fin.close()
    cascade_id2metadata={}
    line_num = -1
    for line in lines:
        line_num += 1
        line=line.replace('\n','')
        item=eval(line)
        cascade_id2metadata[item[0]]=item[1]

    # Change format / Sort for size   
    cascade_data = {'FALSE':{},'TRUE': {},'MIXED':{}}
    cascade_data_categories = {'Politics':{},'Other': {},'Urban':{}}

    for cascade_name in cascade_id2metadata.keys() :
        
        # convenietly name cascade data
        this_cascade = cascade_id2metadata[cascade_name]
        
        
        
        # Define size and veracity variables
        veracity = this_cascade['veracity']
        size = this_cascade['size'] +0
        


        # If no cascades with this size have been saved so far, create key
        if (size not in cascade_data[veracity].keys()) :
            cascade_data[veracity][size] = {}
        
        # Create key for this particular cascade
        cascade_data[veracity][size][cascade_name] = {}
        for key in cascade_id2metadata[cascade_name].keys() :
            
            cascade_data[veracity][size][cascade_name][key] = cascade_id2metadata[cascade_name][key]
            if (key == 'virality' and cascade_id2metadata[cascade_name][key] == None) :
                cascade_data[veracity][size][cascade_name][key] = 0

        # Rumor topic
        rumor_category = this_cascade['rumor_category']
        if (rumor_category != 'Politics') :
            rumor_category = 'Other'        
        
        
        # If no cascades with this size have been saved so far, create key
        if (size not in cascade_data_categories[rumor_category].keys()) :
            cascade_data_categories[rumor_category][size] = {}

        # Create key for this particular cascade
        cascade_data_categories[rumor_category][size][cascade_name] = {}
        for key in cascade_id2metadata[cascade_name].keys() :
            
            cascade_data_categories[rumor_category][size][cascade_name][key] = cascade_id2metadata[cascade_name][key]
            if (key == 'virality' and cascade_id2metadata[cascade_name][key] == None) :
                cascade_data_categories[rumor_category][size][cascade_name][key] = 0

        # Rumor topic
        rumor_category = this_cascade['rumor_category']
        if (rumor_category == 'Viral Photos/Stories/Urban Legends') :
            rumor_category = 'Urban'        
        
        
            # If no cascades with this size have been saved so far, create key
            if (size not in cascade_data_categories[rumor_category].keys()) :
                cascade_data_categories[rumor_category][size] = {}

            # Create key for this particular cascade
            cascade_data_categories[rumor_category][size][cascade_name] = {}
            for key in cascade_id2metadata[cascade_name].keys() :

                cascade_data_categories[rumor_category][size][cascade_name][key] = cascade_id2metadata[cascade_name][key]
                if (key == 'virality' and cascade_id2metadata[cascade_name][key] == None) :
                    cascade_data_categories[rumor_category][size][cascade_name][key] = 0        


    return cascade_data, cascade_data_categories




# print KS statistics.
def KS_minimum_and_quantile(p_vec,features,quantile) :
    
    print("\n KS statistics :")
    print("----------------")
    for feature in features :
        print(feature,"\n\tminimum KS value:",min(p_vec[feature]),"\n\t%s percent of values above:"%((1-quantile)*100),np.quantile(p_vec[feature],quantile),"\n\n\tmaximum KS value:",max(p_vec[feature]),"\n\t%s percent of values below:"%((1-quantile)*100),np.quantile(p_vec[feature],1-quantile))
        
    print("\n----------------")        


def KS_make_table(cascades_unsampled,p_vec_sampled,features,quantile, Goel=False) :

    keys = list(cascades_unsampled.keys())

    # Make result dic
    result_dic = {keys[0]:{},keys[1]:{}}
    for key in result_dic.keys() :
        #result_dic[key] = {'size':[],'depth':[],'max_breadth':[],'virality':[]}
        if (Goel == True) :
            result_dic[key] = ['size','max.depth','virality']
        else : 
            result_dic[key] = {'size':[],'depth':[],'max_breadth':[],'virality':[]}        
    
    # Now extract interesting features from randomly selected cascades..
    for size in cascades_unsampled[keys[0]].keys():
        for cascade_name in cascades_unsampled[keys[0]][size] :
            for topological_feature in result_dic[keys[0]].keys():
                result_dic[keys[0]][topological_feature].append(cascades_unsampled[keys[0]][size][cascade_name][topological_feature])

    for size in cascades_unsampled[keys[1]].keys():
            for cascade_name in cascades_unsampled[keys[1]][size] :
                for topological_feature in result_dic[keys[1]].keys():
                    result_dic[keys[1]][topological_feature].append(cascades_unsampled[keys[1]][size][cascade_name][topological_feature])
                    

    print("\\begin{table}\n")
    print("\t\\begin{tabular}{l c c c c c}\\\\ \n")
    print("\t\t \\hline")
    print("\t\t Quantity & $p$ (US) & $\min{p}$ (S)  & $\max{p}$ (S) & 95pct of $p$ above (S) & 95pct of $p$ below (S)\\\\ ")
    print("\t\t \\hline")
    for feature in features :
        KSstat,KS_unsampled_pvalue = scipy.stats.ks_2samp(result_dic[keys[0]][feature],result_dic[keys[1]][feature])
        print(("\t\t %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\"%((feature.replace('_','-')).title(),KS_unsampled_pvalue,min(p_vec_sampled[feature]),max(p_vec_sampled[feature]),np.quantile(p_vec_sampled[feature],quantile),np.quantile(p_vec_sampled[feature],1-quantile))))
        #print("\t %s & %.2f & %.2f & %.2f & %.2f & %.2f")%(feature,KS_unsampled_pvalue,min(p_vec_sampled[feature]),np.quantile(p_vec_sampled[feature],quantile),max(p_vec_sampled[feature]),np.quantile(p_vec_sampled[feature],1-quantile))
    print("\t\t\\hline")
    print("\t\\end{tabular}")
    print("\t\\caption{INSERT CAPTION \\label{sup:tab:INSERT FIGNAME}}")
    print("\\end{table}")


def KS_make_table_degrees(KS_unsampled_pvalue,p_vec_sampled,quantile) :

    print("\\begin{table}\n")
    print("\t\\begin{tabular}{l c c c c c}\\\\ \n")
    print("\t\t \\hline")
    print("\t\t Quantity & $p$ (US) & $\min{p}$ (S)  & $\max{p}$ (S) & 95pct of $p$ above (S) & 95pct of $p$ below (S)\\\\ ")
    print("\t\t \\hline")
    #for feature in features :
    process_dict = {'TrueFalse':'True/False','SIR':'SIR','IC':'IC'}
    for process in ['TrueFalse','SIR','IC']:
    #KSstat,KS_unsampled_pvalue = scipy.stats.ks_2samp(result_dic[keys[0]][feature],result_dic[keys[1]][feature])
        print(("\t\t %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\"%("Out-degrees %s"%process_dict[process],KS_unsampled_pvalue[process],min(p_vec_sampled[process]),max(p_vec_sampled[process]),np.quantile(p_vec_sampled[process],quantile),np.quantile(p_vec_sampled[process],1-quantile))))
        #print("\t %s & %.2f & %.2f & %.2f & %.2f & %.2f")%(feature,KS_unsampled_pvalue,min(p_vec_sampled[feature]),np.quantile(p_vec_sampled[feature],quantile),max(p_vec_sampled[feature]),np.quantile(p_vec_sampled[feature],1-quantile))
    print("\t\t\\hline")
    print("\t\\end{tabular}")
    print("\t\\caption{INSERT CAPTION \\label{sup:tab:INSERT FIGNAME}}")
    print("\\end{table}")    

def get_cascade_scatter (cascade_data) :
    cascade_scatter = {'TRUE':{'size':[],'depth':[],'max_breadth':[],'virality':[]},'FALSE':{'size':[],'depth':[],'max_breadth':[],'virality':[]}}


    veracity_arr = ['TRUE','FALSE']


    for veracity in  veracity_arr :

        for size in cascade_data[veracity].keys() :
            size_arrays = {'depth':[],'max_breadth':[],'virality':[]}
            
            for cascade in cascade_data[veracity][size].keys() :
                for feature in cascade_scatter[veracity].keys() :
                    cascade_scatter[veracity][feature].append(cascade_data[veracity][size][cascade][feature])

                    
        for feature in cascade_scatter[veracity].keys() :
            cascade_scatter[veracity][feature] = np.array(cascade_scatter[veracity][feature])
                                                    
    return cascade_scatter

def get_cascade_lines(cascade_scatter,cascade_data) :
    cascade_lines = {'TRUE':{'depth':{},'max_breadth':{},'virality':{}},'FALSE':{'depth':{},'max_breadth':{},'virality':{}}}

    veracity_arr = ['TRUE','FALSE']

    plot_lines = ['mean',5,50,95]
    for veracity in cascade_lines.keys() :
        for feature in cascade_lines[veracity].keys() :
            for linetype in plot_lines :
                cascade_lines[veracity][feature][linetype]= []


    for veracity in veracity_arr : 
        cascade_lines[veracity]['size'] = []
        
        for size in sorted(cascade_data[veracity].keys()) :
            cascade_lines[veracity]['size'].append(size)
            size_arr = (cascade_scatter[veracity]['size'] == size)            
            for feature in cascade_scatter[veracity].keys() :         
                if (feature != 'size') :
                    for linetype in plot_lines :
                        if (linetype == 'mean') :
                            cascade_lines[veracity][feature][linetype].append(np.mean(cascade_scatter[veracity][feature][size_arr]))    
                        else : 
                            cascade_lines[veracity][feature][linetype].append(sorted(cascade_scatter[veracity][feature][size_arr])[int(len(cascade_scatter[veracity][feature][size_arr])*linetype/100)])
    return cascade_lines


'''
# ------------------------
#  SIR model
# ------------------------
'''


def model_topology_and_size(cascade_data_SIR,result_dic_key,filename_cascades,filename_measures,N_cascades):
    # Import SIR DATA
    cascade_data_SIR[result_dic_key] = {}

    SIR_First_dataset = [[],[],[]] #Depth, breadth, SV
    SIR_Second_dataset = [[],[],[]]

    SIR_First_cascade_sizes = []
    SIR_Second_cascade_sizes = []


    # GET SIZES OF CASCADES
    f = open(filename_cascades,'r')

    f.readline()
    line_num = -1

    sizes80 = {}

    for line in f : 

        line_num +=1
        if (line_num>=N_cascades) :
            break    

        line.strip()
        columns = line.split(' ')
        columns1 = [x for x in columns if not '\n' in x]
        columns2 = [x for x in columns1 if not 'k' in x]    
        size = len(columns)-1
        if (size not in cascade_data_SIR[result_dic_key].keys()) :
            cascade_data_SIR[result_dic_key][size] = {}
            
        cascade_data_SIR[result_dic_key][size][line_num] = {}
        
        sizes80[line_num] = size
        
    f.close()

    # GET CASCADE TOPOLOGIES

    f = open(filename_measures,'r')
    f.readline()
    line_num = -1
    for line in f :
        line_num += 1
        if (line_num>=N_cascades) :
            break
        line.strip()
        columns = line.split()
        
        cascade_data_SIR[result_dic_key][sizes80[line_num]][line_num] = {'size':sizes80[line_num],'depth':int(columns[0]),'max_breadth':int(columns[1]),'virality':float(columns[2])}
    f.close()

    return cascade_data_SIR

'''
# ------------------------
#  Network data
# ------------------------
'''

def get_network_data(cascades_IC_networks,cascades_IC_networks_sizes,IC_networks_files,IC_networks_files_measure,keys,school) :
    # ON NETWORKS


    # IC MODEL

    R0s = IC_networks_files[school].keys()

    cascades_IC_networks[school] = {keys[0]:{},keys[1]:{}} #= {'Cornell':{35:{},40:{}},'MSU':{35:{},40:{}},'Texas':{35:{},40:{}}}

    cascades_IC_networks_sizes[school] = {keys[0]:{},keys[1]:{}}#{'Cornell':{35:{},40:{}},'MSU':{35:{},40:{}},'Texas':{35:{},40:{}}}


    # GET SIZES OF CASCADES WITH SMALL R0
    for school in cascades_IC_networks.keys() :
        for R0 in cascades_IC_networks[school].keys() :
            filename = IC_networks_files[school][R0]
            f = open(filename,'r')

            f.readline()
            line_num = -1

            #sizes80 = {}

            for line in f : 

                line_num +=1
                if (line_num>=29999) :
                    break    

                line.strip()
                columns = line.split(' ')
                columns1 = [x for x in columns if not '\n' in x]
                columns2 = [x for x in columns1 if not 'k' in x]    
                size = len(columns)-1

                if (size not in cascades_IC_networks[school][R0].keys()) :
                    cascades_IC_networks[school][R0][size] = {}

                cascades_IC_networks[school][R0][size][line_num] = {}

                #sizes80[line_num] = size
                cascades_IC_networks_sizes[school][R0][line_num] = size

            f.close()



    # GET CASCADE TOPOLOGIES OF CASCADES WITH SMALL R0
    for school in cascades_IC_networks.keys() :
        for R0 in cascades_IC_networks[school].keys() :
            filename = IC_networks_files_measure[school][R0]
            f = open(filename,'r')

            f.readline()
            line_num = -1
            
            for line in f : 

                line_num +=1
                if (line_num>=29999) :
                    break    
                line.strip()
                columns = line.split()                                                            
                this_size = cascades_IC_networks_sizes[school][R0][line_num]+0
                cascades_IC_networks[school][R0][this_size][line_num] = {'size':this_size+0,'depth':int(columns[0]),'max_breadth':int(columns[1]),'virality':float(columns[2])}
            f.close()    

    return cascades_IC_networks,cascades_IC_networks_sizes

'''
# --------------------
# PLOTS
# --------------------
'''

def import_pictures() :
    im_size = plt.imread('drawings/size.png')
    im_depth = plt.imread('drawings/depth.png')
    im_max_breadth = plt.imread('drawings/max_breadth.png')
    im_virality = plt.imread('drawings/virality2.png')
    im_temporal = plt.imread('drawings/temporal.png')
    return im_size,im_depth,im_max_breadth,im_virality,im_temporal

# convert cm to inches
def cm_to_inch(cm) :
    return 0.3937007874*cm

divide_h = 40/3

def h_text(xlim,width='narrow') :

    if ( width == 'wide') :
        return (xlim[0] + (xlim[1]-xlim[0])/divide_h)
    
    elif ( width == 'medium') :
        return (xlim[0] + (xlim[1]-xlim[0])/divide_h*1.5)    
        
    else :     
    
        return (xlim[0] + np.log((xlim[1]-xlim[0]))/divide_h)#*3)
    
    
    
def v_text(ylim,width='narrow',logy = False) :
    if (logy == False) :
        if ( width == 'wide') :
            return(ylim[1] - (ylim[1]-ylim[0])/7)#10)

        else : 
            return(ylim[1] - (ylim[1]-ylim[0])/5)#7)
        
    else : 
        if ( width == 'wide') :
            return(ylim[1] - (ylim[1]-ylim[0])/7)#10)

        else : 
            return  (ylim[1] - np.log((ylim[1]-ylim[0]))*4050)#7)        

# TRUE / FALSE NEWS

def FalseTrue_ccdfs(topology_TRUE,topology_FALSE,example_cascades,median_timesingle,median_time,time05,time95,destinations) :



    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(17.8/6*5.6),cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,5)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    #CCDF_yticks = [0.00001,0.0001,0.001,0.01,0.1,1]

    CCDF_yticks = [0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = [1,10,100,1000,10000]
    col1_xticks = col0_xticks
    col2_xticks = [1,10]#[1,10,100]
    coldepth_xticks = [1,10]#[1,10,100]


    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxx_col0 = 100000
    maxx_col1 = 100000
    maxx_col2 = 14#100
    maxx_depth = 30#100

    #maxx_col3 = 70000

    maxy_CCDF = 1.3
    miny_CCDF = 0.000007
    #maxy_col1 = 1
    #maxy_col2 = 1
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    coldepth_xlabel = 'Cascade Depth'

    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Virality'
    col3_xlabel = 'Unique Users'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()


    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(topology_TRUE['%s'%'TRUE']['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth,label='True')


    plotthis = sorted(topology_FALSE['%s'%'FALSE']['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth,label='False')


    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([1,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    ax00.set_ylabel(CCDF_ylabel)
    #ax00.set_xlabel(col0_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        
    plt.text(0.04, 0.77, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)

    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False)

    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    # (0,1) Depth
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(topology_TRUE['%s'%'TRUE']['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(topology_FALSE['%s'%'FALSE']['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_depth])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col2_xticks)

    ax01.set_ylabel(CCDF_ylabel)
    #ax01.set_xlabel(coldepth_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
        

    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    



    # (0,2) BREADTH
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(topology_TRUE['%s'%'TRUE']['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(topology_FALSE['%s'%'FALSE']['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax02.set_xscale('log')
    ax02.set_yscale('log')
    ax02.set_xlim([1,maxx_col1])
    ax02.set_ylim([miny_CCDF,maxy_CCDF])

    ax02.set_yticks(CCDF_yticks)
    ax02.set_xticks(col1_xticks)

    ax02.set_ylabel(CCDF_ylabel)
    #ax02.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax02.set_yticklabels([])
        ax02.set_ylabel('')
        ax02.set_xlabel('')
        
    
    ax02.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    



    # (0,3) VIRALITY
    ax03 = plt.subplot2grid(plot_dimension,(0,3),rowspan=1,colspan=1)
    feature = 'virality'

    plotthis = sorted(topology_TRUE['%s'%'TRUE']['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(topology_FALSE['%s'%'FALSE']['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax03.set_xscale('log')
    ax03.set_yscale('log')
    ax03.set_xlim([1,maxx_col2])
    ax03.set_ylim([miny_CCDF,maxy_CCDF])
    ax03.set_yticks(CCDF_yticks)
    ax03.set_xticks(col2_xticks)

    ax03.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax03.set_yticklabels([])
        ax03.set_ylabel('')
        ax03.set_xlabel('')
        
        
 
    ax03.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # (0,4) TIME
    ax04 = plt.subplot2grid(plot_dimension,(0,4),rowspan=1,colspan=1)

    ax04.plot(median_timesingle['TRUEsizes'],median_timesingle['TRUE'],'g',linewidth=FT_linewidth)
    ax04.plot(median_timesingle['FALSEsizes'],median_timesingle['FALSE'],'r',linewidth=FT_linewidth)

    ax04.set_yscale('log')
    ax04.set_ylim([miny_time,maxy_time])

    ax04.set_xticks(ax03_xticks)
    ax04.set_yticks(time_yticks)

    #ax04.set_xlabel(col3_xlabel)
    ax04.set_ylabel(col3_ylabel)
 
    ax04.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  
    plt.text(0.807, 0.74, 'Minutes',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)


    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(example_cascades['%s'%'TRUE']['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades['%s'%'FALSE']['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([1,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)
    ax10.set_ylabel(CCDF_ylabel)

    ax10.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.04, 0.29+0.03, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)



    # INSERT SMALL IMAGE
    im_ax10 = inset_axes(ax10,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax10.imshow(im_size)
    im_ax10.axis('off')

    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(example_cascades['%s'%'TRUE']['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades['%s'%'FALSE']['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_depth])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col2_xticks)

    ax11.set_xlabel(coldepth_xlabel)

    ax11.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')

    ax11.set_ylabel('G',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax11 = inset_axes(ax11,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax11.imshow(im_depth)
    im_ax11.axis('off')

    # (1,2)
    ax12 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades['%s'%'TRUE']['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades['%s'%'FALSE']['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_xlim([1,maxx_col1])
    ax12.set_ylim([miny_CCDF,maxy_CCDF])
    ax12.set_yticks(CCDF_yticks)
    ax12.set_xticks(col1_xticks)

    ax12.set_xlabel(col1_xlabel)

    ax12.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax12.set_yticklabels([])
        ax12.set_ylabel('')
 
    ax12.set_ylabel('H',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax12 = inset_axes(ax12,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax12.imshow(im_max_breadth)
    im_ax12.axis('off')

    # (1,3)
    ax13 = plt.subplot2grid(plot_dimension,(1,3),rowspan=1,colspan=1)
    feature = 'virality'

    plotthis = sorted(example_cascades['%s'%'TRUE']['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades['%s'%'FALSE']['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax13.set_xscale('log')
    ax13.set_yscale('log')
    ax13.set_xlim([1,maxx_col2])
    ax13.set_ylim([miny_CCDF,maxy_CCDF])
    ax13.set_yticks(CCDF_yticks)
    ax13.set_xticks(col2_xticks)


    ax13.set_xlabel(col2_xlabel)

    ax13.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax13.set_yticklabels([])
        ax13.set_ylabel('')

            
    ax13.set_ylabel('I',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight) 

    # INSERT SMALL IMAGE
    im_ax13 = inset_axes(ax13,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax13.imshow(im_virality)
    im_ax13.axis('off')

    # (1,4)
    ax14 = plt.subplot2grid(plot_dimension,(1,4),rowspan=1,colspan=1)
    feature = 'virality'

    plotthis = sorted(example_cascades['%s'%'TRUE']['%s'%feature])
    ax14.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g')


    plotthis = sorted(example_cascades['%s'%'FALSE']['%s'%feature])
    ax14.plot(median_time['sizes'],median_time['TRUE'],color='g',linewidth=FT_linewidth)
    ax14.plot(median_time['sizes'],median_time['FALSE'],color='r',linewidth=FT_linewidth)
    ax14.fill_between(median_time['sizes'],time05['FALSE'],time95['FALSE'],color='r',alpha=alpha90)

    ax14.set_yscale('log')
    ax14.set_ylim([miny_time,maxy_time])
    ax14.set_xticks(ax13_xticks)
    ax14.set_yticks(time_yticks)

    ax14.set_xlabel(col3_xlabel)
    ax14.set_ylabel(col3_ylabel)
 
    ax14.set_ylabel('J',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax14 = inset_axes(ax14,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax14.imshow(im_temporal)
    im_ax14.axis('off')

    plt.text(0.807, 0.26+0.02, 'Minutes',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)




    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()

    for destination in destinations :
        plt.savefig(destination+'FigFalseTrue.png',dpi=400)

# Loop through all cascades..
def offspring_dist(dataset) :
    degrees = []
    degrees_nonroot = []
    degrees_root = []
    for size in dataset.keys() :
        for cascade in dataset[size].keys() :
            PS = dataset[size][cascade]['prufer_sequence']
            if (PS != None) :
                degree_sequence = nx.from_prufer_sequence(PS).degree()
                nodenum = 0
                for node in degree_sequence :
                    if (nodenum == 0) :
                        degrees.append(node[1]+1)
                        degrees_root.append(node[1]+1)
                    else : 
                        degrees.append(node[1])
                        degrees_nonroot.append(node[1])
                    nodenum += 1
            else : 
                degrees.append(0+1)
                degrees_root.append(0+1)
                
    return degrees, degrees_root, degrees_nonroot


# Loop through all cascades..
def offspring_dist_subsampled(dataset) :
    degrees = []
    degrees_nonroot = []
    degrees_root = []    
    
    sequences = dataset['prufer_sequence']
    for PS in sequences :
        if (PS != None) :
            degree_sequence = nx.from_prufer_sequence(PS).degree()
            nodenum = 0
            for node in degree_sequence :
        
                if (nodenum == 0) :
                    degrees.append(node[1]+1)
                    degrees_root.append(node[1]+1)
                    
                else : 
                    degrees.append(node[1])
                    degrees_nonroot.append(node[1])
                    
                nodenum += 1        
        else : 
            degrees.append(0+1)
            degrees_root.append(0+1)
        
                
    return degrees,  degrees_root, degrees_nonroot


def get_dataset_and_pruferset_models(filepath) :
    

    dataset = {}
    file = open(filepath,'r')
    cascade_num = 0
    for line in file :
        cascade_num +=1
        line.strip()
        column = line.split(' ')
        
        G = nx.Graph()
        G.add_node(0)
        
        next_node = 1
        for event in column :
            if (event[:-1]!='') :
                node = int(event[:-1])
                suffix = event[-1]


                if (suffix == 's') :
                    G.add_edge(node,next_node)
                    next_node +=1

        size = G.size()+1
        if (size>=2) :
            prufer_seq = nx.to_prufer_sequence(G)
        else : 
            prufer_seq = None
        
        if (size not in dataset.keys()) :
            dataset[size] = {}
        dataset[size][cascade_num] = {'size':size,'prufer_sequence':prufer_seq}
    
    
    return dataset

# TRUE / FALSE JOINT DISTRIBUTIONS
def FalseTrue_joint(cascade_scatter,cascade_lines,destinations) :
    figsize2 = (cm_to_inch(14/4*5),cm_to_inch(6))
    figsize_joint = (figsize2[1]/2*2.3,figsize2[0]/5*2.5)

    markersize = .5
    marker_color = {'TRUE':'g','FALSE':'r'}
    marker_alpha = 0.3


    line_color = {'TRUE':'g','FALSE':'r'}
    line_alpha=0.5

    plt.rcParams["legend.borderaxespad"] = 0.5
    legend_text = {'TRUE':'True','FALSE':'False'}



    xlabel = 'Size'
    row0_ylabel = 'Depth'
    row1_ylabel = 'Max breadth'
    row2_ylabel = 'Virality'




   # General: 
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # LEGEND
    legendhandlelength = .8#1
    legendhandletextpad = .5

    horizontal_space_between_subplots = -.5#1.
    vertical_space_between_subplots = .2

    inset_pad = 1


    # INSET NUMBERS
    inset_fontweight='bold'




    plot_dimension_joint = (3,2)

    fig = plt.figure(figsize=figsize_joint)


    ax00 = plt.subplot2grid(plot_dimension_joint,(0,0),rowspan=1,colspan=1)
    ax01 = plt.subplot2grid(plot_dimension_joint,(0,1),rowspan=1,colspan=1)


    plot_this = 'depth'
    for veracity in ['FALSE','TRUE']:
        ax00.scatter(cascade_scatter[veracity]['size'],cascade_scatter[veracity][plot_this],c=marker_color[veracity],s=markersize,alpha=marker_alpha)
        
        
    ax00.set_xscale('log')


    ax00.set_ylabel(row0_ylabel,labelpad=inset_pad)
    ax00.legend((legend_text['FALSE'],legend_text['TRUE']),loc=6,handlelength=legendhandlelength,handletextpad=legendhandletextpad,frameon=False)

    xlim = ax00.get_xlim()
    ylim = ax00.get_ylim()

    print(xlim)
    ax00.text(h_text(xlim),v_text(ylim),'A',fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)
    ax01.text(h_text(xlim),v_text(ylim),'B',fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    ylim_row0 = ax00.get_ylim()


    for veracity in ['FALSE','TRUE']:    
        ax01.plot(cascade_lines[veracity]['size'],cascade_lines[veracity][plot_this]['mean'],ls='-',linewidth=.5,c=marker_color[veracity],alpha=line_alpha)
        
    ax01.set_xscale('log')

    ax01.set_ylim(ylim_row0)

    ax10 = plt.subplot2grid(plot_dimension_joint,(1,0),rowspan=1,colspan=1)
    ax11 = plt.subplot2grid(plot_dimension_joint,(1,1),rowspan=1,colspan=1)

    plot_this = 'max_breadth'
    for veracity in ['FALSE','TRUE']:
        ax10.scatter(cascade_scatter[veracity]['size'],cascade_scatter[veracity][plot_this],c=marker_color[veracity],s=markersize,alpha=marker_alpha)
    ax10.set_xscale('log')
    ax10.set_yscale('log')

    ax10.set_ylabel(row1_ylabel,labelpad=inset_pad)



    xlim = ax10.get_xlim()
    ylim = ax10.get_ylim()

    ax10.text(h_text(xlim),v_text(ylim,logy=True),'C',fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)
    ylim_row1 = ax10.get_ylim()


    for veracity in ['FALSE','TRUE']:    
        ax11.plot(cascade_lines[veracity]['size'],cascade_lines[veracity][plot_this]['mean'],ls='-',linewidth=.5,c=marker_color[veracity],alpha=line_alpha)
        
    ax11.set_xscale('log')
    ax11.set_yscale('log')

    ax11.text(h_text(xlim),v_text(ylim,logy=True),'D',fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    ax11.set_ylim(ylim_row1)



    ax20 = plt.subplot2grid(plot_dimension_joint,(2,0),rowspan=1,colspan=1)
    ax21 = plt.subplot2grid(plot_dimension_joint,(2,1),rowspan=1,colspan=1)



    plot_this = 'virality'
    for veracity in ['FALSE','TRUE']:
        ax20.scatter(cascade_scatter[veracity]['size'],cascade_scatter[veracity][plot_this],c=marker_color[veracity],s=markersize,alpha=marker_alpha)
    ax20.set_xscale('log')


    ax20.set_xlabel(xlabel,labelpad=inset_pad)
    ax20.set_ylabel(row2_ylabel,labelpad=inset_pad)

    xlim = ax20.get_xlim()
    ylim = ax20.get_ylim()
    ax20.text(h_text(xlim),v_text(ylim),'E',fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    ylim_row2 = ax20.get_ylim()


    for veracity in ['FALSE','TRUE']:    
        ax21.plot(cascade_lines[veracity]['size'],cascade_lines[veracity][plot_this]['mean'],ls='-',linewidth=.5,c=marker_color[veracity],alpha=line_alpha)
        
    ax21.set_xscale('log')



    ax21.set_xlabel(xlabel,labelpad=inset_pad)
    ax21.text(h_text(xlim),v_text(ylim),'F',fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)


    ax21.set_ylim(ylim_row2)


    fig.align_ylabels([ax00,ax10,ax20])
    fig.align_xlabels([ax20,ax21])


    plt.tight_layout(h_pad=vertical_space_between_subplots,w_pad = 1)
    for destination in destinations :
        plt.savefig(destination+'FigFalseTrue_joint.png',dpi=400)


def rootedtrees(plot_vecs,keys,colors,p_combined,p_arr,model_name,destinations) :

    destination = 'Outputs/network_drawings/'
    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(17.5),cm_to_inch(6)*1.1)

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,4)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    EB_elinewidth = 15

    zoom_level = 0.035

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    CCDF_yticks = [0.0001,0.001,0.01,0.1]

    time_yticks = [10,1000,100000]

    col0_xticks = [1,10,100,1000,10000]
    col1_xticks = col0_xticks
    col2_xticks = [1,10]#[1,10,100]
    col4_xticks = [1,10,100]#[1,10,100]


    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    ymin = -.1
    ymax = 1.1

    # LABELS
    col0_xlabel = 'Cascade Size'
    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Virality'
    col3_xlabel = 'Unique Users Reached'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'

    # LEGEND
    legendhandlelength = .8#1
    legendhandletextpad = .5

    horizontal_space_between_subplots = -.5#1.
    vertical_space_between_subplots = .2

    inset_pad = 1


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr


    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)

    labels = [r'$R_0=$'+'%s'%keys[0],r'$R_0=$'+'%s'%keys[1]]
    if (keys[0] == 'TRUE') :
        labels = ['True','False']

    size = 3
    ax00.errorbar(plot_vecs[keys[0]][size][0],plot_vecs[keys[0]][size][1],yerr=plot_vecs[keys[0]][size][4],elinewidth=EB_elinewidth,color=colors[0],marker='',markersize=.5,linestyle='',alpha=0.5,label=labels[0])
    ax00.errorbar(plot_vecs[keys[1]][size][0],plot_vecs[keys[1]][size][1],yerr=plot_vecs[keys[1]][size][4],elinewidth=EB_elinewidth,color=colors[1],marker='',markersize = .5,linestyle='',alpha=0.5,label=labels[1])

    ax00.set_ylim(ymin,ymax)
    ax00.set_xticks(np.arange(3,5,1))
    ax00.set_xlim(2,5)
    ax00.text(0.9, 0.9,'Size 3'+ '\n$p_{same}=%.4f$'%p_arr[0],ha='right',va='top',transform=ax00.transAxes)
    ax00.set_ylabel('Probability')

    ax00.legend(loc=3,frameon=False)

    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    size = 4
    ax01.errorbar(plot_vecs[keys[0]][size][0],plot_vecs[keys[0]][size][1],yerr=plot_vecs[keys[0]][size][4],elinewidth=EB_elinewidth,color=colors[0],marker='',linestyle='',alpha=0.5)
    ax01.errorbar(plot_vecs[keys[1]][size][0],plot_vecs[keys[1]][size][1],yerr=plot_vecs[keys[1]][size][4],elinewidth=EB_elinewidth,color=colors[1],marker='',linestyle='',alpha=0.5)

    ax01.set_yticklabels([])
    ax01.set_ylim(ymin,ymax)
    ax01.set_xlim(4,9)
    ax01.set_xticks(np.arange(5,9,1))
    ax01.text(0.9, 0.9,'Size 4'+ '\n$p_{same}=%.4f$'%p_arr[1],ha='right',va='top',transform=ax01.transAxes)




    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=2)
    size = 5
    ax02.errorbar(plot_vecs[keys[0]][size][0],plot_vecs[keys[0]][size][1],yerr=plot_vecs[keys[0]][size][4],elinewidth=EB_elinewidth,color=colors[0],marker='',linestyle='',alpha=0.5)
    ax02.errorbar(plot_vecs[keys[1]][size][0],plot_vecs[keys[1]][size][1],yerr=plot_vecs[keys[1]][size][4],elinewidth=EB_elinewidth,color=colors[1],marker='',linestyle='',alpha=0.5)

    ax02.set_yticklabels([])
    ax02.set_ylim(ymin,ymax)
    ax02.set_xlim(8,18)
    ax02.set_xticks(np.arange(9,18,1))
    ax02.text(0.95, 0.90,'Size 5'+ '\n$p_{same}=%.4f$'%p_arr[2],ha='right',va='top',transform=ax02.transAxes)




    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=4)
    size = 6
    ax10.errorbar(plot_vecs[keys[0]][size][0],plot_vecs[keys[0]][size][1],yerr=plot_vecs[keys[0]][size][4],elinewidth=EB_elinewidth,color=colors[0],marker='',linestyle='',alpha=0.5)
    ax10.errorbar(plot_vecs[keys[1]][size][0],plot_vecs[keys[1]][size][1],yerr=plot_vecs[keys[1]][size][4],elinewidth=EB_elinewidth,color=colors[1],marker='',linestyle='',alpha=0.5)

    ax10.set_ylim(ymin,ymax)
    ax10.set_xlim(17,38)
    ax10.set_xticks(np.arange(18,38,1))
    ax10.text(0.975, 0.9,'Size 6'+ '\n$p_{same}=%.4f$'%p_arr[3]+ '\n$p_{combined}=%.4f$'%p_combined,ha='right',va='top',transform=ax10.transAxes)
    ax10.set_xlabel(' ')
    ax10.set_ylabel('Probability')


    # -----
    # Graph ticks..

    for graph_num in range (3,5) :

        # Annotate the 2nd position with an image
        arr_img = plt.imread(destination+'%i.png'%graph_num, format='png')

        imagebox = OffsetImage(arr_img, zoom=zoom_level)
        imagebox.image.axes = ax00

        ab = AnnotationBbox(imagebox, (graph_num,0),
                            xybox=(0, -7),
                            xycoords=("data", "axes fraction"),
                            boxcoords="offset points",
                            box_alignment=(.5, 1),
                            bboxprops={"edgecolor" : "none"})

        ax00.add_artist(ab)


    for graph_num in range (5,9) :

        # Annotate the 2nd position with an image
        arr_img = plt.imread(destination+'%i.png'%graph_num, format='png')

        imagebox = OffsetImage(arr_img, zoom=zoom_level)
        imagebox.image.axes = ax01

        ab = AnnotationBbox(imagebox, (graph_num,0),
                            xybox=(0, -7),
                            xycoords=("data", "axes fraction"),
                            boxcoords="offset points",
                            box_alignment=(.5, 1),
                            bboxprops={"edgecolor" : "none"})

        ax01.add_artist(ab)    

    for graph_num in range (9,18) :

        # Annotate the 2nd position with an image
        arr_img = plt.imread(destination+'%i.png'%graph_num, format='png')

        imagebox = OffsetImage(arr_img, zoom=zoom_level)
        imagebox.image.axes = ax02

        ab = AnnotationBbox(imagebox, (graph_num,0),
                            xybox=(0, -7),
                            xycoords=("data", "axes fraction"),
                            boxcoords="offset points",
                            box_alignment=(.5, 1),
                            bboxprops={"edgecolor" : "none"})

        ax02.add_artist(ab)        

    for graph_num in range (18,38) :

        # Annotate the 2nd position with an image
        arr_img = plt.imread(destination+'%i.png'%graph_num, format='png')

        imagebox = OffsetImage(arr_img, zoom=zoom_level)
        imagebox.image.axes = ax10

        ab = AnnotationBbox(imagebox, (graph_num,0),
                            xybox=(0, -6),
                            xycoords=("data", "axes fraction"),
                            boxcoords="offset points",
                            box_alignment=(.5, 1),
                            bboxprops={"edgecolor" : "none"})

        ax10.add_artist(ab)          

    plt.tight_layout()
    fig.subplots_adjust(hspace=.5)

    for destinationx in destinations : 
        plt.savefig(destinationx + model_name+'cascade_counts.png',dpi=400)

def rootedtrees_sharetick(plot_vecs,keys,colors,p_combined,p_arr,model_name,destinations) :

    destination = 'Outputs/network_drawings/'
    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(17.5),cm_to_inch(6)*1.1)

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,4)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    EB_elinewidth = 15/10*4

    zoom_level = 0.035

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    CCDF_yticks = [0.0001,0.001,0.01,0.1]

    time_yticks = [10,1000,100000]

    col0_xticks = [1,10,100,1000,10000]
    col1_xticks = col0_xticks
    col2_xticks = [1,10]#[1,10,100]
    col4_xticks = [1,10,100]#[1,10,100]


    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    ymin = -.1
    ymax = 1.1

    # LABELS
    col0_xlabel = 'Cascade Size'
    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Virality'
    col3_xlabel = 'Unique Users Reached'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'

    # LEGEND
    legendhandlelength = .8#1
    legendhandletextpad = .5

    horizontal_space_between_subplots = -.5#1.
    vertical_space_between_subplots = .2

    inset_pad = 1

    delta_tick = 0.16

    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr


    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)

    labels = [r'$R_0=$'+'%s'%keys[0],r'$R_0=$'+'%s'%keys[1]]
    if (keys[0] == 'TRUE') :
        labels = ['True','False']

    size = 3
    ax00.errorbar(plot_vecs[keys[0]][size][0]-np.ones(len(plot_vecs[keys[0]][size][0]))*(delta_tick-0.06),plot_vecs[keys[0]][size][1],yerr=plot_vecs[keys[0]][size][4],elinewidth=EB_elinewidth,color=colors[0],marker='',markersize=.5,linestyle='',alpha=0.5,label=labels[0])
    ax00.errorbar(plot_vecs[keys[1]][size][0]+np.ones(len(plot_vecs[keys[0]][size][0]))*(delta_tick-0.06),plot_vecs[keys[1]][size][1],yerr=plot_vecs[keys[1]][size][4],elinewidth=EB_elinewidth,color=colors[1],marker='',markersize = .5,linestyle='',alpha=0.5,label=labels[1])

    ax00.set_ylim(ymin,ymax)
    ax00.set_xticks(np.arange(3,5,1))
    ax00.set_xlim(2,5)
    ax00.text(0.9, 0.9,'Size 3'+ '\n$p_{same}=%.4f$'%p_arr[0],ha='right',va='top',transform=ax00.transAxes)
    ax00.set_ylabel('Probability')

    ax00.legend(loc=3,frameon=False)

    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    size = 4
    ax01.errorbar(plot_vecs[keys[0]][size][0]-np.ones(len(plot_vecs[keys[0]][size][0]))*delta_tick,plot_vecs[keys[0]][size][1],yerr=plot_vecs[keys[0]][size][4],elinewidth=EB_elinewidth,color=colors[0],marker='',linestyle='',alpha=0.5)
    ax01.errorbar(plot_vecs[keys[1]][size][0]+np.ones(len(plot_vecs[keys[0]][size][0]))*delta_tick,plot_vecs[keys[1]][size][1],yerr=plot_vecs[keys[1]][size][4],elinewidth=EB_elinewidth,color=colors[1],marker='',linestyle='',alpha=0.5)

    ax01.set_yticklabels([])
    ax01.set_ylim(ymin,ymax)
    ax01.set_xlim(4,9)
    ax01.set_xticks(np.arange(5,9,1))
    ax01.text(0.9, 0.9,'Size 4'+ '\n$p_{same}=%.4f$'%p_arr[1],ha='right',va='top',transform=ax01.transAxes)




    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=2)
    size = 5
    ax02.errorbar(plot_vecs[keys[0]][size][0]-np.ones(len(plot_vecs[keys[0]][size][0]))*delta_tick,plot_vecs[keys[0]][size][1],yerr=plot_vecs[keys[0]][size][4],elinewidth=EB_elinewidth,color=colors[0],marker='',linestyle='',alpha=0.5)
    ax02.errorbar(plot_vecs[keys[1]][size][0]+np.ones(len(plot_vecs[keys[0]][size][0]))*delta_tick,plot_vecs[keys[1]][size][1],yerr=plot_vecs[keys[1]][size][4],elinewidth=EB_elinewidth,color=colors[1],marker='',linestyle='',alpha=0.5)

    ax02.set_yticklabels([])
    ax02.set_ylim(ymin,ymax)
    ax02.set_xlim(8,18)
    ax02.set_xticks(np.arange(9,18,1))
    ax02.text(0.95, 0.90,'Size 5'+ '\n$p_{same}=%.4f$'%p_arr[2],ha='right',va='top',transform=ax02.transAxes)




    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=4)
    size = 6
    ax10.errorbar(plot_vecs[keys[0]][size][0]-np.ones(len(plot_vecs[keys[0]][size][0]))*delta_tick,plot_vecs[keys[0]][size][1],yerr=plot_vecs[keys[0]][size][4],elinewidth=EB_elinewidth,color=colors[0],marker='',linestyle='',alpha=0.5)
    ax10.errorbar(plot_vecs[keys[1]][size][0]+np.ones(len(plot_vecs[keys[0]][size][0]))*delta_tick,plot_vecs[keys[1]][size][1],yerr=plot_vecs[keys[1]][size][4],elinewidth=EB_elinewidth,color=colors[1],marker='',linestyle='',alpha=0.5)

    ax10.set_ylim(ymin,ymax)
    ax10.set_xlim(17,38)
    ax10.set_xticks(np.arange(18,38,1))
    ax10.text(0.975, 0.9,'Size 6'+ '\n$p_{same}=%.4f$'%p_arr[3]+ '\n$p_{combined}=%.4f$'%p_combined,ha='right',va='top',transform=ax10.transAxes)
    ax10.set_xlabel(' ')
    ax10.set_ylabel('Probability')


    # -----
    # Graph ticks..

    for graph_num in range (3,5) :

        # Annotate the 2nd position with an image
        arr_img = plt.imread(destination+'%i.png'%graph_num, format='png')

        imagebox = OffsetImage(arr_img, zoom=zoom_level)
        imagebox.image.axes = ax00

        ab = AnnotationBbox(imagebox, (graph_num,0),
                            xybox=(0, -7),
                            xycoords=("data", "axes fraction"),
                            boxcoords="offset points",
                            box_alignment=(.5, 1),
                            bboxprops={"edgecolor" : "none"})

        ax00.add_artist(ab)


    for graph_num in range (5,9) :

        # Annotate the 2nd position with an image
        arr_img = plt.imread(destination+'%i.png'%graph_num, format='png')

        imagebox = OffsetImage(arr_img, zoom=zoom_level)
        imagebox.image.axes = ax01

        ab = AnnotationBbox(imagebox, (graph_num,0),
                            xybox=(0, -7),
                            xycoords=("data", "axes fraction"),
                            boxcoords="offset points",
                            box_alignment=(.5, 1),
                            bboxprops={"edgecolor" : "none"})

        ax01.add_artist(ab)    

    for graph_num in range (9,18) :

        # Annotate the 2nd position with an image
        arr_img = plt.imread(destination+'%i.png'%graph_num, format='png')

        imagebox = OffsetImage(arr_img, zoom=zoom_level)
        imagebox.image.axes = ax02

        ab = AnnotationBbox(imagebox, (graph_num,0),
                            xybox=(0, -7),
                            xycoords=("data", "axes fraction"),
                            boxcoords="offset points",
                            box_alignment=(.5, 1),
                            bboxprops={"edgecolor" : "none"})

        ax02.add_artist(ab)        

    for graph_num in range (18,38) :

        # Annotate the 2nd position with an image
        arr_img = plt.imread(destination+'%i.png'%graph_num, format='png')

        imagebox = OffsetImage(arr_img, zoom=zoom_level)
        imagebox.image.axes = ax10

        ab = AnnotationBbox(imagebox, (graph_num,0),
                            xybox=(0, -6),
                            xycoords=("data", "axes fraction"),
                            boxcoords="offset points",
                            box_alignment=(.5, 1),
                            bboxprops={"edgecolor" : "none"})

        ax10.add_artist(ab)          

    plt.tight_layout()
    fig.subplots_adjust(hspace=.5)

    for destinationx in destinations : 
        plt.savefig(destinationx + model_name+'cascade_counts_sharetick.png',dpi=400)        

# TRUE / FALSE NEWS
def FalseTrue_degrees(offsprings,offsprings_examples,destinations) :



    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(17.8/6*2.3),cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,2)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    CCDF_yticks = [0.000001,0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = [1,10,100,1000,10000]
    col1_xticks = col0_xticks
    col2_xticks = [1,10,100,1000,10000]#[1,10,100]
    coldepth_xticks = [1,10,100,1000,10000]#[1,10,100]


    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxx_col0 = 40000
    maxx_col1 = 100000
    maxx_col2 = 14#100
    maxx_depth = 40000#100

    #maxx_col3 = 70000

    maxy_CCDF = 1.3
    miny_CCDF = 0.0000001
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    coldepth_xlabel = 'Cascade Depth'

    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Virality'
    col3_xlabel = 'Unique Users'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()


    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'


    ax00.hist(offsprings['TRUE'],bins=max(offsprings['TRUE']),density=True,facecolor='g',alpha=0.50,label='True')
    ax00.hist(offsprings['FALSE'],bins=max(offsprings['FALSE']),density=True,facecolor='r',alpha=0.50, label='False')    

    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([1,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    ax00.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        
    plt.text(0.04, 0.82, 'Density',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure,verticalalignment='center')

    ax00.legend(loc=1,handlelength=legendhandlelength,frameon=False)

    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    # (0,1) Depth
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(offsprings['TRUE'])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(offsprings['FALSE'])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_depth])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col2_xticks)

    ax01.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
        

    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    
    
    plt.text(0.525, 0.82, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure,verticalalignment='center')

    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'


    
    ax10.hist(offsprings_examples['TRUE'],bins=max(offsprings_examples['TRUE']),density=True,facecolor='g',alpha=0.50,label='True')
    ax10.hist(offsprings_examples['FALSE'],bins=max(offsprings_examples['FALSE']),density=True,facecolor='r',alpha=0.50, label='False')    
    
    
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([1,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)
    ax10.set_ylabel(CCDF_ylabel)

    ax10.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.04, 0.37, 'Density',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure,verticalalignment='center')


    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(offsprings_examples['TRUE'])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(offsprings_examples['FALSE'])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_depth])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col2_xticks)

    ax11.set_xlabel(coldepth_xlabel)

    ax11.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')
    ax11.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.525, 0.37, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure,verticalalignment='center')



    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()

    for destination in destinations :
        plt.savefig(destination+'FigDegrees.png',dpi=400)

# TRUE / FALSE NEWS
def Plot_degrees(offsprings,offsprings_examples,destinations) :



    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(17.8/6*3.3),cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,3)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    CCDF_yticks = [0.000001,0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = [1,10,100,1000,10000]
    col1_xticks = [1,10]
    col2_xticks = [1,10]
    coldepth_xticks = [1,10]

    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxx_col0 = 40000
    maxx_col1 = 20
    maxx_col2 = 20
    maxx_depth = 20

    maxy_CCDF = 1.3
    miny_CCDF = 0.0000001
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Out-degree+1'
    coldepth_xlabel = 'Out-degree+1'

    col1_xlabel = 'Out-degree+1'
    col2_xlabel = 'Out-degree+1'
    col3_xlabel = 'Out-degree+1'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()


    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'


    plotthis = sorted(offsprings['TRUE'])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth,label='True')


    plotthis = sorted(offsprings['FALSE'])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth,label='False')    
    
    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([1,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    ax00.set_ylabel(CCDF_ylabel)


    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        
    plt.text(0.033, 0.77, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)
    
    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False)

    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    # (0,1) Depth
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(offsprings['SIR80'])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth,label=r'$R_0=0.8$')


    plotthis = sorted(offsprings['SIR90'])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C1',linewidth=FT_linewidth,label=r'$R_0=0.9$')
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_depth])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col2_xticks)

    ax01.set_ylabel(CCDF_ylabel)


    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
        
    ax01.legend(loc=3,handlelength=legendhandlelength,frameon=False)

        
    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    
    
    
    # (0,2) Depth
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(offsprings['SIR80'])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth,label='SIR')


    plotthis = sorted(offsprings['IC80'])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth,label='IC')
    ax02.set_xscale('log')
    ax02.set_yscale('log')
    ax02.set_xlim([1,maxx_depth])
    ax02.set_ylim([miny_CCDF,maxy_CCDF])

    ax02.set_yticks(CCDF_yticks)
    ax02.set_xticks(col2_xticks)

    ax02.set_ylabel(CCDF_ylabel)


    if (hide_CCDF_yticklabels==True) :
        ax02.set_yticklabels([])
        ax02.set_ylabel('')
        ax02.set_xlabel('')
        
    ax02.legend(loc=3,handlelength=legendhandlelength,frameon=False)

    ax02.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    
        
    

    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'


    plotthis = sorted(offsprings_examples['TRUE'])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='g',linewidth=FT_linewidth)


    plotthis = sorted(offsprings_examples['FALSE'])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)    
    
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([1,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)
    ax10.set_ylabel(CCDF_ylabel)

    ax10.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  


    plt.text(0.033, 0.29+0.03, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)


    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(offsprings_examples['SIR80'])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth)


    plotthis = sorted(offsprings_examples['SIR90'])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C1',linewidth=FT_linewidth)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_depth])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col2_xticks)

    ax11.set_xlabel(coldepth_xlabel)

    ax11.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')
    ax11.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  


  
    # (1,2)
    ax12 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(offsprings_examples['SIR_SIR80vsIC80'])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth)


    plotthis = sorted(offsprings_examples['IC_SIR80vsIC80'])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth)
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_xlim([1,maxx_depth])
    ax12.set_ylim([miny_CCDF,maxy_CCDF])
    ax12.set_yticks(CCDF_yticks)
    ax12.set_xticks(col2_xticks)

    ax12.set_xlabel(coldepth_xlabel)

    ax12.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax12.set_yticklabels([])
        ax12.set_ylabel('')
    ax12.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  
    
    
    
    

    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()

    for destination in destinations :
        plt.savefig(destination+'FigDegrees_3col.png',dpi=400)                
        

# Political / other false
def PoliticalOther_ccdfs(topology_Other,topology_Politics,example_cascades_categories,median_timesingle_categories,median_time_categories,time05_categories,time95_categories,destinations,plotcategories=['Other','Politics']) :

    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(14),cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,4)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    CCDF_yticks = [0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = [1,10,100,1000,10000]
    col1_xticks = col0_xticks
    col2_xticks = [1,10]#[1,10,100]

    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxx_col0 = 100000
    maxx_col1 = 100000
    maxx_col2 = 14

    maxy_CCDF = 1.3
    miny_CCDF = 0.000007
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Depth'#'Cascade Virality'
    col3_xlabel = 'Unique Users'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()

    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(topology_Other['%s'%plotcategories[0]]['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth,label=plotcategories[0])


    plotthis = sorted(topology_Politics['%s'%plotcategories[1]]['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth,label=plotcategories[1])


    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([1,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        

    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False)


    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    plt.text(0.03620, 0.5+0.275, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)


    # (0,1) BREADTH
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(topology_Other['%s'%plotcategories[0]]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(topology_Politics['%s'%plotcategories[1]]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_col1])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col1_xticks)

    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
        

    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    

    # (0,2) VIRALITY
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(topology_Other['%s'%plotcategories[0]]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(topology_Politics['%s'%plotcategories[1]]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax02.set_xscale('log')
    ax02.set_yscale('log')
    ax02.set_xlim([1,maxx_col2])
    ax02.set_ylim([miny_CCDF,maxy_CCDF])
    ax02.set_yticks(CCDF_yticks)
    ax02.set_xticks(col2_xticks)

    if (hide_CCDF_yticklabels==True) :
        ax02.set_yticklabels([])
        ax02.set_ylabel('')
        ax02.set_xlabel('')
        
        
    ax02.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # (0,3) TIME
    ax03 = plt.subplot2grid(plot_dimension,(0,3),rowspan=1,colspan=1)

    ax03.plot(median_timesingle_categories['Othersizes'],median_timesingle_categories['Other'],'b',linewidth=FT_linewidth)
    ax03.plot(median_timesingle_categories['Politicssizes'],median_timesingle_categories['Politics'],'r',linewidth=FT_linewidth)

    ax03.set_yscale('log')
    ax03.set_ylim([miny_time,maxy_time])

    ax03.set_xticks(ax03_xticks)
    ax03.set_yticks(time_yticks)

    ax03.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  


    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(example_cascades_categories['%s'%plotcategories[0]]['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_categories['%s'%plotcategories[1]]['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([1,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)

    ax10.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.03620, 0.30, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)
    # INSERT SMALL IMAGE
    im_ax10 = inset_axes(ax10,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax10.imshow(im_size)
    im_ax10.axis('off')




    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_categories['%s'%plotcategories[0]]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_categories['%s'%plotcategories[1]]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_col1])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col1_xticks)

    ax11.set_xlabel(col1_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')
    ax11.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax11 = inset_axes(ax11,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax11.imshow(im_max_breadth)
    im_ax11.axis('off')



    # (1,2)
    ax12 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(example_cascades_categories['%s'%plotcategories[0]]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_categories['%s'%plotcategories[1]]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_xlim([1,maxx_col2])
    ax12.set_ylim([miny_CCDF,maxy_CCDF])
    ax12.set_yticks(CCDF_yticks)
    ax12.set_xticks(col2_xticks)


    ax12.set_xlabel(col2_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax12.set_yticklabels([])
        ax12.set_ylabel('')

            
    ax12.set_ylabel('G',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight) 

    # INSERT SMALL IMAGE
    im_ax12 = inset_axes(ax12,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax12.imshow(im_depth)
    im_ax12.axis('off')




    # (1,3)
    ax13 = plt.subplot2grid(plot_dimension,(1,3),rowspan=1,colspan=1)

    ax13.plot(median_time_categories['sizes'],median_time_categories[plotcategories[0]],color='b',linewidth=FT_linewidth)
    ax13.plot(median_time_categories['sizes'],median_time_categories[plotcategories[1]],color='r',linewidth=FT_linewidth)
    ax13.fill_between(median_time_categories['sizes'],time05_categories[plotcategories[1]],time95_categories[plotcategories[1]],color='r',alpha=alpha90)

    ax13.set_yscale('log')
    ax13.set_ylim([miny_time,maxy_time])
    ax13.set_xticks(ax13_xticks)
    ax13.set_yticks(time_yticks)

    ax13.set_xlabel(col3_xlabel)
    ax13.set_ylabel(col3_ylabel)

    ax13.set_ylabel('H',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  
    # INSERT SMALL IMAGE
    im_ax13 = inset_axes(ax13,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax13.imshow(im_virality)
    im_ax13.axis('off')

    # INCLUDE SIR




    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()

    for destination in destinations :
        plt.savefig(destination+'Fig'+'%s%s'%(plotcategories[1],plotcategories[0])+'.png',dpi=400)

# Political / other false
def PoliticalOther_5ccdfs(topology_Other,topology_Politics,example_cascades_categories,median_timesingle_categories,median_time_categories,time05_categories,time95_categories,destinations,plot_categories = ['Other','Politics']) :

    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(14)/4.*5,cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,5)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS

    CCDF_yticks = [0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = [1,10,100,1000,10000]
    col1_xticks = col0_xticks
    col2_xticks = [1,10]#[1,10,100]
    col3_xticks = [1,10]#[1,10,100]

    ax04_xticks = [0,20000,40000]
    ax14_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxx_col0 = 100000
    maxx_col1 = 100000
    maxx_col2 = 32
    maxx_col3 = 14



    maxy_CCDF = 1.3
    miny_CCDF = 0.000007
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Depth'
    col3_xlabel = 'Cascade Virality'

    col4_xlabel = 'Unique Users'

    CCDF_ylabel = 'CCDF'
    col4_ylabel = 'Minutes'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()

    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(topology_Other['%s'%plot_categories[0]]['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth,label=plot_categories[0])


    plotthis = sorted(topology_Politics['%s'%plot_categories[1]]['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth,label=plot_categories[1])


    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([1,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    ax00.set_xlabel(col0_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        

    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False)


    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    plt.text(0.03620, 0.5+0.275, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)


    # (0,1) BREADTH
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(topology_Other['%s'%plot_categories[0]]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(topology_Politics['%s'%plot_categories[1]]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_col1])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col1_xticks)

    ax01.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
        

    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    


    # (0,2) Depth
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(topology_Other['%s'%plot_categories[0]]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(topology_Politics['%s'%plot_categories[1]]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax02.set_xscale('log')
    ax02.set_yscale('log')
    ax02.set_xlim([1,maxx_col2])
    ax02.set_ylim([miny_CCDF,maxy_CCDF])
    ax02.set_yticks(CCDF_yticks)
    ax02.set_xticks(col2_xticks)

    ax02.set_xlabel(col2_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax02.set_yticklabels([])
        ax02.set_ylabel('')
        ax02.set_xlabel('')
        
        
    ax02.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  





    # (0,3) VIRALITY
    ax03 = plt.subplot2grid(plot_dimension,(0,3),rowspan=1,colspan=1)
    feature = 'virality'

    plotthis = sorted(topology_Other['%s'%plot_categories[0]]['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(topology_Politics['%s'%plot_categories[1]]['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax03.set_xscale('log')
    ax03.set_yscale('log')
    ax03.set_xlim([1,maxx_col3])
    ax03.set_ylim([miny_CCDF,maxy_CCDF])
    ax03.set_yticks(CCDF_yticks)
    ax03.set_xticks(col3_xticks)

    ax03.set_xlabel(col3_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax03.set_yticklabels([])
        ax03.set_ylabel('')
        ax03.set_xlabel('')
        
        
    ax03.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  


    # (0,4) TIME
    ax04 = plt.subplot2grid(plot_dimension,(0,4),rowspan=1,colspan=1)

    ax04.plot(median_timesingle_categories['%s'%plot_categories[0]+'sizes'],median_timesingle_categories['%s'%plot_categories[0]],'b',linewidth=FT_linewidth)
    ax04.plot(median_timesingle_categories['%s'%plot_categories[1]+'sizes'],median_timesingle_categories['%s'%plot_categories[1]],'r',linewidth=FT_linewidth)

    ax04.set_yscale('log')
    ax04.set_ylim([miny_time,maxy_time])

    ax04.set_xticks(ax04_xticks)
    ax04.set_yticks(time_yticks)

    ax04.set_xlabel(col4_xlabel)

    ax04.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  


    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(example_cascades_categories['%s'%plot_categories[0]]['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_categories['%s'%plot_categories[1]]['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([1,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)

    ax10.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.03620, 0.30, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)
    # INSERT SMALL IMAGE
    im_ax10 = inset_axes(ax10,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax10.imshow(im_size)
    im_ax10.axis('off')




    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_categories['%s'%plot_categories[0]]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_categories['%s'%plot_categories[1]]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_col1])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col1_xticks)

    ax11.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')
    ax11.set_ylabel('G',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax11 = inset_axes(ax11,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax11.imshow(im_max_breadth)
    im_ax11.axis('off')



    # (1,2)
    ax12 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(example_cascades_categories['%s'%plot_categories[0]]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_categories['%s'%plot_categories[1]]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_xlim([1,maxx_col2])
    ax12.set_ylim([miny_CCDF,maxy_CCDF])
    ax12.set_yticks(CCDF_yticks)
    ax12.set_xticks(col2_xticks)


    ax12.set_xlabel(col2_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax12.set_yticklabels([])
        ax12.set_ylabel('')

            
    ax12.set_ylabel('H',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight) 

    # INSERT SMALL IMAGE
    im_ax12 = inset_axes(ax12,
                        height="30%", # set height
                        width="30%", # and width
                        loc=3) # center, you can check the different codes in plt.legend?
    im_ax12.imshow(im_depth)
    im_ax12.axis('off')

    # (1,2)
    ax13 = plt.subplot2grid(plot_dimension,(1,3),rowspan=1,colspan=1)
    feature = 'virality'

    plotthis = sorted(example_cascades_categories['%s'%plot_categories[0]]['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='b',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_categories['%s'%plot_categories[1]]['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='r',linewidth=FT_linewidth)
    ax13.set_xscale('log')
    ax13.set_yscale('log')
    ax13.set_xlim([1,maxx_col3])
    ax13.set_ylim([miny_CCDF,maxy_CCDF])
    ax13.set_yticks(CCDF_yticks)
    ax13.set_xticks(col3_xticks)


    ax13.set_xlabel(col3_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax13.set_yticklabels([])
        ax13.set_ylabel('')

            
    ax13.set_ylabel('I',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight) 

    # INSERT SMALL IMAGE
    im_ax13 = inset_axes(ax13,
                        height="30%", 
                        width="30%", 
                        loc=3) 
    im_ax13.imshow(im_virality)
    im_ax13.axis('off')


    # (1,4)
    ax14 = plt.subplot2grid(plot_dimension,(1,4),rowspan=1,colspan=1)

    ax14.plot(median_time_categories['sizes'],median_time_categories[plot_categories[0]],color='b',linewidth=FT_linewidth)
    ax14.plot(median_time_categories['sizes'],median_time_categories['%s'%plot_categories[1]],color='r',linewidth=FT_linewidth)
    ax14.fill_between(median_time_categories['sizes'],time05_categories['%s'%plot_categories[1]],time95_categories['%s'%plot_categories[1]],color='r',alpha=alpha90)

    ax14.set_yscale('log')
    ax14.set_ylim([miny_time,maxy_time])
    ax14.set_xticks(ax14_xticks)
    ax14.set_yticks(time_yticks)

    ax14.set_xlabel(col4_xlabel)
    ax14.set_ylabel(col4_ylabel)

    ax14.set_ylabel('J',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  
    # INSERT SMALL IMAGE
    im_ax14 = inset_axes(ax14,
                        height="30%", 
                        width="30%", 
                        loc=3) 
    im_ax14.imshow(im_temporal)
    im_ax14.axis('off')

    # INCLUDE SIR
    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()

    for destination in destinations :
        plt.savefig(destination+'Fig'+'%s%s'%(plot_categories[1],plot_categories[0])+'.png',dpi=400)


# General 4-column CCDF plot

def four_column_CCDF(example_cascades_full_data,example_cascades_subsampled,keys,plot_colors,model_name,fig_name,destinations,school=None) :


    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(14),cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,4)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


    maxx_col0 = int(max(max(example_cascades_full_data[keys[0]]['size'],example_cascades_full_data[keys[1]]['size'])))
    maxx_col0 = maxx_col0+2*10**len(str(maxx_col0))
    
    maxx_col1 = int(max(max(example_cascades_full_data[keys[0]]['max_breadth'],example_cascades_full_data[keys[1]]['max_breadth'])))
    maxx_col1 = maxx_col1+2*10**len(str(maxx_col1))

    maxx_col2 = int(max(max(example_cascades_full_data[keys[0]]['depth'],example_cascades_full_data[keys[1]]['depth'])))
    maxx_col2 = maxx_col2+2*10**len(str(maxx_col2))

    maxx_col3 = int(max(max(example_cascades_full_data[keys[0]]['virality'],example_cascades_full_data[keys[1]]['virality'])))
    maxx_col3 = maxx_col3+2*10**len(str(maxx_col3))        


    # TICKS
    CCDF_yticks = [0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = 10**np.arange(len(str(int(maxx_col0))))
    col1_xticks = 10**np.arange(len(str(int(maxx_col1))))
    col2_xticks = 10**np.arange(len(str(int(maxx_col2))))
    col3_xticks = 10**np.arange(len(str(int(maxx_col3))))

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxy_CCDF = 1.3
    miny_CCDF = 0.000007
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Depth'#'Cascade Virality'
    col3_xlabel = 'Cascade Virality'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()



    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(example_cascades_full_data[keys[0]]['%s'%feature])
    if (isfloat(keys[0])) :
        if (school == None) :
            ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth,label=model_name+', '+r'$R_0=%s$'%(keys[0]))
        else : 
            ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth,label=model_name+', '+r'$R_0^{\rm Tree}=%s$'%(keys[0]))

    else : 
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth,label=keys[0])


    plotthis = sorted(example_cascades_full_data[keys[1]]['%s'%feature])
    if (isfloat(keys[1])) :
        if (school==None) :

            ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth,label=model_name+', '+r'$R_0=%s$'%(keys[1]))
        else : 
            ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth,label=model_name+', '+r'$R_0^{\rm Tree}=%s$'%(keys[1]))

    else : 
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth,label=keys[1])


    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([1,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)


    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        

    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False)


    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)
    plt.text(0.03620, 0.5+0.275, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)
    if (school != None) :
        ax00.text(maxx_col0**(0.95),miny_CCDF**(1-0.95),school,verticalalignment = 'top',horizontalalignment='right')    


    # (0,1) BREADTH
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_full_data[keys[0]]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_full_data[keys[1]]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_col1])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col1_xticks)


    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
        

    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    

    # (0,2) DEPTH
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(example_cascades_full_data[keys[0]]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)

    plotthis = sorted(example_cascades_full_data[keys[1]]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax02.set_xscale('log')
    ax02.set_yscale('log')
    ax02.set_xlim([1,maxx_col2])
    ax02.set_ylim([miny_CCDF,maxy_CCDF])
    ax02.set_yticks(CCDF_yticks)
    ax02.set_xticks(col2_xticks)


    if (hide_CCDF_yticklabels==True) :
        ax02.set_yticklabels([])
        ax02.set_ylabel('')
        ax02.set_xlabel('')
        
        
    ax02.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # (0,3) VIRALITY
    ax03 = plt.subplot2grid(plot_dimension,(0,3),rowspan=1,colspan=1)


    feature = 'virality'

    plotthis = sorted(example_cascades_full_data[keys[0]]['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_full_data[keys[1]]['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax03.set_xscale('log')
    ax03.set_yscale('log')
    ax03.set_xlim([1,maxx_col3])
    ax03.set_ylim([miny_CCDF,maxy_CCDF])
    ax03.set_yticks(CCDF_yticks)
    ax03.set_xticks(col3_xticks)


    if (hide_CCDF_yticklabels==True) :
        ax03.set_yticklabels([])
        ax03.set_ylabel('')
        ax03.set_xlabel('')

    ax03.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  


    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(example_cascades_subsampled[keys[0]]['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_subsampled[keys[1]]['%s'%feature])

    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([1,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)

    ax10.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.03620, 0.30+0.025, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)
    
    # INSERT SMALL IMAGE
    im_ax10 = inset_axes(ax10,
                        height="30%", 
                        width="30%", 
                        loc=3) 
    im_ax10.imshow(im_size)
    im_ax10.axis('off')


    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_subsampled[keys[0]]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_subsampled[keys[1]]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_col1])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col1_xticks)

    ax11.set_xlabel(col1_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')

    ax11.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  
    im_ax11 = inset_axes(ax11,
                        height="30%", 
                        width="30%", 
                        loc=3) 
    im_ax11.imshow(im_max_breadth)
    im_ax11.axis('off')


    # (1,2)
    ax12 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(example_cascades_subsampled[keys[0]]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_subsampled[keys[1]]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_xlim([1,maxx_col2])
    ax12.set_ylim([miny_CCDF,maxy_CCDF])
    ax12.set_yticks(CCDF_yticks)
    ax12.set_xticks(col2_xticks)


    ax12.set_xlabel(col2_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax12.set_yticklabels([])
        ax12.set_ylabel('')

    ax12.set_ylabel('G',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight) 

    im_ax12 = inset_axes(ax12,
                        height="30%",
                        width="30%", 
                        loc=3) 
    im_ax12.imshow(im_depth)
    im_ax12.axis('off')

    # (1,3)
    ax13 = plt.subplot2grid(plot_dimension,(1,3),rowspan=1,colspan=1)
    feature = 'virality'

    plotthis = sorted(example_cascades_subsampled[keys[0]]['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)

    plotthis = sorted(example_cascades_subsampled[keys[1]]['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax13.set_xscale('log')
    ax13.set_yscale('log')
    ax13.set_xlim([1,maxx_col3])
    ax13.set_ylim([miny_CCDF,maxy_CCDF])
    ax13.set_yticks(CCDF_yticks)
    ax13.set_xticks(col3_xticks)


    ax13.set_xlabel(col3_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax13.set_yticklabels([])
        ax13.set_ylabel('')

    ax13.set_ylabel('H',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  
    im_ax13 = inset_axes(ax13,
                        height="30%", 
                        width="30%", 
                        loc=3) 
    im_ax13.imshow(im_virality)
    im_ax13.axis('off')




    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()

    for destination in destinations : 
        plt.savefig(destination+fig_name+model_name+'.png',dpi=400)

def four_column_CCDF_network(IC_networks_sampling,school,keys,plot_colors,model_name,fig_name,destinations) :


    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(14),cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,4)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    CCDF_yticks = [0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = [1,10,100,1000]
    col1_xticks = [1,10,100]
    col2_xticks = [1,10,100]

    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxx_col0 = 4000
    maxx_col1 = 1200
    maxx_col2 = 100

    maxy_CCDF = 1.3
    miny_CCDF = 0.000007

    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Depth'
    col3_xlabel = 'Cascade Virality'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)

    def ccdf_x(arr) :
        return arr
    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()

    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(IC_networks_sampling[school][keys[0]]['example_cascade'][keys[0]]['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth,label=r'$R_0^{\rm tree}=%s$'%keys[0])


    plotthis = sorted(IC_networks_sampling[school][keys[1]]['example_cascade'][keys[1]]['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth,label=r'$R_0^{\rm tree}=%s$'%keys[1])


    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([1,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')


    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False)


    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)
    plt.text(0.03620, 0.5+0.275, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)

    ax00.text(maxx_col0**(0.95),miny_CCDF**(1-0.95),school,verticalalignment = 'top',horizontalalignment='right')    



    # (0,1) BREADTH
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(IC_networks_sampling[school][keys[0]]['example_cascade'][keys[0]]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)

    plotthis = sorted(IC_networks_sampling[school][keys[1]]['example_cascade'][keys[1]]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_col1])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col1_xticks)

    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')



    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    

    # (0,2) VIRALITY
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(IC_networks_sampling[school][keys[0]]['example_cascade'][keys[0]]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)

    plotthis = sorted(IC_networks_sampling[school][keys[1]]['example_cascade'][keys[1]]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax02.set_xscale('log')
    ax02.set_yscale('log')
    ax02.set_xlim([1,maxx_col2])
    ax02.set_ylim([miny_CCDF,maxy_CCDF])
    ax02.set_yticks(CCDF_yticks)
    ax02.set_xticks(col2_xticks)


    if (hide_CCDF_yticklabels==True) :
        ax02.set_yticklabels([])
        ax02.set_ylabel('')
        ax02.set_xlabel('')

    ax02.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # (0,3) TIME
    ax03 = plt.subplot2grid(plot_dimension,(0,3),rowspan=1,colspan=1)

    feature = 'virality'#'virality'

    plotthis = sorted(IC_networks_sampling[school][keys[0]]['example_cascade'][keys[0]]['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)

    plotthis = sorted(IC_networks_sampling[school][keys[1]]['example_cascade'][keys[1]]['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax03.set_xscale('log')
    ax03.set_yscale('log')
    ax03.set_xlim([1,maxx_col2])
    ax03.set_ylim([miny_CCDF,maxy_CCDF])
    ax03.set_yticks(CCDF_yticks)
    ax03.set_xticks(col2_xticks)

    if (hide_CCDF_yticklabels==True) :
        ax03.set_yticklabels([])
        ax03.set_ylabel('')
        ax03.set_xlabel('')

    ax03.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  


    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(IC_networks_sampling[school]['compare']['example_cascade'][keys[0]]['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)

    plotthis = sorted(IC_networks_sampling[school]['compare']['example_cascade'][keys[1]]['%s'%feature])

    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([1,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)

    ax10.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.03620, 0.30, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)
    # INSERT SMALL IMAGE
    im_ax10 = inset_axes(ax10,
                          height="30%", 
                          width="30%", 
                          loc=3) 
    im_ax10.imshow(im_size)
    im_ax10.axis('off')


    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(IC_networks_sampling[school]['compare']['example_cascade'][keys[0]]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)


    plotthis = sorted(IC_networks_sampling[school]['compare']['example_cascade'][keys[1]]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_col1])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col1_xticks)

    ax11.set_xlabel(col1_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')
    ax11.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  
    im_ax11 = inset_axes(ax11,
                          height="30%", 
                          width="30%", 
                          loc=3) 
    im_ax11.imshow(im_max_breadth)
    im_ax11.axis('off')


    # (1,2)
    ax12 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    feature = 'depth'

    plotthis = sorted(IC_networks_sampling[school]['compare']['example_cascade'][keys[0]]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)

    plotthis = sorted(IC_networks_sampling[school]['compare']['example_cascade'][keys[1]]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_xlim([1,maxx_col2])
    ax12.set_ylim([miny_CCDF,maxy_CCDF])
    ax12.set_yticks(CCDF_yticks)
    ax12.set_xticks(col2_xticks)


    ax12.set_xlabel(col2_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax12.set_yticklabels([])
        ax12.set_ylabel('')


    ax12.set_ylabel('G',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight) 

    im_ax12 = inset_axes(ax12,
                          height="30%", 
                          width="30%", 
                          loc=3) 
    im_ax12.imshow(im_depth)
    im_ax12.axis('off')

    # (1,3)
    ax13 = plt.subplot2grid(plot_dimension,(1,3),rowspan=1,colspan=1)
    feature = 'virality'

    plotthis = sorted(IC_networks_sampling[school]['compare']['example_cascade'][keys[0]]['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[0],linewidth=FT_linewidth)


    plotthis = sorted(IC_networks_sampling[school]['compare']['example_cascade'][keys[1]]['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=plot_colors[1],linewidth=FT_linewidth)
    ax13.set_xscale('log')
    ax13.set_yscale('log')
    ax13.set_xlim([1,maxx_col2])
    ax13.set_ylim([miny_CCDF,maxy_CCDF])
    ax13.set_yticks(CCDF_yticks)
    ax13.set_xticks(col2_xticks)


    ax13.set_xlabel(col3_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax13.set_yticklabels([])
        ax13.set_ylabel('')

    ax13.set_ylabel('H',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  
    im_ax13 = inset_axes(ax13,
                          height="30%",
                          width="30%", 
                          loc=3) 
    im_ax13.imshow(im_virality)
    im_ax13.axis('off')

    # INCLUDE SIR
    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()
    for destination in destinations :
        plt.savefig(destination+fig_name+model_name+'network_%s.png'%school,dpi=400)

def six_column_CCDF(example_cascades_SIR08only,example_cascades_SIR09only,example_cascades_IC08only,example_cascades_IC09only,example_cascades_SIR,example_cascades_IC,example_cascades_ICSIR08,destinations) :
    # NOW DO MIXED...


    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(17.8),cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,6)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    CCDF_yticks = [0.0001,0.01,1]

    time_yticks = []

    col0_xticks = [1,10,100,1000]
    col1_xticks = [1,10,100]
    col2_xticks = [1,10,100]

    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = True

    # AXIS LIMS
    maxx_col0 = 3200
    maxx_col1 = 120
    maxx_col2 = 100

    maxy_CCDF = 1.3
    miny_CCDF = 0.000007
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    col1_xlabel = 'Cascade Max-Breadth'
    col2_xlabel = 'Cascade Depth'
    col3_xlabel = 'Cascade Virality'

    CCDF_ylabel = 'CCDF'
    col3_ylabel = 'Minutes'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'


    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()


    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(example_cascades_SIR08only[0.8]['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth,label=r'$R_0=0.8$')


    plotthis = sorted(example_cascades_SIR09only[0.9]['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C1',linewidth=FT_linewidth,label=r'$R_0=0.9$')


    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([1,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    ax00.set_ylabel(CCDF_ylabel)
    ax00.set_xlabel(col0_xlabel)


    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        

    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False)

    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)


    # (0,1) BREADTH
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_SIR08only[0.8]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_SIR09only[0.9]['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C1',linewidth=FT_linewidth)
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_col1])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col1_xticks)

    ax01.set_ylabel(CCDF_ylabel)
    ax01.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
        

    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    

    # IC 
    # SIZE 
    # ----------
    # define
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(example_cascades_IC08only[0.8]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth,label=r'$R_0=0.8$')


    plotthis = sorted(example_cascades_IC09only[0.9]['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C4',linewidth=FT_linewidth,label=r'$R_0=0.9$')


    ax02.set_xscale('log')
    ax02.set_yscale('log')

    ax02.set_xlim([1,maxx_col0])
    ax02.set_ylim([miny_CCDF,maxy_CCDF])

    ax02.set_yticks(CCDF_yticks)
    ax02.set_xticks(col0_xticks)

    ax02.set_xlabel(col0_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax02.set_yticklabels([])
        ax02.set_xlabel(col0_xlabel)
        ax02.set_xlabel('')
        

    ax02.legend(loc=3,handlelength=legendhandlelength,frameon=False)


    ax02.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)


    # (0,1) BREADTH
    ax03 = plt.subplot2grid(plot_dimension,(0,3),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_IC08only[0.8]['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth)

    plotthis = sorted(example_cascades_IC09only[0.9]['%s'%feature])
    ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C4',linewidth=FT_linewidth)
    ax03.set_xscale('log')
    ax03.set_yscale('log')
    ax03.set_xlim([1,maxx_col1])
    ax03.set_ylim([miny_CCDF,maxy_CCDF])

    ax03.set_yticks(CCDF_yticks)
    ax03.set_xticks(col1_xticks)

    ax03.set_ylabel(CCDF_ylabel)
    ax03.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax03.set_yticklabels([])
        ax03.set_ylabel('')
        ax03.set_xlabel('')
         
    ax03.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    

    # MIXED 
    # SIZE 
    # ----------
    # define
    ax04 = plt.subplot2grid(plot_dimension,(0,4),rowspan=1,colspan=1)
    feature = 'size'


    plotthis = sorted(example_cascades_SIR08only[0.8]['%s'%feature])
    ax04.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth,label=r'$SIR$')
    
    plotthis = sorted(example_cascades_IC08only[0.8]['%s'%feature])
    ax04.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth,label=r'$IC$')




    ax04.set_xscale('log')
    ax04.set_yscale('log')

    ax04.set_xlim([1,maxx_col0])
    ax04.set_ylim([miny_CCDF,maxy_CCDF])

    ax04.set_yticks(CCDF_yticks)
    ax04.set_xticks(col0_xticks)

    ax04.set_ylabel(CCDF_ylabel)
    ax04.set_xlabel(col0_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax04.set_yticklabels([])
        ax04.set_xlabel(col0_xlabel)
        ax04.set_xlabel('')
        

    ax04.legend(loc=3,handlelength=legendhandlelength,frameon=False)

    ax04.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)


    # (0,1) BREADTH
    ax05 = plt.subplot2grid(plot_dimension,(0,5),rowspan=1,colspan=1)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_IC08only[0.8]['%s'%feature])
    ax05.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_SIR08only[0.8]['%s'%feature])
    ax05.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth)
    ax05.set_xscale('log')
    ax05.set_yscale('log')
    ax05.set_xlim([1,maxx_col1])
    ax05.set_ylim([miny_CCDF,maxy_CCDF])

    ax05.set_xticks(col1_xticks)


    ax05.set_ylabel(CCDF_ylabel)
    ax05.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax05.set_yticklabels([])
        ax05.set_ylabel('')
        ax05.set_xlabel('')
        
    ax05.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    


    # ----
    # BOTTOM
    # ----

    # SIZE 
    # ----------
    # define
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    ax10.set_facecolor('white')
    ax10.set_zorder(100)    
    
    feature = 'size'

    plotthis = sorted(example_cascades_SIR[0.8]['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth,label=r'$R_0=0.8$')


    plotthis = sorted(example_cascades_SIR[0.9]['%s'%feature])
    ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C1',linewidth=FT_linewidth,label=r'$R_0=0.9$')


    ax10.set_xscale('log')
    ax10.set_yscale('log')

    ax10.set_xlim([1,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])

    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)

    ax10.set_ylabel(CCDF_ylabel)
    ax10.set_xlabel(col0_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax10.set_xlabel(col0_xlabel)

    ax10.set_ylabel('G',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    im_ax10 = inset_axes(ax10,height="30%", width="30%", loc=3)
    im_ax10.imshow(im_size)
    im_ax10.axis('off')
    im_ax10.set_zorder(100)        

    # (0,1) BREADTH
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    ax11.set_facecolor('white')
    ax11.set_zorder(100)    
    
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_SIR[0.8]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth)

    plotthis = sorted(example_cascades_SIR[0.9]['%s'%feature])
    ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C1',linewidth=FT_linewidth)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_col1])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])

    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col1_xticks)

    ax11.set_ylabel(CCDF_ylabel)
    ax11.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')

    ax11.set_ylabel('H',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    

    im_ax11 = inset_axes(ax11, height="30%", width="30%",loc=3)
    im_ax11.imshow(im_max_breadth)
    im_ax11.axis('off')
    im_ax11.set_zorder(100)    


    # IC 
    # SIZE 
    # ----------
    # define
    ax12 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    ax12.set_facecolor('white')    
    ax12.set_zorder(100)        
    feature = 'size'

    plotthis = sorted(example_cascades_IC[0.8]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth,label=r'$R_0=0.8$')


    plotthis = sorted(example_cascades_IC[0.9]['%s'%feature])
    ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C4',linewidth=FT_linewidth,label=r'$R_0=0.8$')


    ax12.set_xscale('log')
    ax12.set_yscale('log')

    ax12.set_xlim([1,maxx_col0])
    ax12.set_ylim([miny_CCDF,maxy_CCDF])

    ax12.set_yticks(CCDF_yticks)
    ax12.set_xticks(col0_xticks)

    ax12.set_xlabel(col0_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax12.set_yticklabels([])
        ax12.set_xlabel(col0_xlabel)

    ax12.set_ylabel('I',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)
    im_ax12 = inset_axes(ax12,height="30%",width="30%", loc=3)
    im_ax12.imshow(im_size)
    im_ax12.axis('off')
    im_ax12.set_zorder(100)        

    # (0,1) BREADTH
    ax13 = plt.subplot2grid(plot_dimension,(1,3),rowspan=1,colspan=1)
    ax13.set_facecolor('white')    
    ax13.set_zorder(100)
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_IC[0.8]['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_IC[0.9]['%s'%feature])
    ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C4',linewidth=FT_linewidth)
    ax13.set_xscale('log')
    ax13.set_yscale('log')
    ax13.set_xlim([1,maxx_col1])
    ax13.set_ylim([miny_CCDF,maxy_CCDF])

    ax13.set_yticks(CCDF_yticks)
    ax13.set_xticks(col1_xticks)

    ax13.set_ylabel(CCDF_ylabel)
    ax13.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax13.set_yticklabels([])
        ax13.set_ylabel('')
         
    ax13.set_ylabel('J',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    

    im_ax13 = inset_axes(ax13,height="30%", width="30%",loc=3) 
    im_ax13.imshow(im_max_breadth)
    im_ax13.axis('off')
    im_ax13.set_zorder(100)    


    # MIXED 
    # SIZE 
    # ----------
    # define
    ax14 = plt.subplot2grid(plot_dimension,(1,4),rowspan=1,colspan=1)
    ax14.set_facecolor('white')    
    ax14.set_zorder(100)    
    feature = 'size'

    plotthis = sorted(example_cascades_ICSIR08['IC']['%s'%feature])
    ax14.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth,label=r'$R_0=0.8$')


    plotthis = sorted(example_cascades_ICSIR08['SIR']['%s'%feature])
    ax14.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth,label=r'$R_0=0.8$')


    ax14.set_xscale('log')
    ax14.set_yscale('log')

    ax14.set_xlim([1,maxx_col0])
    ax14.set_ylim([miny_CCDF,maxy_CCDF])

    ax14.set_yticks(CCDF_yticks)
    ax14.set_xticks(col0_xticks)

    ax14.set_xlabel(col0_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax14.set_yticklabels([])
        
        ax14.set_xlabel(col0_xlabel)

    ax14.set_xlabel(col0_xlabel)

    ax14.set_ylabel('K',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    # INSERT SMALL IMAGE
    im_ax14 = inset_axes(ax14,height="30%", width="30%", loc=3)
    im_ax14.imshow(im_size)
    im_ax14.axis('off')
    im_ax14.set_zorder(100)    

    # (0,1) BREADTH
    ax15 = plt.subplot2grid(plot_dimension,(1,5),rowspan=1,colspan=1)
    ax15.set_facecolor('white')    
    ax15.set_zorder(100)    
    feature = 'max_breadth'

    plotthis = sorted(example_cascades_ICSIR08['IC']['%s'%feature])
    ax15.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C2',linewidth=FT_linewidth)


    plotthis = sorted(example_cascades_ICSIR08['SIR']['%s'%feature])
    ax15.plot(ccdf_x(plotthis),ccdf_y(plotthis),color='C0',linewidth=FT_linewidth)
    ax15.set_xscale('log')
    ax15.set_yscale('log')
    ax15.set_xlim([1,maxx_col1])
    ax15.set_ylim([miny_CCDF,maxy_CCDF])

    ax15.set_yticks(CCDF_yticks)
    ax15.set_yticklabels([])

    ax15.set_xlabel(col1_xlabel)

    ax15.set_xticks(col1_xticks)

    ax15.set_ylabel(CCDF_ylabel)
    ax15.set_xlabel(col1_xlabel)

    if (hide_CCDF_yticklabels==True) :
        ax15.set_yticklabels([])
        ax15.set_ylabel('')
        
  
    ax15.set_ylabel('L',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    

    # INSERT SMALL IMAGE
    im_ax15 = inset_axes(ax15,height="30%", width="30%",loc=3) 
    im_ax15.imshow(im_max_breadth)
    im_ax15.axis('off')
    im_ax15.set_zorder(100)    

    plt.text(0.013, 0.78, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)
    plt.text(0.013, 0.32, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)

    fig.align_ylabels([ax00,ax10])
    plt.tight_layout(h_pad = 0.0)

    
    
    
    
    
    
    
    
    
    # MAKE COLORED BOXES
    # define

    ax00.set_facecolor('white')
    ax00.set_zorder(100)
    
    ax01.set_facecolor('white')
    ax01.set_zorder(100)

    ax02.set_facecolor('white')
    ax02.set_zorder(100)

    ax03.set_facecolor('white')
    ax03.set_zorder(100)

    ax04.set_facecolor('white')
    ax04.set_zorder(100)

    ax05.set_facecolor('white')
    ax05.set_zorder(100)



    separation = 0.0025
    first_border = 0.355
    second_border = first_border - separation/2+ (1-first_border-separation)/2

    boxcolor = 'b'
    boxalpha = 0.06

    import matplotlib.gridspec as gridspec
    outergs = gridspec.GridSpec(1, 1)
    outergs.update(bottom=0.03,left=0+separation, right=first_border-0.5*separation,top=1.05)
    outerax = fig.add_subplot(outergs[0])
    outerax.tick_params(axis='both',which='both',bottom=0,left=0,
                        labelbottom=0, labelleft=0)
    outerax.set_facecolor(boxcolor)
    
    outerax.patch.set_alpha(boxalpha)
    
    outerax.set_zorder(1)




    outergs1 = gridspec.GridSpec(1, 1)
    outergs1.update(bottom=0.03,left=first_border+separation, right=second_border-separation,top=1.05)
    outerax1 = fig.add_subplot(outergs1[0])
    outerax1.tick_params(axis='both',which='both',bottom=0,left=0,
                        labelbottom=0, labelleft=0)
    outerax1.set_facecolor(boxcolor)
    
    outerax1.patch.set_alpha(boxalpha)
    
    outerax1.set_zorder(1)



    outergs2 = gridspec.GridSpec(1, 1)
    outergs2.update(bottom=0.03,left=second_border+separation, right=1.00-3*separation,top=1.05)
    outerax2 = fig.add_subplot(outergs2[0])
    outerax2.tick_params(axis='both',which='both',bottom=0,left=0,
                        labelbottom=0, labelleft=0)
    outerax2.set_facecolor(boxcolor)
    
    outerax2.patch.set_alpha(boxalpha)
    
    outerax2.set_zorder(1)

    outerax.text(0.50,0.95,'SIR vs. SIR',color='k',alpha=0.60,fontsize=12,horizontalalignment='center',verticalalignment='center')

    outerax1.text(0.50,0.95,'IC vs. IC',color='k',alpha=0.60,fontsize=12,horizontalalignment='center',verticalalignment='center')

    outerax2.text(0.50,0.95,'SIR vs. IC',color='k',alpha=0.60,fontsize=12,horizontalalignment='center',verticalalignment='center')


    
    
    
    
    
    
    
    
    
    
    for destination in destinations : 
        plt.savefig(destination+'FigSIRIC_mix_colors.png',dpi=400,bbox_inches = 'tight')
        plt.savefig(destination+'FigSIRIC_mix_colors.svg',dpi=400,bbox_inches = 'tight')        



        
'''
# ------------------------
#  GOEL DATA
# ------------------------
'''

# IMPORTING DATA ..




# CCDF
def Goel_ccdfs(Goel_topology,Goel_example_cascades,destinations,figname='',UKvsUS=False,image=False) :



    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(17.8/2+.7),cm_to_inch(6))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,3)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    LEGEND_SIZE = 4.5
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS

    CCDF_yticks = [0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = [100,1000,10000]
    col1_xticks = col0_xticks
    col2_xticks = [1,10,100]
    col3_xticks = [1,10,100]
    coldepth_xticks = [1,10,100]


    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxx_col0 = 100000
    maxx_col1 = 100000
    maxx_col2 = 1
    maxx_col3 = 100
    
    maxx_depth = 200
    maxy_CCDF = 1.3
    miny_CCDF = 0.0000007
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    coldepth_xlabel = 'Cascade Depth'

    col1_xlabel = 'Cascade Max Broadcast'
    col2_xlabel = 'Cascade Max Broadcast'
    col3_xlabel = 'Cascade Virality'

    CCDF_ylabel = 'CCDF'

    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'

    sbs_colors = {'blue':(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), 
                  'orange':(1.0, 0.4980392156862745, 0.054901960784313725), 
                  'green':(0.17254901960784313, 0.6274509803921569, 0.17254901960784313), 
                  'red':(0.8392156862745098, 0.15294117647058825, 0.1568627450980392), 
                  'purple':(0.5803921568627451, 0.403921568627451, 0.7411764705882353), 
                  'brown':(0.5490196078431373, 0.33725490196078434, 0.29411764705882354), 
                  'pink':(0.8901960784313725, 0.4666666666666667, 0.7607843137254902), 
                  'grey':(0.4980392156862745, 0.4980392156862745, 0.4980392156862745), 
                  'yellow':(0.7372549019607844, 0.7411764705882353, 0.13333333333333333), 
                  'cyan':(0.09019607843137255, 0.7450980392156863, 0.8117647058823529)
                 }
    
    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()

    
    
    

    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'

    if (UKvsUS == False  and image == False) :
        plotthis = sorted(Goel_topology['%s'%'news']['%s'%feature])
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='News')
        
        plotthis = sorted(Goel_topology['%s'%'videos']['%s'%feature])
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='Videos')

        plotthis = sorted(Goel_topology['%s'%'pictures']['%s'%feature])
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='Pictures')


        plotthis = sorted(Goel_topology['%s'%'petitions']['%s'%feature])
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='Petitions')           
        
    elif(UKvsUS == True) :
        plotthis = sorted(Goel_topology['%s'%'UK']['%s'%feature])
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='BBC+Guardian')

        plotthis = sorted(Goel_topology['%s'%'US']['%s'%feature])
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='NYTimes+HuffPost')

    elif (image == True) :
        plotthis = sorted(Goel_topology['%s'%'instagram.com']['%s'%feature])
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['pink'],linewidth=FT_linewidth,label='instagram.com')

        plotthis = sorted(Goel_topology['%s'%'twitpic.com']['%s'%feature])
        ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='twitpic.com')
        
        
    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([80,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    ax00.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        
    plt.text(0.033, 0.77, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)

    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False,fontsize=LEGEND_SIZE)

    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    # (0,1) Depth
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'max.depth'

    if (UKvsUS == False  and image == False) :
        plotthis = sorted(Goel_topology['%s'%'videos']['%s'%feature])
        ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='Videos')

        plotthis = sorted(Goel_topology['%s'%'pictures']['%s'%feature])
        ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='Pictures')

        plotthis = sorted(Goel_topology['%s'%'news']['%s'%feature])
        ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='News')

        plotthis = sorted(Goel_topology['%s'%'petitions']['%s'%feature])
        ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='Petitions')    

    elif (image == True) :
        plotthis = sorted(Goel_topology['%s'%'instagram.com']['%s'%feature])
        ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['pink'],linewidth=FT_linewidth,label='instagram.com')

        plotthis = sorted(Goel_topology['%s'%'twitpic.com']['%s'%feature])
        ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='twitpic.com')
        
        
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([0.9,maxx_depth])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col2_xticks)

    ax01.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
        
    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    


    # (0,2) BREADTH
    ax03 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'avg.dist'
    
    if (UKvsUS == False and image == False) :
        plotthis = sorted(Goel_topology['%s'%'videos']['%s'%feature])
        ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='Videos')

        plotthis = sorted(Goel_topology['%s'%'pictures']['%s'%feature])
        ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='Pictures')

        plotthis = sorted(Goel_topology['%s'%'news']['%s'%feature])
        ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='News')

        plotthis = sorted(Goel_topology['%s'%'petitions']['%s'%feature])
        ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='Petitions')    

    elif (UKvsUS == True) :
        plotthis = sorted(Goel_topology['%s'%'UK']['%s'%feature])
        ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='BBC+Guardian')

        plotthis = sorted(Goel_topology['%s'%'US']['%s'%feature])
        ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='NYTimes+HuffPost')
        
    elif (image == True) :
        plotthis = sorted(Goel_topology['%s'%'instagram.com']['%s'%feature])
        ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['pink'],linewidth=FT_linewidth,label='instagram.com')

        plotthis = sorted(Goel_topology['%s'%'twitpic.com']['%s'%feature])
        ax03.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='twitpic.com')
                
        
    ax03.set_xscale('log')
    ax03.set_yscale('log')
    ax03.set_xlim([1,maxx_col3])
    ax03.set_ylim([miny_CCDF,maxy_CCDF])

    ax03.set_yticks(CCDF_yticks)
    ax03.set_xticks(col3_xticks)

    ax03.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax03.set_yticklabels([])
        ax03.set_ylabel('')
        ax03.set_xlabel('')
         
    ax03.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    
    
    

    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'

    if (UKvsUS == False and image == False):
        if ('videos' in Goel_example_cascades.keys()) :
            plotthis = sorted(Goel_example_cascades['%s'%'videos']['%s'%feature])
            ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='Videos')
        if ('pictures' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'pictures']['%s'%feature])
            ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='Pictures')
        if ('news' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'news']['%s'%feature])
            ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='News')
        if ('petitions' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'petitions']['%s'%feature])
            ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='Petitions')    
    elif (UKvsUS == True) : 
        plotthis = sorted(Goel_example_cascades['%s'%'UK']['%s'%feature])
        ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='BBC+Guardian')

        plotthis = sorted(Goel_example_cascades['%s'%'US']['%s'%feature])
        ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='NYTimes+HuffPost')

    elif (image == True) :
        plotthis = sorted(Goel_example_cascades['%s'%'instagram.com']['%s'%feature])
        ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['pink'],linewidth=FT_linewidth,label='instagram.com')

        plotthis = sorted(Goel_example_cascades['%s'%'twitpic.com']['%s'%feature])        
        ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='twitpic.com')
                        
        
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([80,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)
    ax10.set_ylabel(CCDF_ylabel)

    ax10.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.033, 0.29+0.03, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)



    # INSERT SMALL IMAGE
    im_ax10 = inset_axes(ax10,
                        height="30%",
                        width="30%",
                        loc=3) 
    im_ax10.imshow(im_size)
    im_ax10.axis('off')

    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'max.depth'

    if (UKvsUS == False and image == False) :
        if ('videos' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'videos']['%s'%feature])
            ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='Videos')
        if ('pictures' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'pictures']['%s'%feature])
            ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='Pictures')
        if ('news' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'news']['%s'%feature])
            ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='News')
        if ('petitions' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'petitions']['%s'%feature])
            ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='Petitions')    
    
    elif (UKvsUS == True) : 
        plotthis = sorted(Goel_example_cascades['%s'%'UK']['%s'%feature])
        ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='BBC+Guardian')

        plotthis = sorted(Goel_example_cascades['%s'%'US']['%s'%feature])
        ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='NYTimes+HuffPost')
            
    elif (image == True) :
        plotthis = sorted(Goel_example_cascades['%s'%'instagram.com']['%s'%feature])
        ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['pink'],linewidth=FT_linewidth,label='instagram.com')

        plotthis = sorted(Goel_example_cascades['%s'%'twitpic.com']['%s'%feature])        
        ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='twitpic.com')
                      
    
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([0.9,maxx_depth])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col2_xticks)

    ax11.set_xlabel(coldepth_xlabel)

    ax11.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')
  
    ax11.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax11 = inset_axes(ax11,
                        height="30%", 
                        width="30%", 
                        loc=3) 
    im_ax11.imshow(im_depth)
    im_ax11.axis('off')

    
    
    # (1,2)
    ax13 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    feature = 'avg.dist'
    
    if (UKvsUS == False and image == False) :
        if ('videos' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'videos']['%s'%feature])
            ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='Videos')
        if ('pictures' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'pictures']['%s'%feature])
            ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='Pictures')
        if ('news' in Goel_example_cascades.keys()) :

            plotthis = sorted(Goel_example_cascades['%s'%'news']['%s'%feature])
            ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='News')

        if ('petitions' in Goel_example_cascades.keys()) :


            plotthis = sorted(Goel_example_cascades['%s'%'petitions']['%s'%feature])
            ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='Petitions')    
    
    elif (UKvsUS == True) : 
        plotthis = sorted(Goel_example_cascades['%s'%'UK']['%s'%feature])
        ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='BBC+Guardian')

        plotthis = sorted(Goel_example_cascades['%s'%'US']['%s'%feature])
        ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='NYTimes+HuffPost')

    elif (image == True) :
        plotthis = sorted(Goel_example_cascades['%s'%'instagram.com']['%s'%feature])
        ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['pink'],linewidth=FT_linewidth,label='instagram.com')

        plotthis = sorted(Goel_example_cascades['%s'%'twitpic.com']['%s'%feature])        
        ax13.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='twitpic.com')
                          
        
    
    ax13.set_xscale('log')
    ax13.set_yscale('log')
    ax13.set_xlim([1,maxx_col3])
    ax13.set_ylim([miny_CCDF,maxy_CCDF])
    ax13.set_yticks(CCDF_yticks)
    ax13.set_xticks(col3_xticks)

    ax13.set_xlabel(col3_xlabel)

    ax13.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax13.set_yticklabels([])
        ax13.set_ylabel('')
 
    ax13.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax13 = inset_axes(ax13,
                        height="30%", 
                        width="30%", 
                        loc=3) 
    im_ax13.imshow(im_virality)
    im_ax13.axis('off')    


    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()
    
    for destination in destinations :
        plt.savefig(destination+figname+'Goel.png',dpi=400)
    

def news_ccdfs(Goel_topology,Goel_example_cascades,destinations,figname='') :



    # NOW PLOT.
    plt.style.use('default')
    figsize2 = (cm_to_inch(8)/2*3,cm_to_inch(8))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (2,3)
    alpha50 = .50
    alpha90 = .25
    FT_linewidth = .75

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'left'
    valignment='center'

    # FONTS
    EVEN_SMALLER_SIZE = 5.5
    SMALL_SIZE = 6
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # TICKS
    CCDF_yticks = [0.0001,0.01,1]

    time_yticks = [10,1000,100000]

    col0_xticks = [100,1000,10000]
    col1_xticks = col0_xticks
    col2_xticks = [1,10]
    coldepth_xticks = [1,10,100]


    ax03_xticks = [0,20000,40000]
    ax13_xticks = [0,1000,2000]

    hide_CCDF_yticklabels = False

    # AXIS LIMS
    maxx_col0 = 100000
    maxx_col1 = 100000
    maxx_col2 = 100
    maxx_depth = 200

    maxy_CCDF = 1.3
    miny_CCDF = 0.000007
    maxy_time = 800000
    miny_time = 1.

    # LABELS
    col0_xlabel = 'Cascade Size'
    coldepth_xlabel = 'Cascade Depth'
    col2_xlabel = 'Cascade Virality'
    CCDF_ylabel = 'CCDF'


    # LEGEND
    legendhandlelength = 1

    horizontal_space_between_subplots = 0#1.
    vertical_space_between_subplots = .2


    # INSET NUMBERS
    inset_fontweight='bold'

    sbs_colors = {'blue':(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), 
                  'orange':(1.0, 0.4980392156862745, 0.054901960784313725), 
                  'green':(0.17254901960784313, 0.6274509803921569, 0.17254901960784313), 
                  'red':(0.8392156862745098, 0.15294117647058825, 0.1568627450980392), 
                  'purple':(0.5803921568627451, 0.403921568627451, 0.7411764705882353), 
                  'brown':(0.5490196078431373, 0.33725490196078434, 0.29411764705882354), 
                  'pink':(0.8901960784313725, 0.4666666666666667, 0.7607843137254902), 
                  'grey':(0.4980392156862745, 0.4980392156862745, 0.4980392156862745), 
                  'yellow':(0.7372549019607844, 0.7411764705882353, 0.13333333333333333), 
                  'cyan':(0.09019607843137255, 0.7450980392156863, 0.8117647058823529)
                 }
    
    def ccdf_y(arr) :
        return 1-np.arange(0,len(arr),1)/len(arr)
        
    def ccdf_x(arr) :
        return arr

    im_size,im_depth,im_max_breadth,im_virality,im_temporal = import_pictures()

    
    
    

    # ---
    # TOP: Original
    # ---

    # SIZE 
    # ----------
    # define
    ax00 = plt.subplot2grid(plot_dimension,(0,0),rowspan=1,colspan=1)
    feature = 'size'

    plotthis = sorted(Goel_topology['%s'%'bbc.co.uk']['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='bbc.co.uk')

    plotthis = sorted(Goel_topology['%s'%'nytimes.com']['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='nytimes.com')

    plotthis = sorted(Goel_topology['%s'%'guardian.co.uk']['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='guardian.co.uk')
    
    plotthis = sorted(Goel_topology['%s'%'huffingtonpost.com']['%s'%feature])
    ax00.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='huffingtonpost.com')    


    ax00.set_xscale('log')
    ax00.set_yscale('log')

    ax00.set_xlim([100,maxx_col0])
    ax00.set_ylim([miny_CCDF,maxy_CCDF])

    ax00.set_yticks(CCDF_yticks)
    ax00.set_xticks(col0_xticks)

    ax00.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax00.set_xlabel(col0_xlabel)
        ax00.set_xlabel('')
        
    plt.text(0.033, 0.77, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)

    ax00.legend(loc=3,handlelength=legendhandlelength,frameon=False)

    ax00.set_ylabel('A',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)

    # (0,1) Depth
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=1,colspan=1)
    feature = 'max.depth'

    
    plotthis = sorted(Goel_topology['%s'%'bbc.co.uk']['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='bbc.co.uk')

    plotthis = sorted(Goel_topology['%s'%'nytimes.com']['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='nytimes.com')

    plotthis = sorted(Goel_topology['%s'%'guardian.co.uk']['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='guardian.co.uk')
    
    plotthis = sorted(Goel_topology['%s'%'huffingtonpost.com']['%s'%feature])
    ax01.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='huffingtonpost.com')    

    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax01.set_xlim([1,maxx_depth])
    ax01.set_ylim([miny_CCDF,maxy_CCDF])

    ax01.set_yticks(CCDF_yticks)
    ax01.set_xticks(col2_xticks)

    ax01.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax01.set_yticklabels([])
        ax01.set_ylabel('')
        ax01.set_xlabel('')
          
    ax01.set_ylabel('B',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    


    
    # (0,2) BREADTH
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=1,colspan=1)
    feature = 'avg.dist'

    plotthis = sorted(Goel_topology['%s'%'bbc.co.uk']['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='bbc.co.uk')

    plotthis = sorted(Goel_topology['%s'%'nytimes.com']['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='nytimes.com')

    plotthis = sorted(Goel_topology['%s'%'guardian.co.uk']['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='guardian.co.uk')
    
    plotthis = sorted(Goel_topology['%s'%'huffingtonpost.com']['%s'%feature])
    ax02.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='huffingtonpost.com')    

    ax02.set_xscale('log')
    ax02.set_yscale('log')
    ax02.set_xlim([1,maxx_col2])
    ax02.set_ylim([miny_CCDF,maxy_CCDF])

    ax02.set_yticks(CCDF_yticks)
    ax02.set_xticks(col2_xticks)

    ax02.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax02.set_yticklabels([])
        ax02.set_ylabel('')
        ax02.set_xlabel('')
           
    ax02.set_ylabel('C',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)    


    # ---
    # BOTTOM: SAMPLED
    # ---

    # (1,0)
    ax10 = plt.subplot2grid(plot_dimension,(1,0),rowspan=1,colspan=1)
    feature = 'size'

    if ('bbc.co.uk' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'bbc.co.uk']['%s'%feature])
        ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='bbc.co.uk')
    if ('nytimes.com' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'nytimes.com']['%s'%feature])
        ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='nytimes.com')
    if ('guardian.co.uk' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'guardian.co.uk']['%s'%feature])
        ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='guardian.co.uk')
    if ('huffingtonpost.com' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'huffingtonpost.com']['%s'%feature])
        ax10.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='huffingtonpost.com')    
    
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.set_xlim([100,maxx_col0])
    ax10.set_ylim([miny_CCDF,maxy_CCDF])
    ax10.set_yticks(CCDF_yticks)
    ax10.set_xticks(col0_xticks)
    ax10.set_xlabel(col0_xlabel)
    ax10.set_ylabel(CCDF_ylabel)
  
    ax10.set_ylabel('D',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    plt.text(0.033, 0.29, 'CCDF',rotation=90, fontsize=SMALL_SIZE, transform=plt.gcf().transFigure)



    # INSERT SMALL IMAGE
    im_ax10 = inset_axes(ax10,
                        height="30%",
                        width="30%", 
                        loc=3) 
    im_ax10.imshow(im_size)
    im_ax10.axis('off')

    # (1,1)
    ax11 = plt.subplot2grid(plot_dimension,(1,1),rowspan=1,colspan=1)
    feature = 'max.depth'

    if ('bbc.co.uk' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'bbc.co.uk']['%s'%feature])
        ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='bbc.co.uk')
    if ('nytimes.com' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'nytimes.com']['%s'%feature])
        ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='nytimes.com')
    if ('guardian.co.uk' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'guardian.co.uk']['%s'%feature])
        ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='guardian.co.uk')
    if ('huffingtonpost.com' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'huffingtonpost.com']['%s'%feature])
        ax11.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='huffingtonpost.com')    
    
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.set_xlim([1,maxx_depth])
    ax11.set_ylim([miny_CCDF,maxy_CCDF])
    ax11.set_yticks(CCDF_yticks)
    ax11.set_xticks(col2_xticks)

    ax11.set_xlabel(coldepth_xlabel)

    ax11.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax11.set_yticklabels([])
        ax11.set_ylabel('')

    ax11.set_ylabel('E',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax11 = inset_axes(ax11,
                        height="30%",
                        width="30%",
                        loc=3)
    im_ax11.imshow(im_depth)
    im_ax11.axis('off')

    # (1,2)
    ax12 = plt.subplot2grid(plot_dimension,(1,2),rowspan=1,colspan=1)
    feature = 'avg.dist'
    if ('bbc.co.uk' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'bbc.co.uk']['%s'%feature])
        ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['purple'],linewidth=FT_linewidth,label='bbc.co.uk')
    if ('nytimes.com' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'nytimes.com']['%s'%feature])
        ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['blue'],linewidth=FT_linewidth,label='nytimes.com')
    if ('guardian.co.uk' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'guardian.co.uk']['%s'%feature])
        ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['orange'],linewidth=FT_linewidth,label='guardian.co.uk')
    if ('huffingtonpost.com' in Goel_example_cascades.keys()) :

        plotthis = sorted(Goel_example_cascades['%s'%'huffingtonpost.com']['%s'%feature])
        ax12.plot(ccdf_x(plotthis),ccdf_y(plotthis),color=sbs_colors['green'],linewidth=FT_linewidth,label='huffingtonpost.com')    
    
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_xlim([1,maxx_col2])
    ax12.set_ylim([miny_CCDF,maxy_CCDF])
    ax12.set_yticks(CCDF_yticks)
    ax12.set_xticks(col2_xticks)

    ax12.set_xlabel(col2_xlabel)

    ax12.set_ylabel(CCDF_ylabel)

    if (hide_CCDF_yticklabels==True) :
        ax12.set_yticklabels([])
        ax12.set_ylabel('') 
    ax12.set_ylabel('F',verticalalignment='top',y=1,rotation=0,fontsize=MEDIUM_SIZE,fontweight=inset_fontweight)  

    # INSERT SMALL IMAGE
    im_ax12 = inset_axes(ax12,
                        height="30%",
                        width="30%", 
                        loc=3)
    im_ax12.imshow(im_virality)
    im_ax12.axis('off')


    fig.align_ylabels([ax00,ax10])
    plt.tight_layout()
    
    for destination in destinations :
        plt.savefig(destination+figname+'Goelnews.png',dpi=400)     
