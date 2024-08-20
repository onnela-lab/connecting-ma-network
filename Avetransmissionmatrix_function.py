# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:21:04 2022
Purpose: Extract the average transmission matrix and at risk nodes for a given sequence of infected ordered
@author: thl902
"""

import numpy as np

from numpy.random import choice


########The first function to extract the transmission matrix for each given node order
def Transmissionmatrix_function(G, initial_infectednodes):
    """"
    Args: G: Given network
          initial_infectednodes: initial infected node
    Outputs: Transmission matrix, column N+1: Number at risk nodes, 
             Column N+2: Sequence order of infected nodes 
             
    """
    
    N = len(list(G.nodes () ))
    risk_matrix = np.zeros((N, (N+1)))
    i_nodes1 = np.array(initial_infectednodes) #make another id for i_nodes
    i_nodes = list(i_nodes1) # assig infected nodes
    s_nodes = list (set (G.nodes ()) - set (i_nodes))
    order_infected = np.zeros(N)
    order_infected[0] = initial_infectednodes[0]
    
# Loop over all time steps 
# =============================================================================
    for step in range(1, N):
        
        
       at_risk_nodes = []     
       #Risk matrix
       for i in i_nodes:
           tmp = set(G.neighbors(i))
           tmp1 = tmp & set(s_nodes)
           at_risk_nodes = at_risk_nodes + list(G.neighbors(i))
           risk_matrix[(step-1),int(i)] = len(tmp1)/(len(tmp)+1) #susceptible neighbors/(neighbors+1)
       
       at_risk_nodes = set(np.unique(at_risk_nodes)) - set(i_nodes)  # update at risk nodes
       at_risk_nodes = list(at_risk_nodes) #make set of updated at_risk_nodes becomes an array
        
       risk_matrix[(step-1),N] = len(at_risk_nodes)  #number at risk nodes
       
       risk_rate_sequence = []
       for i in at_risk_nodes:
           tmp = set(G.neighbors(i)) & set(i_nodes)
           risk_rate_sequence.append(len(tmp))
            
        
       prob_sequence = risk_rate_sequence/np.sum(risk_rate_sequence) 
       new_infected_nodes = []
       if len(at_risk_nodes)>0:
           new_infected_nodes = choice(at_risk_nodes, 1, replace=False, p= prob_sequence) 
            #fill order infected node   
           order_infected[step] = new_infected_nodes 
      
       
      # new_infected_nodes1 = [new_infected_nodes]
       
       ##update infected nodes and recover nodes
       for node in new_infected_nodes:
           s_nodes.remove(node)
           i_nodes.append(node) #      
    order =   list(map(int, list(order_infected)))   
    risk_matrix[:,0:N] =  risk_matrix[:,order] 
              
        ###########
    return(risk_matrix)
   

########The seccond function to extract the average transmission matrix

def Avetransmissionmatrix_function(G, Ini1, Nrep):
    """"
    Args: G: Given network
          Ini1: initial infected node, just a number as the label
          Nrep: Number replications average
          
    Outputs: Average transmission matrix
             
    """
    initial_infectednodes = [Ini1]
    N = len(list(G.nodes()))
    risk_matrix_sum = np.zeros((N,(N+1)))
    

    for i in range(Nrep):
        risk_matrix = Transmissionmatrix_function(G, initial_infectednodes)      
        risk_matrix_sum = risk_matrix  + risk_matrix_sum
    
    risk_matrix_ave = risk_matrix_sum/Nrep
        
    return(risk_matrix_ave)      
   


# # Example usages
#####################
# import networkx as nx
# import matplotlib.pyplot as plt
# import numpy as np

# N = 5
# p=.5
# G = nx.Graph()
# G= nx.erdos_renyi_graph(N,p)

# ##Add edges 0,1,2,..,N-1 to make one component
# for i in range(0,N-1):
#       G.add_edge(i,i+1)

# G = nx.balanced_tree(1,(N-1)) #LINE GRAPH OF N NODES    
# #G= nx.complete_graph(N)
# nx.draw(G,with_labels=True)
# plt.show()
# Ini1 = 0
# initial_infectednodes = [Ini1]

# Transmissionmatrix_Merged(G, initial_infectednodes)
# Nrep = 100    
# Avetransmissionmatrix_function(G, Ini1, Nrep)