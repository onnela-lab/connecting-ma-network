# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 15:30:55 2022
The proposed spreading function, conventional spreading function, and classical massaction spreading function
@author: thl902
"""

import numpy as np
import random
from numpy.random import choice

#######Our proposed spreading rule that adopts the Gilesspie

def SI_Spreadfunction_Network(G, theta, timesteps, initial_infectednodes):
    """"
    Args: G: Given network
          theta: parameter
          timesteps: number time steps
          initial_infectednodes: initial infected node
    Outputs: Status matrix,  columns are number of S,I each timestep, 
             
    """
    status_matrix = np.zeros((timesteps, 2))
    infected_status_matrix = np.zeros((len(G.nodes()),timesteps))
    i_nodes1 = np.array(initial_infectednodes) #make another id for i_nodes
    i_nodes = list(i_nodes1) # assig infected nodes
    s_nodes = list ( set (G.nodes ()) - set (i_nodes ))
   
    for node_order in range(len(i_nodes)):
        infected_status_matrix[i_nodes[node_order],0] = 1 # assign infected nodes to the first time step of the infected matrix
       
# Loop over all time steps .
    for step in range(1, timesteps):
       
        #new infections at each step and update s_nodes, i_nodes
        #Step1: Extract all at risk nodes and calculate their risk rates
#1a. Extract all at risk nodes
        at_risk_nodes = []
        for i in range(len(i_nodes)):
            node = i_nodes[i]
            at_risk_nodes = at_risk_nodes + list(G.neighbors(node))
            
        at_risk_nodes = set(np.unique(at_risk_nodes)) - set(i_nodes)  # update at risk nodes
        
        at_risk_nodes = list(at_risk_nodes) # make set of updated at_risk_nodes becomes an array
        
        #2b. Total transmission rate over infected nodes
        infectedrate_sequence = []
        for i in i_nodes:
            tmp = set(G.neighbors(i))
            tmp1 = tmp & set(s_nodes)
            infectedrate_sequence.append(len(tmp1)/(len(tmp)+1)) #susceptible neighbors/(neighbors+1)
            
        total_infected_rates = np.sum(infectedrate_sequence)    
        #1c. Risk rates of at risk nodes, number of infected neigbors
        risk_rate_sequence = []
        for i in at_risk_nodes:
            tmp = set(G.neighbors(i))
            tmp1 = tmp & set(i_nodes)
            risk_rate_sequence.append(len(tmp1))
            
        
        prob_sequence = risk_rate_sequence/np.sum(risk_rate_sequence) 
        #Step 2: Assigned new infected nodes to the at_risk_nodes
        #2a. Generate number of new infected nodes
        #The reproduction number Rt: p_si*\sum_{i in infectednodes}S_i/(N_i + 1)
        transmission_rate = theta[0]*total_infected_rates
        
        #assign the number of infected to at_risk_nodes
        
        new_infected_nodes = []
        if len(at_risk_nodes)>0:
            
            if transmission_rate/len(at_risk_nodes)>= 1:
                 new_infected_nodes = np.array(at_risk_nodes)
            else:
                new_infected = np.random.binomial(len(at_risk_nodes),transmission_rate/len(at_risk_nodes), 1) # be binomial for a closer approximate for arbitrary network
            
                new_infected_nodes = choice(at_risk_nodes, new_infected, replace=False, p= prob_sequence) 
                
        ##update infected nodes and recover nodes
        for node in new_infected_nodes:
            s_nodes.remove( node)
            i_nodes.append( node) #
        
      
        ##########Update recovered and infected nodes in the matrix
        for node in i_nodes:
            infected_status_matrix[node,step] = 1
        
         
    
    status_matrix[:,1] = infected_status_matrix.sum(axis=0, dtype='float')
   
    status_matrix[:,0]  = np.ones((timesteps))*len(G.nodes()) - status_matrix[:,1]
        
    return (status_matrix)
##############################################

##The convention spreading rule on network, Newman, Pastor

def SI_Spreadfunction_Networkconvention(G, theta, timesteps, initial_infectednodes):
    status_matrix = np.zeros((timesteps, 2))
    infected_status_matrix = np.zeros((len(G.nodes()),timesteps))
    N = len(G.nodes())
    i_nodes1 = np.array(initial_infectednodes) #make another id for i_nodes
    i_nodes = list(i_nodes1) # assig infected nodes
    s_nodes = list ( set (G.nodes ()) - set (i_nodes ))
    for node_order in range(len(i_nodes)):
        infected_status_matrix [i_nodes[node_order],0] = 1 # assign 
# Loop over all time steps .
    for step in range(1, timesteps):
        #new infections at each step and update s_nodes, i_nodes
        new_infections = set()
        for node in i_nodes :
            if G.degree(node) > 0:
               node_neighbors = G.neighbors( node)
               for node1 in  node_neighbors:
                   h1  = random.random ()
                   if (node1 in s_nodes) and (h1 < theta[0]/N):
                      new_infections.add(node1)
        for node in new_infections :
            i_nodes.append( node)
            s_nodes.remove( node) #return newinfection, also update s_nodes and i_nodes
        for node in i_nodes:
            infected_status_matrix [node,step]=1
    status_matrix[:,1] = infected_status_matrix.sum(axis=0, dtype='float')
    status_matrix[:,0]  = np.ones((timesteps))*len(G.nodes()) - status_matrix[:,1]
          
    return (status_matrix)

##########

##Define the SI spreading with a mass action model 

def SI_Spreadfunction_Massaction(theta,timesteps, initial_status):
    """"
    Args: 
          theta: parameter
          timesteps: number time steps
          initial_infectednodes: initial status
          N: network size
    Outputs: Status matrix with 2 columns are number of S,I each timestep
             
    """
    tmp_matrix = np.zeros((timesteps,2)) #number S,I, and R each step
    initial_status1 = list(np.array(initial_status)) # assign another id
    #id(initial_infected) is id(initial_infected1)
    tmp_matrix[0,:] = initial_status1
    N = np.sum(initial_status)
    
    for i in range(1,timesteps):
        x = tmp_matrix[(i-1),:]
        transmission_rate = theta[0]*x[0]*x[1]/N
        #replace by binom
        if x[0]>0:
           y1 = np.random.binomial(x[0], transmission_rate/x[0], 1) #new infected, binom
        else:
            y1 = 0
            
        if y1 <= x[0]:
                tmp0 = x[0] - y1
        else:
                y1 = x[0]
                tmp0 = 0
    
        tmp1 = N - tmp0 
        tmp_matrix[i,0] = tmp0
        tmp_matrix[i,1] = tmp1
     
    return(tmp_matrix)
# ##Usage for mass action
# theta = [.4]
# timesteps = 100
# initial_status = [99,1]

# theta[0]
# SI_Spreadfunction_Massaction(theta,timesteps, initial_status)