# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 15:30:55 2022

@author: thl902
"""

import numpy as np
from numpy.random import choice

def SIR_Spreadfunction_Network(G, theta, timesteps, initial_infectednodes):
    """"
    Args: G: Given network
          theta: parameter
          timesteps: number time steps
          initial_infectednodes: initial infected node
    Outputs: Status matrix, the first 3 columns are number of S,I,R each timestep, 
             
    
    """
    status_matrix = np.zeros((timesteps, 3))
    infected_status_matrix = np.zeros((len(G.nodes()),timesteps))
    recovered_status_matrix = np.zeros((len(G.nodes()),timesteps))
    i_nodes1 = np.array(initial_infectednodes) #make another id for i_nodes
    i_nodes = list(i_nodes1) # assig infected nodes
    r_nodes = list()
    s_nodes = list ( set (G.nodes ()) - set (i_nodes ) - set (r_nodes ) )
   
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
            
        at_risk_nodes = set(np.unique(at_risk_nodes)) - set(i_nodes) - set(r_nodes) # update at risk nodes
        
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
                p_rate = min(transmission_rate/len(at_risk_nodes),1)
                new_infected = np.random.binomial(len(at_risk_nodes),p_rate, 1) # be binomial for a closer approximate for arbitrary network
            
                new_infected_nodes = choice(at_risk_nodes, new_infected, replace=False, p= prob_sequence) 
                
                    
        #########Step 2, Second update recovered node
        if len(i_nodes)==0:
            new_recovered_nodes = []
        else:
            new_recovered = np.random.binomial(len(i_nodes), theta[1], 1) # be binomial for a closer approximate for arbitrary network
    
            new_recovered_nodes = choice(i_nodes, new_recovered, replace=False) 
        #update r nodes and i nodes
        
        
        #assign the number of recovered nodes to infected nodes
        #recovered_rate = theta[1]*len(i_nodes)
        
        
        

        ##update infected nodes and recover nodes
        for node in new_infected_nodes:
            s_nodes.remove( node)
            i_nodes.append( node) #
        
        for node in new_recovered_nodes:
            r_nodes.append( node)
            i_nodes.remove(node) #
        
       
        ##########Update recovered and infected nodes in the matrix
        for node in i_nodes:
            infected_status_matrix[node,step] = 1
            
        for node in r_nodes:
            recovered_status_matrix[node,step] = 1
         
    
    status_matrix[:,1] = infected_status_matrix.sum(axis=0, dtype='float')
    status_matrix[:,2] = recovered_status_matrix.sum(axis=0, dtype='float')
    status_matrix[:,0]  = np.ones((timesteps))*len(G.nodes()) - status_matrix[:,1]-status_matrix[:,2]
        
    return (status_matrix)
