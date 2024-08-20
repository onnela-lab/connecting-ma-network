
import numpy as np

from numpy.random import choice


def SITAD_Spreadfunction_Network(G, theta, timesteps, initial_infectednodes):
    """"
    Args: G:network
          theta = (beta1, beta2, gamma1, delta1, gamma2,  delta2)
          timesteps: number time steps
          initial_status: number initial infected
          topo_seq: transmitted due to topology at each step
    Outputs: Status matrix with 5 columns are number of S,I,T,A,D each timestep, 
             
    
    """
    status_matrix = np.zeros((timesteps, 5))
    S_status_matrix = np.zeros((len(G.nodes()),timesteps))
    I_status_matrix = np.zeros((len(G.nodes()),timesteps))
    T_status_matrix = np.zeros((len(G.nodes()),timesteps))
    A_status_matrix = np.zeros((len(G.nodes()),timesteps))
    D_status_matrix = np.zeros((len(G.nodes()),timesteps))
    i_nodes1 = np.array(initial_infectednodes) #make another id for i_nodes
    i_nodes = list(i_nodes1) # assig infected nodes
    t_nodes = list()
    a_nodes = list()
    d_nodes = list()
    s_nodes = list ( set (G.nodes ()) - set(i_nodes) - set(t_nodes) - set(a_nodes) - set(d_nodes))
       
    for node_order in range(len(i_nodes)):
        I_status_matrix[i_nodes[node_order],0] = 1 # assign infected nodes to the first time step of the infected matrix
       
    # Loop over all time steps .
    for step in range(1, timesteps):
        #new infections at each step and update s_nodes, i_nodes
        #Step1: Extract all at risk nodes and calculate their risk rates
        #1a. Extract all at risk nodes, neigbor of i nodes and a nodes
      
        at_risk_nodes = []
        for i in range(len(i_nodes)):
            node = i_nodes[i]
            at_risk_nodes = at_risk_nodes + list(G.neighbors(node))
        
        for i1 in range(len(a_nodes)):
            node = a_nodes[i1]
            at_risk_nodes = at_risk_nodes + list(G.neighbors(node))
            
        at_risk_nodes = set(np.unique(at_risk_nodes)) - set(i_nodes) - set(t_nodes) - set(a_nodes) - set(d_nodes) # update at risk nodes
        
        at_risk_nodes = list(at_risk_nodes) # make set of updated at_risk_nodes becomes an array
        
        #1b. Total transmission rate over I nodes and A nodes
        infectedrate_sequence_I = []
        for i in i_nodes:
            tmp = set(G.neighbors(i))
            tmp1 = tmp & set(s_nodes)
            infectedrate_sequence_I.append(len(tmp1)/(len(tmp)+1))
        
        infectedrate_sequence_A = []
        for i in a_nodes:
            tmp2 = set(G.neighbors(i))
            tmp3 = tmp2 & set(s_nodes)
            infectedrate_sequence_A.append(len(tmp3)/(len(tmp2)+1))
        
        infected_ratesI = np.sum(infectedrate_sequence_I) 
        infected_ratesA = np.sum(infectedrate_sequence_A)    
        #1c. Check out number I and A neighbors of each at risk nodes
        ## Assign risk weight each at risk node
        risk_rate_sequence = []
        for i in at_risk_nodes:
            tmp = set(G.neighbors(i))
            tmp1 = tmp & set(i_nodes)
            tmp2 = tmp & set(a_nodes)
            risk_rate_sequence.append(theta[0]*len(tmp1)+theta[1]*len(tmp2))
        prob_sequence = risk_rate_sequence/np.sum(risk_rate_sequence) 
        
        #1d. I treated rate,  AIDS progress rate, 
        #AIDS treated rate, Death rate
        transmission_rate = theta[0]*infected_ratesI + theta[1]*infected_ratesA
           
        
        #Step 2a: Generate new infected I
        
        new_I_nodes = []
        if len(at_risk_nodes)>0:
            
            if transmission_rate/len(at_risk_nodes)>= 1:
                 new_I_nodes = np.array(at_risk_nodes)
            else:
                new_I = np.random.binomial(len(at_risk_nodes),transmission_rate/len(at_risk_nodes), 1) # be binomial for a closer approximate for arbitrary network
            
                new_I_nodes = choice(at_risk_nodes, new_I, replace=False, p= prob_sequence) 
        
           
        #Step 2b:   Generate new I treated     
        new_treatednodes_I = []
        if len(i_nodes)>0:
            new_treatedI = np.random.binomial(len(i_nodes), theta[2], 1) # be binomial for a closer approximate for arbitrary network
            new_treatednodes_I = choice(i_nodes, new_treatedI, replace=False) 
              
        #Step 2c: Generate new AIDS ordered based on the infection duration
        
        new_A_nodes = []
        i_nodes_remain = np.setdiff1d(i_nodes, new_treatednodes_I)
        new_A = np.random.binomial(len(i_nodes),theta[3], 1) 
        
            
        if len(i_nodes_remain) > new_A:
            new_A_nodes = choice(i_nodes_remain, new_A, replace=False)
        else:
            new_A_nodes = i_nodes_remain
        
        #Step 2d: Generate new A treated
        if len(a_nodes)==0:
            new_treatednodes_A = []
        else:
            new_treatedA = np.random.binomial(len(a_nodes), theta[4], 1) # be binomial for a closer approximate for arbitrary network
            new_treatednodes_A = choice(a_nodes, new_treatedA, replace=False) 
               
        #Step 2e: Generate new death D by AIDS, ordered based on the AIDS duration
        new_D_nodes = []
        a_nodes_remain = np.setdiff1d(a_nodes, new_treatednodes_A)
        new_D = np.random.binomial(len(a_nodes),theta[5], 1)
        
       
        if len(a_nodes_remain) > new_D:
            new_D_nodes = choice(a_nodes_remain, new_D, replace=False)
        else:
            new_D_nodes = a_nodes_remain
#################################################################
        
        

        #########Step 3, Update  
        ##update infected nodes I
        for node in new_I_nodes:
            s_nodes.remove( node)
            i_nodes.append( node) #
        
        ##update new treated nodes I
        for node in new_treatednodes_I:
            t_nodes.append( node)
            i_nodes.remove( node) #
        
        ##update new AIDS nodes A
        for node in new_A_nodes:
            a_nodes.append( node)
            i_nodes.remove( node) #
            
        ##update new treated nodes A
        for node in new_treatednodes_A:
            t_nodes.append( node)
            a_nodes.remove( node) #    
        ##update death nodes by AIDS
        for node in new_D_nodes:
            d_nodes.append( node)
            a_nodes.remove( node) #  
            
        
           
        ##########Update tracking matrix of each status
        for node in s_nodes:
            S_status_matrix[node,step] = 1
        
        
        for node in i_nodes:
            I_status_matrix[node,step] = 1
        
            
        for node in t_nodes:
            T_status_matrix[node,step] = 1
            
        for node in a_nodes:
            A_status_matrix[node,step] = 1    
            
        for node in d_nodes:
            D_status_matrix[node,step] = 1    
            
        
         
           ###########finalized number infected, recover over time
    status_matrix[:,0] = S_status_matrix.sum(axis=0, dtype='float')
    status_matrix[:,1] = I_status_matrix.sum(axis=0, dtype='float')
    status_matrix[:,2] = T_status_matrix.sum(axis=0, dtype='float')
    status_matrix[:,3] = A_status_matrix.sum(axis=0, dtype='float')
    status_matrix[:,4] = D_status_matrix.sum(axis=0, dtype='float')
    
    

    return(status_matrix)

