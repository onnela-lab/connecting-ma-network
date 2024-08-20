# -*- coding: utf-8 -*-
"""
Check function SIR new approach
Last check 09/03/22, 12h16

@author: thl902
"""



import numpy as np
import math



def SIR_Spreadfunction(theta, timesteps, risk_matrix, initial_status):
    """"
    Args: 
          theta: parameter
          timesteps: number time steps
          risk_matrix: the transmission matrix
          initial_status: the initial status, number S,I,R
    Outputs: Status matrix with 3 columns are number of S,I,R each timestep
             
    """

    tmp_matrix = np.zeros((timesteps,3)) #number S,I, and R each step
    initial_status1 = list(np.array(initial_status)) # assign another id
    risk_matrix1= np.array(risk_matrix)
    N = sum(initial_status)
    tmp_matrix[0,:] = initial_status1
    #update second row of the matrix
    x = tmp_matrix[0,:] 
    total_infected = x[1]+x[2]
    
    transmission_rate = theta[0]*x[1]/total_infected*sum(risk_matrix1[int(total_infected-1),0:N])
         
    y1 = np.random.poisson(transmission_rate, 1) #new infected, poisson(transmission_rate)
    y2 = np.random.binomial(x[1], theta[1], 1)
        #Update S,I,R status
        ##   S
    if y1 <= x[0]:
            tmp0 = x[0] - y1
    else:
            y1 = x[0]
            tmp0 = 0

    tmp1 = x[1] + y1 - y2      
    tmp2 = N - tmp0 - tmp1
        
    tmp_matrix[1,0] = tmp0
    tmp_matrix[1,1] = tmp1
    tmp_matrix[1,2] = tmp2
        
    
    for i in range(2,timesteps):
     
        x1 = tmp_matrix[(i-2),:]
        x2 = tmp_matrix[(i-1),:]
        new_recovered = x2[2] - x1[2]
        total_infect1 = x1[1]+x1[2]
        total_infect2 = x2[1]+x2[2]
                
            
        total_rate1 = (x1[1]-new_recovered)/total_infect1*sum(risk_matrix1[int(total_infect2-1),0:int(total_infect1)])    
        total_rate2 = sum(risk_matrix1[int(total_infect2-1),int(total_infect1):int(total_infect2)])
        
     
        transmission_rate = theta[0]*( total_rate1  + total_rate2 )
        
                
        if min((N - total_infect2), transmission_rate)>0:
            y1 = np.random.poisson(transmission_rate, 1) #new infected, poisson(transmission_rate)
        else: 
            y1 = 0
           
        y2 = np.random.binomial(x2[1], theta[1], 1)
        #Update S,I,R status
        ##   S
        if y1 <= x2[0]:
                tmp0 = x2[0] - y1
        else:
                y1 = x2[0]
                tmp0 = 0
        
        ##   I
        tmp1 = x2[1] + y1 - y2 #as y2 always less than x[1]
                          
        tmp2 = N - tmp0 - tmp1
        
        tmp_matrix[i,0] = tmp0
        tmp_matrix[i,1] = tmp1
        tmp_matrix[i,2] = tmp2
        
    return(tmp_matrix)


#####Example usage
# import os
# os.chdir("C:\\Users\\lemin\\Desktop\\HIV\\1.Manuscript_ClassicalvsNetwork\\SIR_Approximate\\")


# import networkx as nx
# import matplotlib.pyplot as plt
# from  Avetransmissionmatrix_function import Avetransmissionmatrix_function


# N = 100

# theta = [1.1,0.1]
# G = nx.Graph()
# p=.5
# G= nx.erdos_renyi_graph(N,p)


# ##Add edges 0,1,2,..,N-1 to make one component
# for i in range(0,N-1):
#       G.add_edge(i,i+1)
      
# nx.draw(G,with_labels=True)
# plt.show()
      
     
# Ini1 = 0 
# Nrep = 100    
# risk_matrix = Avetransmissionmatrix_function(G, Ini1, Nrep)
# timesteps = 100
# initial_status = [N-1,1,0]
# SIR_Spreadfunction(theta, timesteps, risk_matrix, initial_status)


