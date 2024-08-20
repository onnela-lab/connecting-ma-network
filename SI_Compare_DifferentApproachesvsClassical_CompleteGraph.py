# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 14:55:33 2022
Make figure in the manuscript compare the matching of the conventional approach, our approach, and the classical model 
@author: lemin
"""


import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import os
import random
import time
os.chdir("C:\\Users\\pwr844\\OneDrive - University of Tennessee\\Attachments\\Desktop\\2.Connecting\\SI_fullyconnected\\")

os.getcwd()
plot_dir = "C:\\Users\\pwr844\\OneDrive - University of Tennessee\\Attachments\\Desktop\\2.Connecting\\SI_fullyconnected\\"

from SI_Spreadfunction_Network  import SI_Spreadfunction_Network, SI_Spreadfunction_Networkconvention, SI_Spreadfunction_Massaction

########################################################################################

start_time = time.time()

N = 100 # Graph with N nodes, increase to 1000
nrepeat = 100 #number stochastic realizations, increase to 300
#timesteps = 25
timesteps = 100


G = nx.complete_graph(N)


# nx.draw(G,with_labels=True)
# plt.show()


Ini1 = 0
initial_infectednodes = [Ini1]



initial_status = [N-len(initial_infectednodes),len(initial_infectednodes)]




beta = .12
#beta = .7
theta = [beta] #

thetain = int(theta[0]*100)



###########

##Plotting the average realization to  check the match
#1. Average realization of the three models

realization_infected_matrix_massaction = np.zeros((nrepeat,timesteps)) #keep track total infected at each time step for different realizations
realization_infected_matrix_network_convention = np.zeros((nrepeat,timesteps)) #keep track total infected at each time step for different realizations
realization_infected_matrix_proposednetwork = np.zeros((nrepeat,timesteps)) #keep track total infected at each time step for different realizations
for i in range(nrepeat):
    tmp1 = SI_Spreadfunction_Network(G, theta, timesteps, initial_infectednodes)
    realization_infected_matrix_proposednetwork[i,:] = tmp1[:,1]
    tmp2 = SI_Spreadfunction_Networkconvention(G, theta, timesteps, initial_infectednodes)
    realization_infected_matrix_network_convention[i,:] = tmp2[:,1]
    tmp3 = SI_Spreadfunction_Massaction(theta, timesteps, initial_status)
    realization_infected_matrix_massaction[i,:] = tmp3[:,1]
    print('This is the replicate',i/nrepeat)
########



###################
#Mean
mean_infected_matrix_network_massaction = np.mean(realization_infected_matrix_massaction,axis=0)/N
mean_infected_matrix_network_convention = np.mean(realization_infected_matrix_network_convention,axis=0)/N
mean_infected_matrix_proposednetwork = np.mean(realization_infected_matrix_proposednetwork,axis=0)/N


# index = np.linspace(0, 99, num = 30)
# index = [ int(x) for x in index ]

# mean_infected_matrix_network_massaction = mean_infected_matrix_network_massaction[index]
# mean_infected_matrix_network_convention = mean_infected_matrix_network_convention[index]
# mean_infected_matrix_proposednetwork = mean_infected_matrix_proposednetwork[index]

h1,  = plt.plot (mean_infected_matrix_network_massaction,"-", linewidth=1)
h2,  = plt.plot (mean_infected_matrix_network_convention,"-+", linewidth=1)
h3,  = plt.plot (mean_infected_matrix_proposednetwork,"-x", linewidth=1)


plt.legend ([h1,h2,h3], ["Mass-action" , "Convention",
                                "Proposed"],
            prop={'size': 7})

plt.xlabel ("Time steps")
plt.ylabel ("Infected proportion ")

xname = list(range(0,timesteps,2))
plt.xticks(xname)

#plt.title('Mean total infected for a complete graph of N=%i nodes, theta=%.2f'% (N,beta))

#plt.xticks(np.arange(min(index), max(index)+1, 3))
# plt.xlim ([0 , (timesteps - 1)])


#plot_name = "SI_CompleteGraph_DifferentApproaches_" + str(N) + "nodes.pdf"

#plt.savefig(plot_dir + plot_name)
#plt.show()
# #################


total_time = time.time() -start_time

total_time
# Displaying the array
print('Array:\n', realization_infected_matrix_proposednetwork)
 
# Saving the array in a text file
filename1 = "RealizationsProposed_"+ str(thetain) + "parameter_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

np.savetxt(filename1, realization_infected_matrix_proposednetwork)
# Saving the array in a text file
filename2 = "RealizationsConvention_"+ str(thetain) + "parameter_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
np.savetxt(filename2, realization_infected_matrix_network_convention)
 
# Saving the array in a text file
filename3 = "RealizationsMassaction_"+ str(thetain) + "parameter_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

np.savetxt(filename3, realization_infected_matrix_massaction)
  
time.time() - start_time
