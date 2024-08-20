# -*- coding: utf-8 -*-
"""
Created on Wed 03, 08, 2023 11:56:08 2023
Plots for the Approximate for different spreading processes. 
@author: thl902
"""



import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import os
import random
import time
os.chdir("C:\\Users\\pwr844\\OneDrive - University of Tennessee\\Attachments\\Desktop\\2.Connecting\\SimulationApproximate\\")

os.getcwd()
plot_dir = "C:\\Users\\pwr844\\OneDrive - University of Tennessee\\Attachments\\Desktop\\2.Connecting\\SimulationApproximate\\"




theta = [.2] #

thetain = int(theta[0]*100)


N = 1000 # Graph with N nodes

nrepeat = 1000 #number stochastic realizations
timesteps = 100

#######SI 
# Load txt file
filename1 = "SI_BA_RealizationsRiskMatrix_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

# Saving the array in a text file
filename2 = "SI_BA_RealizationsRiskMatrixAve_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
 
# Saving the array in a text file
filename3 = "SI_BA_RealizationsNetwork_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

# Displaying the contents of the text file
realization_infected_riskmatrix1_SI =  np.loadtxt(filename1)
realization_infected_riskmatrix_ave1_SI = np.loadtxt(filename2)
realization_infected_matrix_network1_SI = np.loadtxt(filename3)


mean_infected_riskmatrix1_SI = np.mean(realization_infected_riskmatrix1_SI,axis=0)/N 
mean_infected_riskmatrix_ave1_SI = np.mean(realization_infected_riskmatrix_ave1_SI,axis=0)/N 
mean_infected_matrix_network1_SI = np.mean(realization_infected_matrix_network1_SI,axis=0)/N


filename4 = "SI_ER_RealizationsRiskMatrix_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

# Saving the array in a text file
filename5 = "SI_ER_RealizationsRiskMatrixAve_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
 
# Saving the array in a text file
filename6 = "SI_ER_RealizationsNetwork_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

# Displaying the contents of the text file
realization_infected_riskmatrix2_SI =  np.loadtxt(filename4)
realization_infected_riskmatrix_ave2_SI = np.loadtxt(filename5)
realization_infected_matrix_network2_SI = np.loadtxt(filename6)


mean_infected_riskmatrix2_SI = np.mean(realization_infected_riskmatrix2_SI,axis=0)/N 
mean_infected_riskmatrix_ave2_SI = np.mean(realization_infected_riskmatrix_ave2_SI,axis=0)/N 
mean_infected_matrix_network2_SI = np.mean(realization_infected_matrix_network2_SI,axis=0)/N

##############SIR

# Load txt file
filename1a = "SIR_BA_Realizationsinfected_Network_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

filename2a = "SIR_BA_Realizationsinfected_Riskmatrix_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

filename3a = "SIR_BA_Realizationsinfected_Riskmatrixave_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

# Displaying the contents of the text file
realization_infected_matrix_network1_SIR = np.loadtxt(filename1a)
realization_infected_riskmatrix1_SIR =  np.loadtxt(filename2a)
realization_infected_riskmatrix_ave1_SIR = np.loadtxt(filename3a)

mean_infected_riskmatrix1_SIR = np.mean(realization_infected_riskmatrix1_SIR,axis=0)/N 
mean_infected_riskmatrix_ave1_SIR = np.mean(realization_infected_riskmatrix_ave1_SIR,axis=0)/N 
mean_infected_matrix_network1_SIR = np.mean(realization_infected_matrix_network1_SIR,axis=0)/N


# Load txt file
filename4a = "SIR_ER_Realizationsinfected_Network_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
filename5a = "SIR_ER_Realizationsinfected_Riskmatrix_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
filename6a= "SIR_ER_Realizationsinfected_Riskmatrixave_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
 

# Displaying the contents of the text file
realization_infected_matrix_network2_SIR  = np.loadtxt(filename4a)
realization_infected_riskmatrix2_SIR =  np.loadtxt(filename5a)
realization_infected_riskmatrix_ave2_SIR = np.loadtxt(filename6a)


mean_infected_matrix_network2_SIR = np.mean(realization_infected_matrix_network2_SIR,axis=0)/N
mean_infected_riskmatrix2_SIR = np.mean(realization_infected_riskmatrix2_SIR,axis=0)/N 
mean_infected_riskmatrix_ave2_SIR = np.mean(realization_infected_riskmatrix_ave2_SIR,axis=0)/N 
#############SITAD


filename1b = "SITAD_BA_Realizationsinfected_Network_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
filename2b = "SITAD_BA_Realizationsinfected_riskmatrix_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
filename3b = "SITAD_BA_Realizationsinfected_avematrix_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

# Displaying the contents of the text file
realization_infected_matrix_network1_SITAD = np.loadtxt(filename1b)
realization_infected_riskmatrix1_SITAD =  np.loadtxt(filename2b)
realization_infected_riskmatrix_ave1_SITAD = np.loadtxt(filename3b)



mean_infected_riskmatrix1_SITAD = np.mean(realization_infected_riskmatrix1_SITAD,axis=0)/N 
mean_infected_riskmatrix_ave1_SITAD = np.mean(realization_infected_riskmatrix_ave1_SITAD,axis=0)/N 
mean_infected_matrix_network1_SITAD = np.mean(realization_infected_matrix_network1_SITAD,axis=0)/N

##################
########################
##############ER saving area###################################

filename4b = "SITAD_ER_Realizationsinfected_Network_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
filename5b = "SITAD_ER_Realizationsinfected_riskmatrix_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"
filename6b = "SITAD_ER_Realizationsinfected_avematrix_" + str(N) + "nodes"+ str(nrepeat)+"rep.txt"

# Displaying the contents of the text file
realization_infected_matrix_network2_SITAD = np.loadtxt(filename4b)
realization_infected_riskmatrix2_SITAD =  np.loadtxt(filename5b)
realization_infected_riskmatrix_ave2_SITAD = np.loadtxt(filename6b)


mean_infected_matrix_network2_SITAD = np.mean(realization_infected_matrix_network2_SITAD,axis=0)/N
mean_infected_riskmatrix2_SITAD = np.mean(realization_infected_riskmatrix2_SITAD,axis=0)/N 
mean_infected_riskmatrix_ave2_SITAD = np.mean(realization_infected_riskmatrix_ave2_SITAD,axis=0)/N 


mean_infected_matrix_network1_SITAD = mean_infected_matrix_network1_SITAD[0:100]
mean_infected_riskmatrix1_SITAD = mean_infected_riskmatrix1_SITAD[0:100]
mean_infected_riskmatrix_ave1_SITAD = mean_infected_riskmatrix_ave1_SITAD[0:100]

mean_infected_matrix_network2_SITAD = mean_infected_matrix_network2_SITAD[0:100]
mean_infected_riskmatrix2_SITAD = mean_infected_riskmatrix2_SITAD[0:100] 
mean_infected_riskmatrix_ave2_SITAD = mean_infected_riskmatrix_ave2_SITAD[0:100]

len(mean_infected_riskmatrix_ave2_SITAD)
##############




#plt.subplot(1, 2, 1)
#fig, ((axs1, axs2), (axs3, axs4)) = plt.subplots(nrows=2, ncols=2)
#fig, ((axs1,axs2),(axs3,axs4),(axs5,axs6)) = plt.subplots(nrows=3, ncols=2,
                       #constrained_layout = True)


fig, ((axs1,axs2),(axs3,axs4),(axs5,axs6)) = plt.subplots(nrows=3, ncols=2)
##SI process
h1,  = axs1.plot (mean_infected_riskmatrix1_SI,"-x", linewidth=1)
h2,  = axs1.plot (mean_infected_riskmatrix_ave1_SI,"-+", linewidth=1)
h3,  = axs1.plot (mean_infected_matrix_network1_SI,"-", linewidth=1)
axs1.legend ([h1,h2,h3], [ "Risk matrix" , "Average risk matrix", "Network"],
            prop={'size': 7})


#axs1.set_ylabel ("Infected proportion")
axs1.get_legend().remove()
xname = list(range(0,timesteps+1,20))
#xname = list(range(0,timesteps,20))
axs1.set_xticks(xname)
axs1.axes.xaxis.set_ticklabels([])
################
#plt.subplot(1, 2, 2)
h4,  = axs2.plot (mean_infected_riskmatrix2_SI,"-x", linewidth=1)
h5,  = axs2.plot (mean_infected_riskmatrix_ave2_SI,"-+", linewidth=1)
h6,  = axs2.plot (mean_infected_matrix_network2_SI,"-", linewidth=1)
axs2.legend ([h4,h5,h6], [ "Risk matrix" , "Average risk matrix", "Network"],loc ='lower right',
            prop={'size': 7})



axs2.get_legend().remove()
axs2.set_xticks(xname)
axs2.axes.yaxis.set_ticklabels([])
axs2.axes.xaxis.set_ticklabels([])
#fig.tight_layout(pad=0.3)
##############Theta change to .7 and Time step is 25
##SIR process
h1a,  = axs3.plot (mean_infected_riskmatrix1_SIR,"-x", linewidth=1)
h2a,  = axs3.plot (mean_infected_riskmatrix_ave1_SIR,"-+", linewidth=1)
h3a,  = axs3.plot (mean_infected_matrix_network1_SIR,"-", linewidth=1)
axs3.legend ([h1a,h2a,h3a], [ "Risk matrix" , "Average risk matrix", "Network"],
            prop={'size': 7})


#axs3.set_ylabel ("Infected proportion")
axs3.get_legend().remove()
axs3.set_xticks(xname)
axs3.axes.xaxis.set_ticklabels([])
h4a,  = axs4.plot (mean_infected_riskmatrix2_SIR,"-x", linewidth=1)
h5a,  = axs4.plot (mean_infected_riskmatrix_ave2_SIR,"-+", linewidth=1)
h6a,  = axs4.plot (mean_infected_matrix_network2_SIR,"-", linewidth=1)
axs4.legend ([h4a,h5a,h6a], [ "Risk matrix" , "Average risk matrix", "Network"],loc ='lower right',
            prop={'size': 5})



axs4.set_xticks(xname)
axs4.axes.yaxis.set_ticklabels([])
axs4.get_legend().remove()
axs4.axes.xaxis.set_ticklabels([])
#####SITAD
h1b,  = axs5.plot (mean_infected_riskmatrix1_SITAD,"-x", linewidth=1)
h2b,  = axs5.plot (mean_infected_riskmatrix_ave1_SITAD,"-+", linewidth=1)
h3b,  = axs5.plot (mean_infected_matrix_network1_SITAD,"-", linewidth=1)
axs5.legend ([h1b,h2b,h3b], [ "Risk matrix" , "Average risk matrix", "Network"],
            prop={'size': 7})

#axs5.set_xlabel ("Time steps")
#axs5.set_ylabel ("Infected proportion")
axs5.get_legend().remove()
axs5.set_xticks(xname)

h4b,  = axs6.plot (mean_infected_riskmatrix2_SITAD,"-x", linewidth=1)
h5b,  = axs6.plot (mean_infected_riskmatrix_ave2_SITAD,"-+", linewidth=1)
h6b,  = axs6.plot (mean_infected_matrix_network2_SITAD,"-", linewidth=1)
axs6.legend ([h4b,h5b,h6b], [ "Modified" , "ATMM", "Proposed"],loc ='lower right',
            prop={'size': 6})


#axs6.set_xlabel ("Time steps")
axs6.set_xticks(xname)
axs6.axes.yaxis.set_ticklabels([])

fig.text(0.01, 0.5, 'Infected proportion', ha='center', va='center', rotation='vertical')

fig.text(0.555, 0.012, 'Time', ha='center', va='center')
#fig.tight_layout(pad=3)
fig.tight_layout(pad=0.5, w_pad=0.1, h_pad=.35)
###############

plot_name = "ApproximateCompare.pdf"

plt.savefig(plot_dir + plot_name)

###############



plt.show()
# #################

