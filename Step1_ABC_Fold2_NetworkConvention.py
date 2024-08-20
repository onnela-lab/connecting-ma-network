
import os
os.chdir("C:\\Users\\pwr844\\OneDrive - University of Tennessee\\Desktop\\2.Connecting\\Realdata\\fold2\\")
#os.getcwd()
import time
import sys
#para1 = int(sys.argv[1]) #NUMBER PARTICLES
#para2 = int(sys.argv[2]) #THRESHOLD FINAL
#para3 = int(sys.argv[3]) #REALIZATIONS
para1 = 100
para2 =  75
para3 = 1

import networkx as nx
from networkx.algorithms import approximation
import random
import scipy.stats as ss
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.random import choice
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors, KNeighborsRegressor

from scipy.stats import rv_discrete, rv_continuous

from Avetransmissionmatrix_function import Avetransmissionmatrix_function, Transmissionmatrix_function
##Network Working on
start_time = time.time()
filename1 = "CNScombined_adjacencymatrix.txt"
CNSMatrix =  np.loadtxt(filename1)
N = len(CNSMatrix) # Network size

##Underlying network
G = nx.from_numpy_array(CNSMatrix) #make graph G from matrix A


#underlying truth need to recover, need from step 1
realization = para3

     
Ini1 = 0
initial_infectednodes = [Ini1]
initial_status = [N-len(initial_infectednodes),len(initial_infectednodes),0]
filename = "SIR_SyntheticEpidata_realization" + str(realization)+".txt"
data = np.loadtxt(filename)
data_obs = data[:,0:3]
timesteps = len(data_obs)

######### ABC set up

alpha = 0.1 # percentage of particles replaced at each step
scale_factor = 2
num_acc_sim = para1
threshold_init = 100000
threshold_final = para2 #make it .001


def compute_summaries(data):
    """ Compute network features, computational times and their nature.            
    """
    dictSums = dict()   # Will store the summary statistic values
                        
    # Daily confirmed and recovered
    IR = pd.DataFrame(data[:,1]+data[:,2])
    R = pd.DataFrame(data[:,2]) 
    dictSums["Infected_daily"] = IR.diff()[1:len(IR)] 
    dictSums["Recovered_daily"] = R.diff()[1:len(R)] 
    resDicts = dictSums
    
    return resDicts

dict_summaries_obs = compute_summaries(data_obs)

df_summaries_obs = pd.DataFrame([dict_summaries_obs])
###################





##define SIR_Spreadfunction using network convention
def SIR_Spreadfunction_Network_Convention(G, beta,gamma, timesteps, initial_infectednodes):
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
    i_nodes = list(i_nodes1) # assign infected nodes
    r_nodes = list()
    s_nodes = list ( set (G.nodes ()) - set (i_nodes ) - set (r_nodes ) )
   
    for node_order in range(len(i_nodes)):
        infected_status_matrix[i_nodes[node_order],0] = 1 # assign infected nodes to the first time step of the infected matrix
       
# Loop over all time steps .
    for step in range(1, timesteps):
               
        #new infections at each step and update s_nodes, i_nodes
             
        new_infected_nodes = set()
        for node in i_nodes :
            
            if G.degree(node) > 0:
               node_neighbors = G.neighbors( node)

               for node1 in  node_neighbors:
                   h1  = random.random ()
                   if (node1 in s_nodes) and (h1 < beta):
                      new_infected_nodes.add(node1)
                     
         
        #########Step 2, Second update recovered node
        if len(i_nodes)==0:
            new_recovered_nodes = []
        else:
            new_recovered = np.random.binomial(len(i_nodes), gamma, 1) # be binomial for a closer approximate for arbitrary network
    
            new_recovered_nodes = choice(i_nodes, new_recovered, replace=False) 
        #update r nodes and i nodes
        
  
        
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


######Stop here at 11h26, Feb 1
#define summary stat

### Compute the summaries on a simulated data

def _merge_dict(dict1, dict2):
    """
    Function to merge two dictionaries
    """
    res = {**dict1, **dict2}
    return res


# Start by defining the function to generate one data
def data_indiv_simulation(model, prior_args_model = None,fixed_args_model = None):
   
    
    
    if prior_args_model is None:
        prior_args_model = dict()
   
    if fixed_args_model is None:
        fixed_args_model = dict()
    
    
 
    # Generate the mechanism parameters from the priors
    dict_params = dict()
    for (key, value) in prior_args_model.items():
        if isinstance(value, ss._distn_infrastructure.rv_frozen):
            dict_params[key] = value.rvs(random_state=ss.randint(0,4294967296).rvs()) #2**32 = 4294967296
        else:
            raise ValueError('Invalid value for the prior simulation object.')

    # Fuse the dictionaries of simulated and fixed parameters together
    args_model = _merge_dict(dict_params, fixed_args_model)
    
    # Simulate from the model.
    data_sim = model( **args_model)
    
    # Compute the summary statistics on the simulated data
    dict_summaries = compute_summaries(data_sim)
    
    return dict_params, dict_summaries


# -- Example of its use


# Define the prior on the parameters of the DMC model
beta_prior = ss.uniform(loc = 0, scale = 1) # for a uniform prior between loc + scale*U(0,1)
gamma_prior = ss.uniform(loc = 0, scale = .5)
prior_args_model = {"beta":beta_prior, "gamma":gamma_prior}
prior_supports = [[0,1],[0,0.5]]
# Define the fixed arguments of the model



###########

##### Implementation of an SMC-ABC algorithm #####

def _perturb_continuous_param_on_support(prior_support, perturb_kernel):
    """ Perturb a continuous parameter thanks to a truncated Gaussian distribution """
    
    # Generate the perturbed value
    perturbed_float = perturb_kernel.rvs()
    while perturbed_float < prior_support[0] or perturbed_float > prior_support[1]:
        perturbed_float = perturb_kernel.rvs()

    return perturbed_float
####################

def distance_func(df_sim_summaries, df_obs_summaries):   
    dist1 = np.sum(np.array( (df_sim_summaries.iloc[0,0] - df_obs_summaries.iloc[0,0])**2 ))   
    dist2 = np.sum(np.array( (df_sim_summaries.iloc[0,1] - df_obs_summaries.iloc[0,1])**2 ))   
    dist = (dist1  + dist2)**.5
    return dist





##########
def abc_RSMCABC(model,prior_supports = prior_supports,
                prior_args_model = None, fixed_args_model = None,
                threshold_init = threshold_init, threshold_final = threshold_final,
                alpha = 0.1, scale_factor = 2,
                perturbation = "Gaussian",
                num_acc_sim = num_acc_sim, df_observed_summaries = None,
                distance_func = distance_func):
    
    """ Implementation of the replenishment SMC ABC algorithm.
    """
         
    # Recover a list for the parameter priors
    list_priors = []
    # For each key, we recover from the prior
    for (key, value) in prior_args_model.items():
        list_priors += [value]
        
    
    # Uniform distribution for MH test
    unif_dist = ss.uniform(loc = 0, scale = 1)
    
    # Number of particles to discard at each step
    num_drop_sim = int(alpha * num_acc_sim)
    if num_drop_sim == 0:
        num_drop_sim = 1
    
    # Identify the summary statistics to keep while simulating
    #cols_to_keep = df_observed_summaries.columns
    
    step_count = 0      # number of sequential steps
    sim_count_total = 0 # total number of simulated data
    
    # To store accepted weights/parameters values, and distances
    df_params = pd.DataFrame()
    df_dist_acc = pd.DataFrame()
    
    # Keep track of the epsilon values
    epsilon_values = [threshold_init]
    
    if step_count == 0:
        
        sim_count = 0 # number of accepted simulations during the current step
        
        ### Initial classic rejection sampling algorithm
        while sim_count < num_acc_sim:
            
            sim_count_total += 1
            
            ### Simulate the parameters
                        
            # Generate parameters from the priors
            sim_values_args_model = dict()
            
            # For each key, we simulate from the prior
            for (key, value) in prior_args_model.items():
                if isinstance(value, ss._distn_infrastructure.rv_frozen):
                    sim_values_args_model[key] = value.rvs()
                else:
                    raise ValueError('Invalid specified value for the parameter prior distribution. Each prior must be a scipy.stats._distn_infrastructure.rv_frozen object.')
            
            # Fuse the dictionaries of simulated and fixed parameters together
            args_model = _merge_dict(sim_values_args_model,
                                          fixed_args_model)
                                
            
            ### Simulate from the model given the simulated params
            data_sim = model(**args_model)    

            # Compute the summary statistics on the simulated data
            dict_summaries = compute_summaries(data_sim)
        
            ### Format the parameter values simulated # could be removed in version for Thien
            dict_params = {}
            for (key, value) in sim_values_args_model.items():
                    dict_params[key] = value

            ### Convert the simulated network summaries and reduce it if necessary
            df_summaries = pd.DataFrame([dict_summaries])
            #df_summaries_reduced = df_summaries[cols_to_keep]

            ### Compute the distance between simulated and observed data summarized
            dist = distance_func(df_summaries, df_observed_summaries)
            
            #dist = distance_func(df_summaries_reduced, df_observed_summaries)
            
            # If the distance is low enough, accept the simulated parameters
            if dist <= threshold_init:
    
                # store also the accepted weights/parameters and resulting distance
               
                #df_params = df_params.append(pd.DataFrame([dict_params]), ignore_index=True)
                #df_dist_acc = df_dist_acc.append(pd.DataFrame([dist]), ignore_index=True)
                df_params = pd.concat([df_params, pd.DataFrame([dict_params])], ignore_index=True)
                
                df_dist_acc = pd.concat([df_dist_acc, pd.DataFrame([dist])], ignore_index=True)
                
                sim_count += 1
        
        step_count += 1
    
    # SMC-ABC core part
    if step_count > 0:
        
        # Determine the order of the distances when sorted in increasing order
        idx_sort = np.argsort(df_dist_acc.iloc[:,0])
        
        # Reorder the parameters and distance with this order
        df_dist_acc = df_dist_acc.iloc[idx_sort,:]
        df_dist_acc = df_dist_acc.reset_index(drop=True)
                
        df_params = df_params.iloc[idx_sort,:]
        df_params = df_params.reset_index(drop=True)
                
        # Compute epsilon_max = the maximal distance
        epsilon_max = df_dist_acc.iloc[num_acc_sim-1,0]
        
        epsilon_values = epsilon_values + [epsilon_max]
        
        # while epsilon_max is greater than epsilon_final
        while epsilon_max > threshold_final:
            
            print(epsilon_max, threshold_final)
            
            # Drop the num_drop_sim (Na) particles with largest distances
            df_params.drop(df_dist_acc.tail(num_drop_sim).index, inplace=True)
            df_dist_acc.drop(df_dist_acc.tail(num_drop_sim).index, inplace=True)
            
            epsilon_next = df_dist_acc.tail(1).iloc[0,0] # the largest distance of the remaining simulations
                        
            std_params = scale_factor * df_params.apply(np.std)
            
            ### Resample num_drop_sim new particles and data that are accepted
            
            num_acc_next = 0
            
            while num_acc_next < num_drop_sim:
            
            #for j in range(num_drop_sim):
                
                ### Sample an old parameter value from the
                ### num_acc_sim - num_drop_sim previously accepted values 
                idx_sel = np.random.choice(df_params.index[:(num_acc_sim-num_drop_sim)])
                sim_count_total += 1
                
                ### Perturb the selected parameter values with a kernel
                              
                # Parameter perturbation
                prev_params = np.array(df_params.iloc[idx_sel,:])
                
                perturbed_params = np.empty(len(prev_params))
                # For each parameter value
                for i in range(len(prev_params)):
                    perturbation_kernel_Gauss = ss.norm(prev_params[i], std_params[i])
                    perturbed_params[i] = _perturb_continuous_param_on_support(prior_supports[i], perturbation_kernel_Gauss)
                    
                # To use the simulated parameters in our data generation function
                # we need a list of dict, with same structure as sim_args_model
                perturbed_params_dict = dict()
                idx_params = 0
                for (key, value) in sim_values_args_model.items():
                    perturbed_params_dict[key] = perturbed_params[idx_params]
                    idx_params += 1
                
                ### Generate a new data given the perturbed parameters
                args_model = _merge_dict(perturbed_params_dict,
                                         fixed_args_model)

                data_sim = model(**args_model)
                dict_summaries = compute_summaries(data_sim)
                
                df_summaries = pd.DataFrame([dict_summaries])
                #df_summaries_reduced = df_summaries[cols_to_keep]
    
                dist_new = distance_func(df_summaries,
                                         df_observed_summaries)
    
                #dist_new = distance_func(df_summaries_reduced,
                   #                      df_observed_summaries)

                if dist_new <= epsilon_next:
                    
                    print("Dist_new: ", dist_new, " Epsilon next: ", epsilon_next)
                                                            
                    # For the parameters
                    list_prior_params_old = []
                    list_prior_params_new = []
                    list_pdf_new_given_old = []
                    list_pdf_old_given_new = []
                    for i in range(len(list_priors)):
                        list_prior_params_old += [list_priors[i].pdf(prev_params[i])]
                        list_prior_params_new += [list_priors[i].pdf(perturbed_params[i])]
                        list_pdf_new_given_old += [ss.norm(prev_params[i], std_params[i]).pdf(perturbed_params[i])]
                        list_pdf_old_given_new += [ss.norm(perturbed_params[i], std_params[i]).pdf(prev_params[i])]

                    prior_ratio_params = np.prod(list_prior_params_new) / np.prod(list_prior_params_old)
                    transition_ratio_params = np.prod(list_pdf_old_given_new) / np.prod(list_pdf_new_given_old)
                    
                    mh_ratio = np.min([1, prior_ratio_params * transition_ratio_params])
                    
                    if unif_dist.rvs() < mh_ratio:
                            
                        if len(perturbed_params) > 0:
                            perturbed_params_df = pd.DataFrame(perturbed_params.reshape(-1, len(perturbed_params)),columns=df_params.columns)
                            #df_params = df_params.append(perturbed_params_df, ignore_index=True)
                            df_params = pd.concat([df_params,perturbed_params_df], ignore_index=True)
                        else:
                            df_params = pd.concat([df_params, pd.DataFrame([], index=[1])], ignore_index=True)
                        df_dist_acc = pd.concat([df_dist_acc, pd.DataFrame([dist_new])], ignore_index=True)
                        
                        num_acc_next += 1
                        
            # Determine the order of the distances when sorted in increasing order
            idx_sort = np.argsort(df_dist_acc.iloc[:,0])
            
            # Reorder the parameters and distance with this order
            df_dist_acc = df_dist_acc.iloc[idx_sort,:]
            df_dist_acc = df_dist_acc.reset_index(drop=True)
                        
            df_params = df_params.iloc[idx_sort,:]
            df_params = df_params.reset_index(drop=True)
            
            # Compute epsilon_max = the maximal distance
            epsilon_max = df_dist_acc.iloc[num_acc_sim-1,0]
            
            epsilon_values = epsilon_values + [epsilon_max]
            
            step_count += 1
            
        threshold_values = np.array(epsilon_values)

        return df_params, df_dist_acc, sim_count_total, threshold_values
    


###Use modified version, Mymodel1 with average of 30 random realizations###############

Mymodel_start_time = time.time()


Mymodel1 =  SIR_Spreadfunction_Network_Convention

initial_status = [N-1,1,0]
fixed_args_model = {"G":G,"timesteps":timesteps, "initial_infectednodes":initial_infectednodes}    


(df_params_RABC_obs, 
 df_dist_acc_RABC_obs, 
 sim_count_total_RABC_obs, 
 threshold_values_RABC_obs) = abc_RSMCABC(model = Mymodel1,prior_supports = prior_supports,
                                          prior_args_model = prior_args_model, 
                                          fixed_args_model = fixed_args_model,
                                          threshold_init = threshold_init, threshold_final = threshold_final,
                                          alpha = alpha, scale_factor = scale_factor,
                                          perturbation = "Gaussian",
                                          num_acc_sim = num_acc_sim,
                                          df_observed_summaries = df_summaries_obs,
                                          distance_func = distance_func)
                                          

        
# Estimate the parameters
#underlyingtruth
#filename1 = "ABCFold1_SIR_AverageApproximateParamterEstimation_realization" + str(realization)+"particles_"+str(para1)+"thre_"+str(para2)+".txt"


#np.savetxt(filename1,df_params_RABC_obs)
df_observed_summaries = df_summaries_obs
dist_vector = np.zeros(len(df_params_RABC_obs))
currentI_data = np.zeros((len(df_params_RABC_obs), len(data))) #a matrix of 100 cols, each col with length of data
totalI_data = np.zeros((len(df_params_RABC_obs), len(data))) 




for i in range(len(df_params_RABC_obs)):
    beta = df_params_RABC_obs.iloc[i,0]
    gamma = df_params_RABC_obs.iloc[i,1]
    initial_status = [N-1,1,0]
    data = SIR_Spreadfunction_Network_Convention(G, beta,gamma, timesteps, initial_infectednodes)
    dict_summaries = compute_summaries(data)
    df_summaries = pd.DataFrame([dict_summaries])
    dist = distance_func(df_summaries, df_observed_summaries)    
    dist_vector[i] = dist
    currentI_data[:,i] = data[:,1]
    totalI_data[:,i] = data[:,1] + data[:,2]
    
dist_vector1 = pd.DataFrame(dist_vector)  
observedkeep = 30
    
idx_sort30 = np.argsort(dist_vector1.iloc[:,0])[0:observedkeep]

# Reorder the parameters and distance with this order
dist_vector_sort30 = dist_vector1.iloc[idx_sort30]    
currentI_data_sort30 = currentI_data[:,idx_sort30]
totalI_data_sort30 = totalI_data[:,idx_sort30]


currentI_mean = np.mean(currentI_data_sort30, axis=1)

totalI_mean = np.mean(totalI_data_sort30, axis=1)

currentI_obs = data_obs[:,1]
totalI_obs = data_obs[:,1] + data_obs[:,2] 
########Plotting


x = range(0,timesteps)
y1 = currentI_mean
y2 = currentI_obs
CI_currentI_mean = 1.96*np.std(currentI_data_sort30, axis=1)/(observedkeep**.5)



plt.fill_between(x, (y1-CI_currentI_mean), (y1+CI_currentI_mean), color='b', alpha=.1)

plt.plot(x, y1, label = "Convention Predicted Current Infected")
plt.plot(x, y2, '-.',  label = "Observed Current Infected")

plt.xlabel("Time")
plt.ylabel("Current Infected")
plt.legend(loc="lower right")
plt.savefig('Fold2_NetworkConvention_CurrentInfectedFit.pdf')
plt.show()




##############TOTAL INFECTED#########

y3 = totalI_mean

y4 = totalI_obs
CI_totalI_mean = 1.96*np.std(totalI_data_sort30, axis=1)/(observedkeep**.5)
plt.fill_between(x, (y3-CI_totalI_mean), (y3+CI_totalI_mean), color='b', alpha=.1)

plt.plot(x, y3, label = "Convention Predicted Total Infected")
plt.plot(x, y4, '-.',  label = "Observed Total Infected")

plt.xlabel("Time")
plt.ylabel("Total Infected")
plt.legend(loc="lower right")
plt.savefig('Fold2_NetworkConvention_TotalInfectedFit.pdf')
plt.show()
avenet_runtime = time.time() - start_time
filename = "ABCFold2_SIR_NetworkConventionSummary_realization" + str(realization)+"particles_"+str(para1)+"thre_"+str(para2)+".txt"
mean_dist_sort30 = np.mean(dist_vector_sort30)
std_dist_sort30 = np.std(dist_vector_sort30)[0]
summary_vals = [mean_dist_sort30,std_dist_sort30,avenet_runtime] 
np.savetxt(filename,summary_vals)

