import os
import numpy as np

os.chdir("C:\\Users\\pwr844\\OneDrive - University of Tennessee\\Desktop\\2.Connecting\\Realdata\\fold1\\")
path = "C:\\Users\\pwr844\\OneDrive - University of Tennessee\\Desktop\\2.Connecting\\Realdata\\fold1\\"

os.getcwd()
Ns = [100,300]#change to 100, 1000, PARTICLES
Ts = [50,60] #change to 30,50, FINAL THRESHOLD
Nrep1 = 1 # change to 101
Nrep2 = 101
#define a cover function, to check if x belongs to (a,b)
def cover_func(a,b,x):
    u = (a-x)*(b-x)
    if u<=0:
        return(1)
    else:
        return(0)
    
cover_func(100,200,122)
size = len(Ns)*len(Ts)*(Nrep2-Nrep1)

##function to obtain outputs for [time, beta_cover_IQR, beta_cover_CI95, beta_IQR ,
#  gamma_cover_IQR, gamma_cover_CI95, gamma_IQR, number files exist]

#1st row for exact
#2nd row for average     
def ABCoutputfunc(particles, finalthreshold):
    matrixcollect_Exact = np.zeros((Nrep2-Nrep1,7))   #1st col 95% cover beta, 2nd 50% cover beta, #3rd col 95% cover gamma, 4th 50% cover gamma, 5th col time 
    matrixcollect_Ave = np.zeros((Nrep2-Nrep1,7))  
    outputmat = np.zeros((2,8))
    for i in range(Nrep1,Nrep2):
        para1 = particles
        para2 = finalthreshold
        realization = i
        
       
        filename1 = "ABCFold1_SIR_ExactNetworkSummaryEstimation_realization" + str(realization)+"particles_"+str(para1)+"thre_"+str(para2)+".txt"
        filename2 = "ABCFold1_SIR_AverageApproximateSummaryEstimation_realization" + str(realization)+"particles_"+str(para1)+"thre_"+str(para2)+".txt"
        
        if(os.path.isfile(filename1)):
            data1 =  np.loadtxt(filename1)
            beta_cover_CI95 = cover_func(data1[0,0],data1[1,0],data1[4,0])
            beta_cover_IQR = cover_func(data1[2,0],data1[3,0],data1[4,0])
            beta_IQR = data1[3,0]-data1[2,0]
            gamma_cover_CI95 = cover_func(data1[0,1],data1[1,1],data1[4,1])
            gamma_cover_IQR = cover_func(data1[2,1],data1[3,1],data1[4,1])
            gamma_IQR = data1[3,1]-data1[2,1]
            time = data1[4,2]
            matrixcollect_Exact[i-1] = [time, beta_cover_IQR, beta_cover_CI95, beta_IQR ,  gamma_cover_IQR, gamma_cover_CI95, gamma_IQR]
        
        if(os.path.isfile(filename2)):
            data2 =  np.loadtxt(filename2)
            beta_cover_CI95_2 = cover_func(data2[0,0],data2[1,0],data2[4,0])
            beta_cover_IQR_2 = cover_func(data2[2,0],data2[3,0],data2[4,0])
            beta_IQR_2 = data2[3,0]-data2[2,0]
            gamma_cover_CI95_2 = cover_func(data2[0,1],data2[1,1],data2[4,1])
            gamma_cover_IQR_2 = cover_func(data2[2,1],data2[3,1],data2[4,1])
            gamma_IQR_2 = data2[3,1]-data2[2,1]
            time_2 = data2[4,2]
            matrixcollect_Ave[i-1] = [time_2, beta_cover_IQR_2, beta_cover_CI95_2, beta_IQR_2 ,  gamma_cover_IQR_2, gamma_cover_CI95_2, gamma_IQR_2]
       

    indices_Exact = np.all( matrixcollect_Exact == 0, axis=1)
    matrixcollect_Exact = np.delete(matrixcollect_Exact, indices_Exact, axis = 0)
    
    indices_Ave = np.all( matrixcollect_Ave == 0, axis=1)
    matrixcollect_Ave = np.delete(matrixcollect_Ave, indices_Ave, axis = 0)
    
    coverrate_exact = np.mean(matrixcollect_Exact, axis=0)
    coverrate_ave = np.mean(matrixcollect_Ave, axis=0)
    outputmat[0,:] = np.append(coverrate_exact,len(matrixcollect_Exact))
    outputmat[1,:] = np.append(coverrate_ave, len(matrixcollect_Ave))
    return(outputmat)
        

            
particles = 100
finalthreshold = 40

myoutput = ABCoutputfunc(particles, finalthreshold)
myoutput

################ARRAY TO LATEX
def ABCoutputfunc_commonindices(particles, finalthreshold):
    matrixcollect_Exact = np.zeros((Nrep2-Nrep1,7))   #1st col 95% cover beta, 2nd 50% cover beta, #3rd col 95% cover gamma, 4th 50% cover gamma, 5th col time 
    matrixcollect_Ave = np.zeros((Nrep2-Nrep1,7))  
    outputmat = np.zeros((2,8))
    for i in range(Nrep1,Nrep2):
        para1 = particles
        para2 = finalthreshold
        realization = i
        
       
        filename1 = "ABCFold1_SIR_ExactNetworkSummaryEstimation_realization" + str(realization)+"particles_"+str(para1)+"thre_"+str(para2)+".txt"
        filename2 = "ABCFold1_SIR_AverageApproximateSummaryEstimation_realization" + str(realization)+"particles_"+str(para1)+"thre_"+str(para2)+".txt"
        
        if(os.path.isfile(filename1)):
            data1 =  np.loadtxt(filename1)
            beta_cover_CI95 = cover_func(data1[0,0],data1[1,0],data1[4,0])
            beta_cover_IQR = cover_func(data1[2,0],data1[3,0],data1[4,0])
            beta_IQR = data1[3,0]-data1[2,0]
            gamma_cover_CI95 = cover_func(data1[0,1],data1[1,1],data1[4,1])
            gamma_cover_IQR = cover_func(data1[2,1],data1[3,1],data1[4,1])
            gamma_IQR = data1[3,1]-data1[2,1]
            time = data1[4,2]
            matrixcollect_Exact[i-1] = [time, beta_cover_IQR, beta_cover_CI95, beta_IQR ,  gamma_cover_IQR, gamma_cover_CI95, gamma_IQR]
        
        if(os.path.isfile(filename2)):
            data2 =  np.loadtxt(filename2)
            beta_cover_CI95_2 = cover_func(data2[0,0],data2[1,0],data2[4,0])
            beta_cover_IQR_2 = cover_func(data2[2,0],data2[3,0],data2[4,0])
            beta_IQR_2 = data2[3,0]-data2[2,0]
            gamma_cover_CI95_2 = cover_func(data2[0,1],data2[1,1],data2[4,1])
            gamma_cover_IQR_2 = cover_func(data2[2,1],data2[3,1],data2[4,1])
            gamma_IQR_2 = data2[3,1]-data2[2,1]
            time_2 = data2[4,2]
            matrixcollect_Ave[i-1] = [time_2, beta_cover_IQR_2, beta_cover_CI95_2, beta_IQR_2 ,  gamma_cover_IQR_2, gamma_cover_CI95_2, gamma_IQR_2]
       

    indices_Exact = np.all( matrixcollect_Exact == 0, axis=1)
    indices_Ave = np.all( matrixcollect_Ave == 0, axis=1)
    indices_Common1 = np.logical_or(indices_Exact,indices_Ave)
    indices_Common = np.array(indices_Common1)

    
    
    matrixcollect_Exact = np.delete(matrixcollect_Exact, indices_Common, axis = 0)
    
    
    matrixcollect_Ave = np.delete(matrixcollect_Ave, indices_Common, axis = 0)
    
    coverrate_exact = np.mean(matrixcollect_Exact, axis=0)
    coverrate_ave = np.mean(matrixcollect_Ave, axis=0)
    outputmat[0,:] = np.append(coverrate_exact,len(matrixcollect_Exact))
    outputmat[1,:] = np.append(coverrate_ave, len(matrixcollect_Ave))
    return(outputmat)
        


################################


particles = 100
finalthreshold = 40

myoutput = ABCoutputfunc(particles, finalthreshold)
myoutput

import array_to_latex as a2l

latex_code = a2l.to_ltx(myoutput, frmt = '{:6.3f}', arraytype = 'bmatrix', print_out=False)

latex_code

myoutput1 = ABCoutputfunc_commonindices(particles, finalthreshold)
myoutput1

import array_to_latex as a2l

latex_code1 = a2l.to_ltx(myoutput1, frmt = '{:6.3f}', arraytype = 'bmatrix', print_out=False)

latex_code1





2833/(36)
777/(36)
