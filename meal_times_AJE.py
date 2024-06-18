# -*- coding: utf-8 -*-
"""
This script outlines the code used to produce the figures in Case study 1 in the manuscript
 "The Wells-Riley Model Revisited: Randomness, Heterogeneity and Transient Behaviours"; 
 Alexander J. Edwards, Marco-Felipe King, Catherine J. Noakes, Daniel Peckham, Martin Lopez-Garcia.


This data used in this script is taken from Table 1 in https://doi.org/10.1016/S0195-6663(03)00109-0 .


CREATED 27/06/2023 AJE
"""
import matplotlib.pyplot as plt
import  math
import numpy as np
from scipy.stats import erlang #for the erland distrbution
import sympy as sym
from round_up import rounding
import seaborn as sns
import pandas


################### DATA ###############################
#The data here is taken from Table 1 in https://doi.org/10.1016/S0195-6663(03)00109-0. 
#The mean and standard deviation values are given in minutes, so for consistency - all parameter values will be standardised to minutes

#values from data
###### CAFETERIA###
mu_cafe = [12.6, 23.0, 33.0, 41.1, 44.0]
sd_cafe = [3.8, 7.9, 11.3, 10.6, 14.2]
tables_cafe = [24, 34, 28, 41, 21]

##### Restaurant
mu_restaurant = [27.6, 44.9, 47.2, 52.3, 58.5] #MEAN
sd_restaurant = [6.7, 10.8, 10.1, 8.5, 13.1] #STANDARD DEVIATION
tables_restaurant = [8, 29, 13, 24, 21] # the total number of tables used to give the results 

##### Fast-food Restaurant
mu_fastfood = [10.7, 18.2, 18.4, 19.7, 21.9]
sd_fastfood = [3.3, 6.0, 6.8, 7.2, 5.8]
tables_fastfood = [22, 33, 23, 28, 18]

###########

######### combining the data sets
n_data = 3
set_type = ['Workplace Cafeteria', 'Restaurant', 'Fast-food Restaurant']
n_groups = 5 #number of data sets in each setting

mu = [mu_cafe, mu_restaurant, mu_fastfood] #MEAN for each setting
sd = [sd_cafe, sd_restaurant, sd_fastfood] #STANDARD DEVIATION for each setting
group_size = ['1', '2', '3', '4', '5+'] #number of people in each groups
group_size_num = [1,2,3,4,5]
#############################################################
#############################################################
#############################################################

##################### CALCULATING THE ERLANGS ###############
#initialise erlang distrbution values to store each one
k_erl = np.empty((n_data,n_groups))
lam_erl = np.empty((n_data,n_groups))

###########Sampling 100 times from each erlang ##########
times_erl_cafe = np.empty((100,n_groups))
times_erl_restaurant = np.empty((100,n_groups))
times_erl_fastfood = np.empty((100,n_groups))
times_erl=[times_erl_cafe,times_erl_restaurant,times_erl_fastfood]

#################################
#######FOR SANITY CHECK#######
####### check area under curve sums to 1
area = np.empty((n_data,n_groups))
#################################


colour = ['tab:blue', 'tab:orange', 'tab:purple', 'tab:red']


############ CALUCLATING ERLANGS
###Calculating the parameter values from the erlang distribution
#plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(n_data):
    plt.figure(dpi=750)#set dots per inch for better quality images
    #to include the party size=1 case, run from 0 insetad of 1: 
    for i in range(1,n_groups):#note we run from 1 here, to exlcude the group with n=1. This doesn't allow for a susceptible and an infectious indv.
        
        #solves simulateous equations TO FIND ERLANG PARAMETERS
        k,lam = sym.symbols('k,lam')
        eq1 = sym.Eq(k-mu[j][i]*lam,0)
        eq2 = sym.Eq(k - (lam*sd[j][i])**2,0)
        result = sym.solve([eq1,eq2],k,lam)
        print("Results from simulatneous equations = %s" %result)
        #####
        #converts results to floats
        k = float(result[1][0])
        #rounds k as k should be an integer for erlang
        k_erl[j,i] = rounding(k, 0)
        lam_erl[j,i] = float(result[1][1])
    
        #print values used for distribtion
        print("Distirbution parameters k, lam = %s, %s" %(k_erl[j,i], lam_erl[j,i]))
        
        ###### PLOTTING ######
        #Plot the pdf for the given distribution
        x= np.linspace(0,100,1000) #x values for the distribtuoin
        pdf = erlang.pdf(x,k_erl[j,i],loc=0, scale=1/lam_erl[j,i])
        #plt.plot(x,pdf, label = 'k=%s, $\lambda$=%s, Group Size=%s' %(k_erl[j,i],round(lam_erl[j,i],4),group_size_num[i]))
        plt.plot(x,pdf, label = 'Party Size=%s' %(group_size_num[i]), color=colour[i-1])
        
        ############
        ###### Sampling from an Erlang for random times for end analysis
        times_erl[j][:,i] = erlang.rvs(k_erl[j,i],scale = 1/lam_erl[j,i], size=100)
        #######################
        
        
        #####################
        #####SANITY CHECK########
        ####checking area under curve sums to 1
        area[j,i] = np.trapz(pdf, x)
        ####################
        
        
    plt.legend(fontsize=10)
    plt.ylabel('$f_T(t)$')
    plt.xlabel('T [min]')
    if j==1:
        plt.yticks(np.arange(0,0.06,0.01))
    #plt.title(set_type[j])
    plt.yticks(np.arange(0,0.08,0.01))
    plt.ylim(0,0.075)
    plt.show()
#############################################################





#############################################################
#What difference do we get from rounding??
mu_erl = np.empty((n_data,n_groups)) #initialise empty matrix for the means produced from the erlang after rounding
sd_erl = np.empty((n_data,n_groups)) #initialise empty matrix for the sd produced from the erlang after rounding
mu_dif = np.empty((n_data,n_groups)) #initialise empty matrix for the difference in means produced from the erlang compared with data
sd_dif = np.empty((n_data,n_groups)) #initialise empty matrix for the difference in sd produced from the erlang compared with data
mu_reldif = np.empty((n_data,n_groups)) #initialise empty matrix for the relative difference in means produced from the erlang compared with data
sd_reldif = np.empty((n_data,n_groups)) #initialise empty matrix for the relative difference in sd produced from the erlang compared with data
for j in range(n_data):
    for i in range(n_groups):
        mu_erl[j,i] = k_erl[j,i]/lam_erl[j,i]
        sd_erl[j,i] = math.sqrt((k_erl[j,i]/(lam_erl[j,i]**2)))
        mu_dif[j,i] = mu[j][i] - mu_erl[j,i]
        sd_dif[j,i] = sd[j][i] - sd_erl[j,i]
        mu_reldif[j,i] = (mu_dif[j,i]/mu[j][i])*100
        sd_reldif[j,i] = (sd_dif[j,i]/sd[j][i])*100
    
print('Relative difference in Means mu_reldif = %s'%mu_reldif)
print('Relative difference in standard deviation sd_reldif =%s'%sd_reldif)
    
                
#############################################################
#############################################################



##################### CALCULATING PROBABILITIES FOR RANDOM TIME ##############
##### EACH OUTBREAK SCENARIO IS THE SAME SET-UP FOR ALL OF THE DIFFERENT SETTINGS

#Set up with initial parameters (based on previous models)
# the data is given as minutes so all parameters will be in minutes also
###########################Ventilation and room parameters ##################
#quanta rate =  quanta per min)
q=6
#pulmonary rate as volume/min)
p=0.01
#Ventilaation rate = m^3/min
# for 3ACH
#Q=[15, 15, 15]  #3ach[15,15,15]
#for 1ACH
Q=[5, 5, 5] 

#volume of indoor space v= m^3
#v=100
v=[300,300,300] #for the different voled[46.4,456,456] #volume based on papers - discussed in PhD working file

############################################################
##### NEW PARAMETERS FOR REMOVAL RATE 
#ranges define first and last, then values in bewteen at equal steps
r_i = 0.00525  # this is viral decay/inactivation rate - range is given as [0,0.63] per hour in paper 
r_d = 0.015  #  This is viral deposition rate - range is given as [0.3,1.5] per hour

##The REMOVAL RATE
R = [(Q[i] + (v[i]*r_d) + (v[i]*r_i)) for i in range(len(Q))] 


################################################
############################ INITIAL CONDITIONS ############################
#inital Population size is N
N= group_size_num
#initial exposed poulation is E0
E0=np.zeros(n_groups)
#inital infected is I
I0 = np.ones(n_groups)
#succeptible population is S
S0= N - E0 - I0



nE = [[[] for i in range(n_groups)] for j in range(n_data)] #each possible value for number of exposed in each group
prob = [[[] for i in range(n_groups)] for j in range(n_data)] #initialise vector to store probability
for aa in range(n_data):
    #the steady state concentration is calculated from dc/dt=0 in vdc/dt = qI0 - QC
    #This is given for each scenario, depending on the number of infectors - in this case they are all the same
    beta_t=p*q*I0/(R[aa]) #I=1
    for kk in range(n_groups):
        S = int(S0[kk])
        prob_x = np.zeros(S+1)
        nE[aa][kk] = list(range(0,S+1)) 
        for ii in range(S+1):
            sum_i =0
            for jj in range(ii+1):#calculating the probability of becoming eposed
                int_sum_i =((-1)**jj) * math.comb(ii,jj) * (((S-ii+jj)*beta_t[kk] + lam_erl[aa,kk])**(-k_erl[aa,kk]))
                sum_i = sum_i + int_sum_i
                #print(int_sum_i)
                #print(sum_i)
            #calculating the probability
            prob_x[ii] = math.comb(S,ii) * (lam_erl[aa,kk]**k_erl[aa,kk]) * sum_i    
        
        #assigning for each exposed possiblility
        prob[aa][kk] = prob_x
                
        cumulutive_prob_x = sum(prob[aa][kk][:])
        print("cumulative for data set %s, group %s = %s"%(aa+1,kk+1,cumulutive_prob_x))
        
        
avg = [np.zeros(n_groups) for j in range(n_data)] #initialise vector to store mean in each scenario
for i in range(n_data):
    for j in range(n_groups):
        S = int(S0[j])
        for k in range(S+1):
            avg[i][j] = (prob[i][j][k])*k + avg[i][j]
        


bar_loc = [1,2,0] #this is to change the order of the scenarios when plotting
for i in range(1,n_groups):
    plt.figure(dpi=750)#set dots per inch for better quality images   
    S = int(S0[i])
    #including a width vector to plot bar charts
    width = np.empty(S+1);
    width.fill(0.1)
    #colour = iter(cm.tab10(np.linspace(0, 1, 10)))    
    for j in [2,0,1]:#this is to change the order of the scenarios when plotting
        #c=next(colour)
        plt.bar(nE[j][i] + bar_loc[j]*width, prob[j][i], width, label = '%s -> mean =%s' %(set_type[j],round(avg[j][i],2)), color=colour[j] )#k=%s, $\lambda$=%s' %(set_type[j], k_erl[j,i],round(lam_erl[j,i],4)), color=c)
        plt.vlines(avg[j][i], 0, 1,color=colour[j])
    plt.xticks(nE[j][i] + width, nE[j][i])
    plt.yticks(np.arange(0,1.1,0.1))
    #plt.title('Group size = %s' %(group_size[i]))
    plt.xlabel('Exposures [people]')
    plt.ylabel('Probability')
    plt.legend( fontsize=8, title='Scenario') 
    plt.ylim(0,1)
    plt.show()



###########################################################################
########## Sampling times for violin plot and per capita analysis##########
#######################################################################
#calculate the wellsriley per capita probability for each time, then
#compare with erlang per capita for unknown time
#also plot the sampled times as a swarmplot

#set up some parameters for use - pick lowest vent rate, to give highest prob
#quanta rate =  quanta per min)
q=6
#pulmonary rate as volume/min)
p=0.01
#Ventilaation rate = m^3/min
#for 1ACH
Q=5
#volume of indoor space v= m^3
v=300 
############################################################
##### NEW PARAMETERS FOR REMOVAL RATE 
#ranges define first and last, then values in bewteen at equal steps
r_i = 0.00525  # this is viral decay/inactivation rate - range is given as [0,0.63] per hour in paper 
r_d = 0.015  #  This is viral deposition rate - range is given as [0.3,1.5] per hour
##The REMOVAL RATE
R = Q + (v*r_d) + (v*r_i) 

######Wells-Riley Per Capita##############################
WR_percap_cafe = np.empty((100,n_groups))
WR_percap_restaurnt = np.empty((100,n_groups))
WR_percap_fastfood = np.empty((100,n_groups))
WR_percap = [WR_percap_cafe, WR_percap_restaurnt, WR_percap_fastfood]

for j in range(n_data):
    for i in range(1,n_groups):
        for k in range(100):
            WR_percap[j][k,i] = 1 - math.exp(-(p*q*times_erl[j][k,i]/R))
#############################################################
            
###### Erlang Per capita #############################
erl_percap = np.empty((n_data, n_groups))
for j in range(n_data):
    for i in range(1,n_groups):# we avoid i=1 since this doesnt result in distribtuion since only party size of 1. start at Party size = 2
        erl_percap[j,i] = 1- ((lam_erl[j,i]/((p*q/R)+lam_erl[j,i]))**k_erl[j,i])
##########################################################


#############Data frames for violin plot ####################

### Wells Riley Per cApita Data frame
dWRPC_cafe={'Party Size = 2' : WR_percap[0][:,1],'Party Size = 3' : WR_percap[0][:,2],'Party Size = 4' : WR_percap[0][:,3],'Party Size = 5' : WR_percap[0][:,4]}
dWRPC_restaurant={'Party Size = 2' : WR_percap[1][:,1],'Party Size = 3' : WR_percap[1][:,2],'Party Size = 4' : WR_percap[1][:,3],'Party Size = 5' : WR_percap[1][:,4]}
dWRPC_fastfood={'Party Size = 2' : WR_percap[2][:,1],'Party Size = 3' : WR_percap[2][:,2],'Party Size = 4' : WR_percap[2][:,3],'Party Size = 5' : WR_percap[2][:,4]}

df_WRPCcafe = pandas.DataFrame(dWRPC_cafe)
df_WRPCrestaurant= pandas.DataFrame(dWRPC_restaurant)
df_WRPCfastfood= pandas.DataFrame(dWRPC_fastfood)

df_WRPC = [df_WRPCcafe, df_WRPCrestaurant, df_WRPCfastfood]

##############Violin plot###############
colour = {'Party Size = 2': 'tab:blue', 'Party Size = 3':'tab:orange', 'Party Size = 4':'tab:purple', 'Party Size = 5':'tab:red'}#,care_type[4]: 'tab:olive'}

for i in range(n_data):
    plt.figure(dpi=750)#set dots per inch for better quality images
    sns.violinplot(data=df_WRPC[i], inner='quartiles', palette=colour)
    plt.ylim(0,1)
    plt.yticks(np.arange(0,1.1,0.1))
    plt.ylabel('Probability P(T) from Equation (1)')
    plt.scatter([0,1,2,3],erl_percap[i,1:], marker='o', color='k', label='Erlang Per Capita P(T)')
    plt.legend()
    plt.show()

