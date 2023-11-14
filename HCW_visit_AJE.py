# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:59:06 2023
This scirpt takes the data from a paper 
Figure 2 in https://doi.org/10.1016/j.jhin.2020.11.016
For the care episode duration for actual and mock activities (for a variety of tasks
                                                              - importnat point is that this is contact time with the patient).
 

The data is computed and then an erlang distribution is formed from this data. 


CREATED 28/06/2023 AJE
"""
import matplotlib.pyplot as plt
import random, math
import statistics as stat
import numpy as np
from scipy.integrate import odeint #for solving odes
import matplotlib.colors as mcolors #colour name package
from matplotlib.pyplot import cm #colour map package
from pywaffle import Waffle #visual package for visuallising icons
from matplotlib.ticker import StrMethodFormatter #to set decimal places in axis
import time #for live run time of code
from scipy.stats import binom # for the binomial distirbtion sampling 
from scipy.stats import erlang #for the erland distrbution
import sympy as sym
from round_up import round_up, rounding
from read_csv_WRdata import ReadCSV
from fitter import Fitter, get_common_distributions, get_distributions
import seaborn as sns
import pandas
##############################################################
######################## DATA IMPORT ########################
#filepath to a csv file of the data
filepath = r"C:\Users\scaje\OneDrive - University of Leeds\UNIV. OF LEEDS\PhD PROJECT\Ward Transmission\Stochastic_models\For Wells-Riley Paper\Data_for_results\HCW_visit-time_data\time_in_hospital_room.csv"

#importing the data for use in script
data = ReadCSV(filepath)

ID = data[1][:,0] # this has the ID of the activity
Type = data[1][:,4] #this contains the type of activity being completed, mock or actual
Time = data[1][:,6] #this has the duration of the activity in seconds

#split into the different types of activity, actual vs mock 
#Act = actual
#mock = mock
#idx101 has the last entry for actual care
#actual care
ActID = ID[:102] # ID data for actual care
ActType = Type[:102] #type data for actual care - should all be actual
ActTime = Time[:102] #time data for all actual care
#mock care
MockID = ID[102:] # ID data for mock  care - should all be actual
MockType = Type[102:] #type data for mock care - should all be mock
MockTime = Time[102:] #time data for all amock care
########################################################################
#########################################################################

########################################################################
################## PROBABILITY AFTER INFECTOR LEAVES ###################
########################################################################
#HERE WE HAVE TWO DIFFERENT PROBABILITies. 
#we explore the probability of infection when considering not just the period
#where the infector is present [0,T] say, but also the probability once the
#infector has left [T,t] say. In many hospital scenarios, the susceptible will
#remain indefinitely after the drug round and so we also consider the t->infinity
#All of the probability is considered here as a per capita rate

############
#Initialise parameters for infection outbreak scenario...
###########################Ventilation and room parameters ##################
#quanta rate =  quanta per min)
q=0.5
#pulmonary rate as volume/min)
p=0.01
#Ventilaation rate = m^3/min
Q=1.42 #for 3ach=1.42 #for1.5ach = 0.71 for 0.5ach = 0.23 #for 6ACH=2.84
#volume of indoor space v= m^3
v=28.57 # volume of a single bed room at St james J12 
#Time that the infector is present in Mins
T=14
###Initial number of infected individuals
I0=1
#Calculate the probability for each different visit time to see how this varies for the individual in the room
#calculate probability over initial period
P_inf = 1 - math.exp((p*q*I0/(Q**2))*(v*(1-math.exp(-Q*T/v))-Q*T))
print('Probability of infection whilst infector is present P_inf = %s' %(P_inf))

#Calculate the probability of infection over period when infector has left for long t.....
P_noinf =(math.exp((p*q*I0/(Q**2))*(v*(1-math.exp(-Q*T/v))-Q*T)))*(p*q*I0*v/(Q**2))*(1-math.exp(-Q*T/v))
print('Probability of infection when infector leaves (t->\inf) P_noinf = %s' %(P_noinf))

#Calculate the probability of infection over total period for long t.....
P_total = 1 - (math.exp((p*q*I0/(Q**2))*(v*(1-math.exp(-Q*T/v))-Q*T)))*(1-(p*q*I0*v/(Q**2))*(1-math.exp(-Q*T/v)))
print('Probability of infection in total P_total = %s' %(P_total))


###################################################
###Probability for all Actual Care Times
###################################################
#########Take the different times from the Actual care
ActTime = ActTime.astype(np.float) #convert times to floats to be used in python
n=len(ActTime)
#ActTime_mins = [ActTime[i]/60 for i in range(n)] #convert the Times from seconds into minutes for consistency with other parameters



#Time that the infector is present in Mins for each time
T=ActTime

#initialse probability vectors
P_inf = np.empty(n)
P_noinf = np.empty(n)
P_total = np.empty(n)


#Calculate the probability for each different visit time to see how this varies for the individual in the room
for i in range(n):
    #calculate probability over initial period
    P_inf[i] = 1 - math.exp((p*q*I0/(Q**2))*(v*(1-math.exp(-Q*T[i]/v))-Q*T[i]))

    #Calculate the probability of infection over period when infector has left for long t.....
    P_noinf[i] =(math.exp((p*q*I0/(Q**2))*(v*(1-math.exp(-Q*T[i]/v))-Q*T[i])))*(p*q*I0*v/(Q**2))*(1-math.exp(-Q*T[i]/v))

    #Calculate the probability of infection over total period for long t.....
    P_total[i] = 1 - (math.exp((p*q*I0/(Q**2))*(v*(1-math.exp(-Q*T[i]/v))-Q*T[i])))*(1-(p*q*I0*v/(Q**2))*(1-math.exp(-Q*T[i]/v)))


#print('Probability of infection whilst infector is present P_inf = %s' %(P_inf))
#print('Probability of infection when infector leaves (t->\inf) P_noinf = %s' %(P_noinf))
#print('Probability of infection in total P_total = %s' %(P_total))  

plt.figure(dpi=750)#set dots per inch for better quality images
plt.scatter(T, P_inf, label='Infector Present')
plt.scatter(T,P_noinf, label='After Infector Leaves')
plt.scatter(T, P_total, label='Total')
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Probability of Infection')
plt.legend(title='Time Period')
plt.title('probabilities for 3ACH')
plt.show()



#######################################################################
###Probability for all Actual Care Times and varying ventialtion rates
########################################################################
#multiple ventialtion rates
#we keep the same parameters other than ventialtion which changes 
 #for 3ach=1.42 #for1.5ach = 0.71 for 0.5ach = 0.23 #for 6ACH=2.84 # room of volume 28.57m^3 as in St James J12 single room 
Q_vary = [0.238, 0.714, 1.428, 2.856]

#initialse probability vectors
P_inf_Qvary = np.empty((n,len(Q_vary)))
P_noinf_Qvary = np.empty((n,len(Q_vary)))
P_total_Qvary = np.empty((n,len(Q_vary)))


plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(Q_vary)):# calculate for each ventialtion time
    
    #Calculate the probability for each different visit time to see how this varies for the individual in the room
    for i in range(n):
        #calculate probability over initial period
        P_inf_Qvary[i,j] = 1 - math.exp((p*q*I0/(Q_vary[j]**2))*(v*(1-math.exp(-Q_vary[j]*T[i]/v))-Q_vary[j]*T[i]))
    
        #Calculate the probability of infection over period when infector has left for long t.....
        P_noinf_Qvary[i,j] =(math.exp((p*q*I0/(Q_vary[j]**2))*(v*(1-math.exp(-Q_vary[j]*T[i]/v))-Q_vary[j]*T[i])))*(p*q*I0*v/(Q_vary[j]**2))*(1-math.exp(-Q_vary[j]*T[i]/v))
    
        #Calculate the probability of infection over total period for long t.....
        P_total_Qvary[i,j] = 1 - (math.exp((p*q*I0/(Q_vary[j]**2))*(v*(1-math.exp(-Q_vary[j]*T[i]/v))-Q_vary[j]*T[i])))*(1-(p*q*I0*v/(Q_vary[j]**2))*(1-math.exp(-Q_vary[j]*T[i]/v)))
    
    #plotting for the  total probability of infection
    plt.scatter(T, P_total_Qvary[:,j], label='Q=%s ACH' %(round(Q_vary[j]*60/v,2)))
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Probability of Infection')
plt.legend(title='Ventilation Rate')
plt.title('Total probability over whole period')
plt.show()

#plotting for the probabiluty when the infector is present
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(Q_vary)):
        plt.scatter(T, P_inf_Qvary[:,j], label='Q=%s ACH' %(round(Q_vary[j]*60/v,2)))
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Probability of Infection')
plt.legend(title='Ventilation Rate')
plt.title('probability whilst infector is present')
plt.show()

#plotting for when the infector leaves
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(Q_vary)):
        plt.scatter(T, P_noinf_Qvary[:,j], label='Q=%s ACH' %(round(Q_vary[j]*60/v,2)))
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Probability of Infection')
plt.legend(title='Ventilation Rate')
plt.title('probability after infector leaves')
plt.show()

########Think about the relative difference at each time point between the 
#additional risk acheieved once the infector leaves (hopefully show benefit of considering this)
P_reldif_Qvary = np.empty((n,len(Q_vary)))#initialise relative difference vectore
for i in range(n):
    for j in range(len(Q_vary)):
        if ActTime[i] == 0:#adressing the activity where infector enters and leaves in <1min (this gives zero event time in data)
            P_reldif_Qvary[i,j] = 0#set to zero to avoid NaN
        else:
            P_reldif_Qvary[i,j] = (P_noinf_Qvary[i,j]/P_total_Qvary[i,j])*100 #caulate the relative difference as a percentage

#Plot the relative difference as a scatter for each different vent rate
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(Q_vary)):
    plt.scatter(T, P_reldif_Qvary[:,j], label='Q=%s ACH' %(round(Q_vary[j]*60/v,2)))
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Relative difference [%]')
plt.legend(title='Ventilation Rate')
plt.title('Rel Dif of including P_noinf in risk assessment')
plt.show()
#######################################################################
#######################################################################
#######################################################################







#######################################################################
#######################################################################
#######################################################################
#######################################################################
###Probability for all Actual Care Times and varying quanta rates
########################################################################
#multiple quanta rates
#we keep the same parameters other than quanta which changes 
#different quanta values incl. 1q/hr=0.016667, 5q/hr=0.083334, 10q/hr=0.16667, 30q/hr=0.5
q_vary = [0.016667, 0.083334, 0.16667, 0.5] #these are in q/min

#initialse probability vectors
P_inf_qvary = np.empty((n,len(q_vary)))
P_noinf_qvary = np.empty((n,len(q_vary)))
P_total_qvary = np.empty((n,len(q_vary)))


plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(q_vary)):# calculate for each quanta rate
    
    #Calculate the probability for each different visit time to see how this varies for the individual in the room
    for i in range(n):
        #calculate probability over initial period
        P_inf_qvary[i,j] = 1 - math.exp((p*q_vary[j]*I0/(Q**2))*(v*(1-math.exp(-Q*T[i]/v))-Q*T[i]))
    
        #Calculate the probability of infection over period when infector has left for long t.....
        P_noinf_qvary[i,j] =(math.exp((p*q_vary[j]*I0/(Q**2))*(v*(1-math.exp(-Q*T[i]/v))-Q*T[i])))*(p*q_vary[j]*I0/(Q**2))*(1-math.exp(-Q*T[i]/v))
    
        #Calculate the probability of infection over total period for long t.....
        P_total_qvary[i,j] = 1 - (math.exp((p*q_vary[j]*I0/(Q**2))*(v*(1-math.exp(-Q*T[i]/v))-Q*T[i])))*(1-(p*q_vary[j]*I0/(Q**2))*(1-math.exp(-Q*T[i]/v)))
    
    #plotting for the  total probability of infection
    plt.scatter(T, P_total_qvary[:,j], label='q=%s q/hr' %(int(q_vary[j]*60)))
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Probability of Infection')
plt.legend(title='Quanta Production Rate')
plt.title('Total probability over whole period')
plt.show()

#plotting for the probabiluty when the infector is present
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(q_vary)):
        plt.scatter(T, P_inf_qvary[:,j], label='q=%s q/hr' %(int(q_vary[j]*60)))
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Probability of Infection')
plt.legend(title='Quanta Production Rate')
plt.title('probability whilst infector is present')
plt.show()

#plotting for when the infector leaves
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(q_vary)):
        plt.scatter(T, P_noinf_qvary[:,j], label='q=%s q/hr' %(int(q_vary[j]*60)))
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Probability of Infection')
plt.legend(title='Quanta Production Rate')
plt.title('probability after infector leaves')
plt.show()

########Think about the relative difference at each time point between the 
#additional risk acheieved once the infector leaves (hopefully show benefit of considering this)
P_reldif_qvary = np.empty((n,len(q_vary)))#initialise relative difference vectore
for i in range(n):
    for j in range(len(q_vary)):
        if ActTime[i] == 0:#adressing the activity where infector enters and leaves in <1min (this gives zero event time in data)
            P_reldif_qvary[i,j] = 0#set to zero to avoid NaN
        else:
            P_reldif_qvary[i,j] = (P_noinf_qvary[i,j]/P_total_qvary[i,j])*100 #caulate the relative difference as a percentage

#Plot the relative difference as a scatter for each different vent rate
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(q_vary)):
    plt.scatter(T, P_reldif_qvary[:,j], label='q=%s q/hr' %(int(q_vary[j]*60)))
plt.xlabel('Time Infector is Present [mins]')
plt.ylabel('Relative difference [%]')
plt.legend(title='Quanta Production Rate')
plt.title('Rel Dif of including P_noinf in risk assessment')
plt.show()
#######################################################################
#######################################################################
#######################################################################






#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
###Probability for all Actual Care Times with altered REMOVAL RATE
########################################################################

######################################################
#Initialise parameters for infection outbreak scenario...
###########################Ventilation and room parameters ##################
#quanta rate =  quanta per min)
q=6
#pulmonary rate as volume/min)
p=0.01
#volume of indoor space v= m^3
v=28.57 # volume of a single bed room at St james J12 
###Initial number of infected individuals
I0=1
###########################################


############################################################
##### NEW PARAMETERS FOR REMOVAL RATE 
#ranges define first and last, then values in bewteen at equal steps
r_i = [0,0.0035,0.007,0.0105] # this is viral decay/inactivation rate - range is given as [0,0.63] per hour in paper 
r_d = [0.005,0.0115,0.01815,0.025] #This is viral deposition rate - range is given as [0.3,1.5] per hour

#multiple ventilation rates 
 #for 3ach=1.42 #for1.5ach = 0.71 for 0.5ach = 0.23 #for 6ACH=2.84 # room of volume 28.57m^3 as in St James J12 single room 
Q_vary = [0.238, 0.714, 1.428, 2.856]

##The REMOVAL RATE
R = [Q_vary[i] + (v*r_d[i]) + (v*r_i[i]) for i in range(len(Q_vary))]
################################################

colour = ['tab:blue', 'tab:orange', 'tab:purple', 'tab:red']

#########
#initialse probability vectors
P_inf_Rvary = np.empty((n,len(R)))
P_noinf_Rvary = np.empty((n,len(R)))
P_total_Rvary = np.empty((n,len(R)))


plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(R)):# calculate for each ventialtion time
    
    #Calculate the probability for each different visit time to see how this varies for the individual in the room
    for i in range(n):
        #calculate probability over initial period
        P_inf_Rvary[i,j] = 1 - math.exp((p*q*I0/(R[j]**2))*(v*(1-math.exp(-R[j]*T[i]/v))-R[j]*T[i]))
    
        #Calculate the probability of infection over period when infector has left for long t.....
        P_noinf_Rvary[i,j] =(math.exp((p*q*I0/(R[j]**2))*(v*(1-math.exp(-R[j]*T[i]/v))-R[j]*T[i])))*(p*q*I0/(R[j]**2))*(1-math.exp(-R[j]*T[i]/v))
    
        #Calculate the probability of infection over total period for long t.....
        P_total_Rvary[i,j] = 1 - (math.exp((p*q*I0/(R[j]**2))*(v*(1-math.exp(-R[j]*T[i]/v))-R[j]*T[i])))*(1-(p*q*I0/(R[j]**2))*(1-math.exp(-R[j]*T[i]/v)))
    
    #plotting for the  total probability of infection
    plt.scatter(T, P_total_Rvary[:,j], label='R=%s ' %(round(R[j],3)), color=colour[j])
plt.xlabel('Time Infector is Present [min]')
plt.ylabel('Probability of Infection')
plt.legend(title='Removal Rate [$m^3 min^{-1}$]')
#plt.title('Total probability over whole period')
plt.xticks(np.arange(0,22,2))
plt.yscale('log')
plt.ylim(7e-06, 10e-01)
plt.show()

#plotting for the probabiluty when the infector is present
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(R)):
        plt.scatter(T, P_inf_Rvary[:,j], label='R=%s' %(round(R[j],3)), color=colour[j])
plt.xlabel('Time Infector is Present [min]')
plt.ylabel('Probability of Infection')
plt.legend(title='Removal Rate [$m^3 min^{-1}$]')
#plt.title('probability whilst infector is present')
plt.xticks(np.arange(0,22,2))
plt.yscale('log')
plt.ylim(7e-06, 10e-01)
plt.show()

#plotting for when the infector leaves
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(R)):
        plt.scatter(T, P_noinf_Rvary[:,j], label='R=%s ' %(round(R[j],3)), color=colour[j])
plt.xlabel('Time Infector is Present [min]')
plt.ylabel('Probability of Infection')
plt.legend(title='Removal Rate [$m^3 min^{-1}$]')
#plt.title('probability after infector leaves')
plt.xticks(np.arange(0,22,2))
plt.yscale('log')
plt.ylim(7e-06, 10e-01)
plt.show()

########Think about the relative difference at each time point between the 
#additional risk acheieved once the infector leaves (hopefully show benefit of considering this)
P_reldif_Rvary = np.empty((n,len(R)))#initialise relative difference vectore
for i in range(n):
    for j in range(len(R)):
        if ActTime[i] == 0:#adressing the activity where infector enters and leaves in <1min (this gives zero event time in data)
            P_reldif_Rvary[i,j] = 0#set to zero to avoid NaN
        else:
            P_reldif_Rvary[i,j] = (P_noinf_Rvary[i,j]/P_total_Rvary[i,j])*100 #caulate the relative difference as a percentage

#Plot the relative difference as a scatter for each different vent rate
plt.figure(dpi=750)#set dots per inch for better quality images
for j in range(len(R)):
    plt.scatter(T, P_reldif_Rvary[:,j],label='R=%s ' %(round(R[j],3)), color=colour[j])
plt.xlabel('Time Infector is Present [min]')
plt.ylabel('Contribution [%]')
plt.legend(fontsize=7,  title='Removal Rate [$m^3 min^{-1}$]')
#plt.title('Rel Dif of including P_noinf in risk assessment')
plt.xticks(np.arange(0,22,2))
plt.yscale('log')
#plt.ylim(10e-06, 10e02)
plt.show()
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################




mu = stat.mean(ActTime)
sd = math.sqrt(stat.variance(ActTime))






#example of colourbar use for scatter from resulst_vs_wth.py script
# =============================================================================
# for i in range(int(np.shape(exposures)[1])):
#     plt.figure(dpi=750)#set dots per inch for better quality images
#     plt.scatter(wind_dir, wind_speed, c=conc_hr[:,i], cmap='Reds')
#     plt.xlabel("Wind Direction [deg]") 
#     plt.ylabel("Wind Speed [$ms^{-1}$]")
#     #plt.legend(loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Zones')
#     plt.colorbar(label = 'Concentration [$qm^{-3}$]')
#     plt.title('Comparing weather conditions in Zone %s' %(i+1))
#     #plt.legend()
#     plt.show()
# =============================================================================

#######################################################################
#######################################################################
#######################################################################
################## VIOLIN PLOTS FOR EACH CARE TYPE ####################
#######################################################################
#######################################################################


#NEED TO SPLIT DATA INTO EACH CARE TYPE (USING ACTUAL CARE ONLY - NOT MOCK)
#RECALL ACTUAL CARE VALUES
#ActID = ID[:102] # ID data for actual care
#ActType = Type[:102] #type data for actual care - should all be actual
#ActTime = Time[:102] #time data for all actual care
#RECALL DATA[1][:,1] GIVES THE CARE TYPE


#CARE TYPES
#set() finds unique values in vector
#list() converts to a list
care_type = list(set(data[1][:102,1]))

#for all care types, split times into each unique care type for actual care only
#### HCW VISIT TIMES 
bloods_time = ActTime[:13]
IV_time = ActTime[13:33]
checks_time = ActTime[33:46]
obs_time =ActTime[46:78]
DRrounds_time = ActTime[78:102]

#INFECTION PROBABILITY
#TOTAL PROBABILITY TO SUSCEPTIBLE IN ROOM 
bloods_totalrisk = P_total_Rvary[:13,:]
IV_totalrisk = P_total_Rvary[13:33,:]
checks_totalrisk = P_total_Rvary[33:46,:]
obs_totalrisk = P_total_Rvary[46:78,:]
DRrounds_totalrisk = P_total_Rvary[78:102,:]





#MAKING A DATAFRAME FOR THE VISIT TIMES AND CARE TYPE
dtime={'Care Type': data[1][:102,1], 'Visit Length': ActTime[:]}
df_time = pandas.DataFrame(dtime)

#MAKING A DF FOR THE INFECTION RISK FOR DIFFERENT REMOVAL RATES, AND CARE TYPES
drisk = {'Care Type': data[1][:102,1], 'R = %s' %(round(R[0],3)): P_total_Rvary[:,0] , 'R = %s' %(round(R[1],3)): P_total_Rvary[:,1] ,'R = %s' %(round(R[2],3)): P_total_Rvary[:,2] ,'R = %s' %(round(R[3],3)): P_total_Rvary[:,3] }
df_risk = pandas.DataFrame(drisk)

######### The plots

colour = {care_type[0]: 'tab:blue', care_type[1]:'tab:orange', care_type[2]:'tab:purple', care_type[3]:'tab:red',care_type[4]: 'tab:olive'}
#my_pal = {"versicolor": "g", "setosa": "b", "virginica": "m"}

#Violin plot for the times
plt.figure(dpi=750)#set dots per inch for better quality images
sns.violinplot(x='Care Type', y='Visit Length', data=df_time, inner="quart", linewidth=1, palette=colour)
plt.xlabel('Care Types')
plt.ylabel('Duration [min]')
plt.show()


#Violin plot for the the different removal rates across each care type
for i in range(len(R)):
    plt.figure(dpi=750)#set dots per inch for better quality images
    sns.violinplot(x='Care Type', y='R = %s' %(round(R[i],3)), data=df_risk, inner="quart", linewidth=1, palette=colour)
    plt.xlabel('Care Types')
    plt.ylabel('Probability')
    #plt.title('Removal rate R = %s'%(round(R[i],3)))
    #plt.yscale('log')
    #plt.ylim(0,)
    plt.ylim(-0.085,0.475)
    #plt.yscale('log')
    plt.show()
#current seaborn version 0.11.0





#Violin plot for the the different removal rates across each care type
#plotting all four removal rates across different care types

#data frames
dbloods = {'Care Type': data[1][:13,1], 'R = %s' %(round(R[0],3)): P_total_Rvary[:13,0] , 'R = %s' %(round(R[1],3)): P_total_Rvary[:13,1] ,'R = %s' %(round(R[2],3)): P_total_Rvary[:13,2] ,'R = %s' %(round(R[3],3)): P_total_Rvary[:13,3] }
dIV = {'Care Type': data[1][13:33,1], 'R = %s' %(round(R[0],3)): P_total_Rvary[13:33,0] , 'R = %s' %(round(R[1],3)): P_total_Rvary[13:33,1] ,'R = %s' %(round(R[2],3)): P_total_Rvary[13:33,2] ,'R = %s' %(round(R[3],3)): P_total_Rvary[13:33,3] }
dchecks = {'Care Type': data[1][33:46,1], 'R = %s' %(round(R[0],3)): P_total_Rvary[33:46,0] , 'R = %s' %(round(R[1],3)): P_total_Rvary[33:46,1] ,'R = %s' %(round(R[2],3)): P_total_Rvary[33:46,2] ,'R = %s' %(round(R[3],3)): P_total_Rvary[33:46,3] }
dobs = {'Care Type': data[1][46:78,1], 'R = %s' %(round(R[0],3)): P_total_Rvary[46:78,0] , 'R = %s' %(round(R[1],3)): P_total_Rvary[46:78,1] ,'R = %s' %(round(R[2],3)): P_total_Rvary[46:78,2] ,'R = %s' %(round(R[3],3)): P_total_Rvary[46:78,3] }
ddr = {'Care Type': data[1][78:102,1], 'R = %s' %(round(R[0],3)): P_total_Rvary[78:102,0] , 'R = %s' %(round(R[1],3)): P_total_Rvary[78:102,1] ,'R = %s' %(round(R[2],3)): P_total_Rvary[78:102,2] ,'R = %s' %(round(R[3],3)): P_total_Rvary[78:102,3] }

df_bloods = pandas.DataFrame(dbloods)
df_IV = pandas.DataFrame(dIV)
df_checks = pandas.DataFrame(dchecks)
df_obs = pandas.DataFrame(dobs)
df_dr = pandas.DataFrame(ddr)

df = [df_dr, df_IV,  df_obs, df_checks,df_bloods]

for i in range(len(df)):
    plt.figure(dpi=750)#set dots per inch for better quality images
    sns.violinplot(data=df[i], inner="quart", linewidth=1)#, palette=colour)
    #plt.xlabel('Care Types')
    #plt.ylabel('Probability')
    plt.title('Care type = %s'%(care_type[i]))
    #plt.yscale('log')
    #plt.ylim(0,)
    #plt.ylim(-0.085,0.475)
    #plt.yscale('log')
    plt.show()


####trying to log data
# =============================================================================
# from matplotlib import ticker as mticker
# 
# log_data = np.empty((len(P_total_Rvary),len(R)))
# for i in range(len(P_total_Rvary)):
#     for j in range(len(R)):
#         log_data[i,j] = np.log10(P_total_Rvary[i,j])
#         
# #log_data = [[[np.log10(d) for d in row] for row in P_total_Rvary]
# #MAKING A DF FOR THE INFECTION RISK FOR DIFFERENT REMOVAL RATES, AND CARE TYPES
# drisk_log = {'Care Type': data[1][:102,1], 'R = %s' %(round(R[0],3)): log_data[:,0] , 'R = %s' %(round(R[1],3)): log_data[:,1] ,'R = %s' %(round(R[2],3)): log_data[:,2] ,'R = %s' %(round(R[3],3)): log_data[:,3] }
# df_risk_log = pandas.DataFrame(drisk)
# 
# for i in range(len(R)):
#     fig, ax = plt.subplots(dpi=750)
#     sns.violinplot(x='Care Type', y='R = %s' %(round(R[i],3)), data=df_risk, inner="quart", linewidth=1, palette=colour, ax=ax)
#     ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
#     ymin, ymax = ax.get_ylim()
#     tick_range = np.arange(np.floor(ymin), ymax)
#     ax.yaxis.set_ticks(tick_range)
#     ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)
#     plt.tight_layout()
#     plt.show()
# =============================================================================



#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################