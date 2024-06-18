# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:59:06 2023

This script outlines the code used to produce the figures in Case study 1 in the manuscript
 "The Wells-Riley Model Revisited: Randomness, Heterogeneity and Transient Behaviours"; 
 Alexander J. Edwards, Marco-Felipe King, Catherine J. Noakes, Daniel Peckham, Martin Lopez-Garcia.

CREATED 28/06/2023 AJE
"""
import matplotlib.pyplot as plt
import math
import statistics as stat
import numpy as np
from read_csv_WRdata import ReadCSV
import seaborn as sns
import pandas
##############################################################
######################## DATA IMPORT ########################
#filepath to a csv file of the data
filepath = r"time_in_hospital_room.csv"

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

###################################################
###Probability for all Actual Care Times
###################################################
#########Take the different times from the Actual care
ActTime = ActTime.astype(np.float) #convert times to floats to be used in python
n=len(ActTime)

#Time that the infector is present in Mins for each time
T=ActTime

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
        P_noinf_Rvary[i,j] =1 - math.exp(-(p*q*I0*v/(R[j]**2))*(1-math.exp(-R[j]*T[i]/v)))
    
        #Calculate the probability of infection over total period for long t.....
        P_total_Rvary[i,j] = 1 - math.exp(-((p*q*I0*T[i])/(R[j])))
        
        
        
    #plotting for the  total probability of infection
    plt.scatter(T, P_total_Rvary[:,j], label='R=%s ' %(round(R[j],3)), color=colour[j])
plt.xlabel('Time Infector is Present [min]')
plt.ylabel('Probability of Infection')
plt.legend(title='Removal Rate [$m^3 min^{-1}$]')
#plt.title('Total probability over whole period')
plt.xticks(np.arange(0,22,2))
plt.yscale('log')
plt.ylim(1e-05, 1.5e0)
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
plt.ylim(1e-05, 1.5e0)
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
plt.ylim(1e-05, 1.5e0)
plt.show()

########Think about the relative difference at each time point between the 
#additional risk acheieved once the infector leaves (hopefully show benefit of considering this)
P_reldif_Rvary = np.empty((n,len(R)))#initialise relative difference vectore
for i in range(n):
    for j in range(len(R)):
        if ActTime[i] == 0:#adressing the activity where infector enters and leaves in <1min (this gives zero event time in data)
            P_reldif_Rvary[i,j] = 0#set to zero to avoid NaN
        else:
            P_reldif_Rvary[i,j] = ((P_total_Rvary[i,j]-P_inf_Rvary[i,j])/P_total_Rvary[i,j])*100 #caulate the relative difference as a percentage

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


#######################################################################
#######################################################################
#######################################################################
################## VIOLIN PLOTS FOR EACH CARE TYPE ####################
#######################################################################
#######################################################################


#CARE TYPES
#set() finds unique values in vector
#list() converts to a list
#sorted() converts to a list but maintains the order of the input
care_type = sorted(set(data[1][:102,1]))

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

colour = {care_type[0]: 'tab:olive', care_type[1]:'tab:red', care_type[2]:'tab:orange', care_type[3]:'tab:purple',care_type[4]: 'tab:blue'}

#Violin plot for the times
plt.figure(dpi=750)#set dots per inch for better quality images
sns.violinplot(x='Care Type', y='Visit Length', data=df_time, inner="quart", linewidth=1, palette=colour, cut=0)
plt.xlabel('Care Types')
plt.ylabel('Duration [min]')
plt.yticks(np.arange(0,25,5))
plt.show()


#Violin plot for the the different removal rates across each care type
for i in range(len(R)):
    plt.figure(dpi=750)#set dots per inch for better quality images
    sns.violinplot(x='Care Type', y='R = %s' %(round(R[i],3)), data=df_risk, inner="quart", linewidth=1, palette=colour, cut=0)#, bw=0.3, cut=0)
    plt.xlabel('Care Types')
    plt.ylabel('Probability')
    #plt.title('Removal rate R = %s'%(round(R[i],3)))
    plt.ylim(-0.05,1)
    plt.show()
#current seaborn version 0.11.0


#######################################################################
#######################################################################