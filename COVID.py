#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import math
import os
import sys
from os.path import exists
import csv

def initial_infected(grid_size,virus_center,virus_size,virus_numb):
    n1 = grid_size[0]
    n2 = grid_size[1]
    x0 = virus_center[0]
    y0 = virus_center[1]
    s1 = virus_size[0]
    s2 = virus_size[1]
    
    a=[0,1]
    initial_virus_prod_rate = virus_numb/(s1*s2)
    M_infected = np.zeros(n1*n2).reshape(n1,n2)
    
    if s1==1 and s2==1:
        M_infected[x0][y0]=1
    elif s1==n1 and s2==n2:
        M_infected=np.random.choice(a,size=n1*n2,p=[1-initial_virus_prod_rate,initial_virus_prod_rate]).reshape(n1,n2)
    else:
        M = np.random.choice(a,size=s1*s2,p=[1-initial_virus_prod_rate,initial_virus_prod_rate]).reshape(s1,s2)
        M_infected[int(x0-s1/2):int(x0+s1/2),int(y0-s2/2):int(y0+s2/2)]=M
        M_infected=M_infected.astype(int)
    return M_infected

def initial_virus_prod(grid_size,virus_center,virus_prods_n,virus_prods_p):
    n1 = grid_size[0]
    n2 = grid_size[1]
    x0 = virus_center[0]
    y0 = virus_center[1]  
    a_virus_prod = np.round(np.random.negative_binomial(virus_prods_n,virus_prods_p,n1*n2))
    M_virus_prod = a_virus_prod.reshape(n1,n2)
    M_virus_prod = M_virus_prod.astype(int)

                                                    
    return M_virus_prod

def initial_ifn_prod(grid_size,a,ifn_prob,ifn_prods):
    n1 = grid_size[0]
    n2 = grid_size[1]
    M_ifn_prod = np.round(np.random.choice(a,size=n1*n2,p=[1-ifn_prob,ifn_prob]).reshape(n1,n2)*(float(ifn_prods)/ifn_prob))
    M_ifn_prod = M_ifn_prod.astype(int)
    return M_ifn_prod

def initial_lifespan(grid_size,lifespan_mu,lifespan_sigma):
    n1 = grid_size[0]
    n2 = grid_size[1]    
    a_lifespan = np.round(np.random.normal(lifespan_mu,lifespan_sigma,size=n1*n2))
    for i in range(n1*n2):
        if a_lifespan[i]<1:
            a_lifespan[i]=0
   
    M_lifespan = a_lifespan.reshape(n1,n2).astype(int)
    return M_lifespan

def plotmatrix(M):
    im = plt.imshow(M, cmap='Greys') 
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar(im, orientation='horizontal')
    plt.show()
    return 0

def plothistogram(h,name):
    plt.hist(h,bins=100) 
    plt.xlabel(name)
    plt.ylabel('counts')
    plt.show()
    return 0

def matrix_iter_stochastic(M,Par):
    M_virus_prod=M[0]
    M_ifn_prod=M[1]
    M_infected=M[2]
    M_protected=M[3]
    M_dead=M[4]
    M_exposure = M[5]
    M_exposure_timer=M[6]
    M_lifespan=M[7]
    M_virus_reduction=M[8]
    M_viruscell = M[9]
    M_virus_contact = M[10]
    M_ifncell = M[11]
    M_ifn_contact = M[12]
    
    prob_infect = Par[0]
    virus_prod_delay = Par[1] 
    ifn_prod_delay = Par[2]
    virus_diff = Par[3]
    ifn_diff = Par[4]
    virus_reduction_factor=Par[5]
    
    n1=len(M_infected)
    n2=len(M_infected[0])
    
    M_virus_prod1 = M_virus_prod.copy()
    M_infected1=np.zeros(n1*n2).reshape(n1,n2)
    M_protected1=np.zeros(n1*n2).reshape(n1,n2)    
    M_virus_reduction1 = np.zeros(n1*n2).reshape(n1,n2)
    M_viruscell1=np.zeros(n1*n2).reshape(n1,n2)
    M_virus_contact1=np.zeros(n1*n2).reshape(n1,n2)
    M_ifncell1=np.zeros(n1*n2).reshape(n1,n2)
    M_ifn_contact1=np.zeros(n1*n2).reshape(n1,n2)
    
    for i in range(n1):
        for j in range(n2):
            #If the cell is protected, its virus production reduces.
            if M_virus_reduction[i][j]>0:
                M_virus_prod1[i][j] = np.round(virus_reduction_factor*M_virus_prod[i][j]).astype(int)
            
            if (M_dead[i][j] > 0):
                M_infected[i][j]=0
            else:
                if (M_virus_prod1[i][j]>0 and M_exposure_timer[i][j] > virus_prod_delay):   
                    M_viruscell1[i][j]=1
                    dt = np.random.exponential(scale=1/prob_infect,size = M_virus_prod1[i][j])
                    r = np.abs(np.random.normal(loc=0,scale=np.sqrt(2*virus_diff*dt)))
        
                    theta = 2*math.pi*np.random.rand(M_virus_prod1[i][j])
                    x=np.round(r*np.sin(theta)).astype(int)
                    y=np.round(r*np.cos(theta)).astype(int)
                    #print(["virus,dt,r,theta:",dt,r,theta])
                    for idx, xi in np.ndenumerate(x):
                        a = i+x[idx]
                        b = j+y[idx] 
                        #print(["virus,idx,i,j,a,b:",idx,i,j,a,b])
                        if (x[idx]!=0 or y[idx]!=0) and a>=0 and a<n1 and b>=0 and b<n2:
                            M_virus_contact1[a,b]=1
                            if M_exposure_timer[a,b]<virus_prod_delay:
                                M_exposure[a,b] = 1 + M_exposure[a,b]
                            if M_infected[a,b]==0 and M_protected[a,b]==0 and M_dead[a,b]==0:
                                M_infected1[a,b]=1
                                                      
                if (M_ifn_prod[i][j] > 0 and M_exposure_timer[i][j] >  ifn_prod_delay):
                    M_ifncell1[i][j]=1
                    
                    dt = np.random.exponential(scale=1/prob_infect,size = M_ifn_prod[i][j])
                    r = np.abs(np.random.normal(loc=0,scale=np.sqrt(2*ifn_diff*dt)))
                    theta = 2*math.pi*np.random.rand(M_ifn_prod[i][j])
                    x=np.round(r*np.sin(theta)).astype(int)
                    y=np.round(r*np.cos(theta)).astype(int) 
                    #print(["ifn,dt,r,theta:",dt,r,theta])
                    for idx, xi in np.ndenumerate(x):
                        a = i+x[idx]
                        b = j+y[idx]  
                        #print(["ifn,idx,i,j,a,b:",idx,i,j,a,b])
                        if (x[idx]!=0 or y[idx]!=0) and a>=0 and a<n1 and b>=0 and b<n2:
                            M_ifn_contact1[a,b]=1

                            #if M_virus_reduction[a,b]==0:
                            M_virus_reduction1[a,b] = 1

                            if M_infected[a,b]==0 and M_protected[a,b]==0 and M_dead[a,b]==0:
                                M_protected1[a,b]=1
                            
                                 
    # the original infected add the new infected
    M_infected = M_infected + M_infected1
    M_infected = 1*(M_infected>0)
    # the original protected add the new protected
    M_protected = M_protected+M_protected1
    M_protected = 1*(M_protected>0)
    # the original virus reduction layer add the new virus reduction layer
    M_virus_reduction = M_virus_reduction+M_virus_reduction1
    M_virus_reduction = 1*(M_virus_reduction>0)
    # add the exposure timer 
    M_exposure_timer = M_exposure_timer + M_infected
    # if the exposure timer larger than lifespan, the cell is dead
    M_dead = M_exposure_timer>M_lifespan
    M_dead = 1*M_dead     

    Mtot = [M_virus_prod,M_ifn_prod,M_infected,M_protected,M_dead,M_exposure,M_exposure_timer,M_lifespan,M_virus_reduction,M_viruscell1,M_virus_contact1,M_ifncell1,M_ifn_contact1]
    return Mtot
    
def matrix_iter_homogenous(M,Par):
    M_virus_prod=M[0]
    M_ifn_prod=M[1]
    M_infected=M[2]
    M_protected=M[3]
    M_dead=M[4]
    M_exposure = M[5]
    M_exposure_timer=M[6]
    M_lifespan=M[7]
    M_virus_reduction=M[8]
    M_viruscell = M[9]
    M_virus_contact = M[10]
    M_ifncell = M[11]
    M_ifn_contact = M[12]
    
    prob_infect = Par[0]
    virus_prod_delay = Par[1] 
    ifn_prod_delay = Par[2]
    virus_diff = Par[3]
    ifn_diff = Par[4]
    virus_reduction_factor=Par[5]
    
    n1=len(M_infected)
    n2=len(M_infected[0])
    
    M_virus_prod1 = M_virus_prod.copy()
    M_infected1=np.zeros(n1*n2).reshape(n1,n2)
    M_protected1=np.zeros(n1*n2).reshape(n1,n2)    
    M_virus_reduction1 = np.zeros(n1*n2).reshape(n1,n2)
    M_viruscell1=np.zeros(n1*n2).reshape(n1,n2)
    M_virus_contact1=np.zeros(n1*n2).reshape(n1,n2)
    M_ifncell1=np.zeros(n1*n2).reshape(n1,n2)
    M_ifn_contact1=np.zeros(n1*n2).reshape(n1,n2)
    for i in range(n1):
        for j in range(n2):
            
            #If the cell is protected, its virus production reduces.
            if M_virus_reduction[i][j]>0:
                M_virus_prod1[i][j] = np.round(virus_reduction_factor*M_virus_prod[i][j]).astype(int)
            
            
            if (M_dead[i][j] > 0):
                M_infected[i][j]=0
            else:
                if (M_virus_prod1[i][j]>0 and M_exposure_timer[i][j] > virus_prod_delay):   
                    M_viruscell1[i][j]=1
                    infprob=np.random.choice(a=[0,1],size=M_virus_prod1[i][j],p=[1-prob_infect, prob_infect])
                    
                    virusx = np.round(np.random.rand(M_virus_prod1[i][j])*n1)
                    virusy = np.round(np.random.rand(M_virus_prod1[i][j])*n2)
                    
                    for idx, xi in np.ndenumerate(virusx):
                        a = np.round(virusx[idx]).astype(int)
                        b = np.round(virusy[idx]).astype(int)
                        infprobi = infprob[idx]
                        if infprobi!=0:
                            if (a!=i or b!=j) and a>=0 and a<n1 and b>=0 and b<n2:
                                M_virus_contact1[a,b]=1
                                if M_exposure_timer[a,b]<virus_prod_delay:
                                    M_exposure[a,b] = 1 + M_exposure[a,b]
                                if M_infected[a,b]==0 and M_protected[a,b]==0 and M_dead[a,b]==0:
                                    M_infected1[a,b]=1
                            
                if (M_ifn_prod[i][j] > 0 and M_exposure_timer[i][j] >  ifn_prod_delay):
                    M_ifncell1[i][j]=1
                    protprob=np.random.choice(a=[0,1],size=M_ifn_prod[i][j],p=[1-prob_infect, prob_infect])
                    ifnx = np.round(np.random.rand(M_ifn_prod[i][j])*n1)
                    ifny = np.round(np.random.rand(M_ifn_prod[i][j])*n2)

                            
                    for idx, xi in np.ndenumerate(ifnx):
                        a = np.round(ifnx[idx]).astype(int)
                        b = np.round(ifny[idx]).astype(int)
                        protprobi =protprob[idx]
                        if protprobi!=0:
                            if (a!=i or b!=j) and a>=0 and a<n1 and b>=0 and b<n2:

                                M_ifn_contact1[a,b]=1
                                if M_virus_reduction[a,b]==0:
                                    M_virus_reduction1[a,b] = 1

                                if M_infected[a,b]==0 and M_protected[a,b]==0 and M_dead[a,b]==0:
                                    M_protected1[a,b]=1
                                 
    # the original infected add the new infected
    M_infected = M_infected + M_infected1
    M_infected = 1*(M_infected>0)
    # the original protected add the new protected
    M_protected = M_protected+M_protected1
    M_protected = 1*(M_protected>0)
    # the original virus reduction layer add the new virus reduction layer
    M_virus_reduction = M_virus_reduction+M_virus_reduction1
    M_virus_reduction = 1*(M_virus_reduction>0)
    # add the exposure timer 
    M_exposure_timer = M_exposure_timer + M_infected
    # if the exposure timer larger than lifespan, the cell is dead
    M_dead = M_exposure_timer>M_lifespan
    M_dead = 1*M_dead     
    Mtot = [M_virus_prod,M_ifn_prod,M_infected,M_protected,M_dead,M_exposure,M_exposure_timer,M_lifespan,M_virus_reduction,M_viruscell1,M_virus_contact1,M_ifncell1,M_ifn_contact1]
    
    return Mtot

def Run(irun, step, Par0, IsImage=0,IsStochastic=1):
    grid_size=Par0[0]
    virus_size=Par0[1]
    virus_numb=Par0[2]
    virus_center=Par0[3]
    a=Par0[4]
    virus_prods=Par0[5]
    virus_prods_n=Par0[6]
    virus_prods_p=Par0[7]
    virus_diff=Par0[8]
    ifn_prods=Par0[9]
    ifn_diff=Par0[10]
    ifn_prob=Par0[11]
    virus_reduction_factor=Par0[12]
    prob_infect=Par0[13]
    virus_prod_delay=Par0[14]
    ifn_prod_delay=Par0[15]
    lifespan_mu=Par0[16]
    lifespan_sigma=Par0[17]
    n1=int(grid_size[0])
    n2=int(grid_size[1])
    s1=int(virus_size[0])
    s2=int(virus_size[1])
    filename = "virusprods{}_virusdiff{}_ifnprods{}_ifndiff{}_ifnprob{}_virusreduct{}_isstochastic{}_initialvirus{}_{}_{}_run{}".format(virus_prods,virus_diff,ifn_prods,ifn_diff,ifn_prob,virus_reduction_factor,IsStochastic,s1,s2,virus_numb,irun)
    Par1 = [prob_infect,virus_prod_delay, ifn_prod_delay, virus_diff,ifn_diff,virus_reduction_factor]
  
    M_infected=initial_infected(grid_size,virus_center,virus_size,virus_numb)
    M_protected=np.zeros(n1*n2).reshape(n1,n2)
    M_dead=np.zeros(n1*n2).reshape(n1,n2)
        
    M_virus_prod=initial_virus_prod(grid_size,virus_center,virus_prods_n,virus_prods_p)
    if(s1==1 and s2==1):
        while M_virus_prod[virus_center[0]][virus_center[1]] == 0 :
            M_virus_prod=initial_virus_prod(grid_size,virus_center,virus_prods_n,virus_prods_p)

    M_ifn_prod=initial_ifn_prod(grid_size,a,ifn_prob,ifn_prods)
    M_lifespan=initial_lifespan(grid_size,lifespan_mu,lifespan_sigma)  
    
    M_exposure=M_infected*2
    M_exposure_timer=M_infected*(virus_prod_delay+1)
    M_virus_reduction = np.zeros(n1*n2).reshape(n1,n2)

    M_ifncell = np.zeros(n1*n2).reshape(n1,n2)       # cell generating ifn
    M_ifn_contact = np.zeros(n1*n2).reshape(n1,n2)   # cell contacting ifn
    M_viruscell = np.zeros(n1*n2).reshape(n1,n2)     # cell generating virus
    M_virus_contact = np.zeros(n1*n2).reshape(n1,n2) # cell contacting virus
 
    Mtot = [M_virus_prod,M_ifn_prod,M_infected,M_protected,M_dead,M_exposure,M_exposure_timer,M_lifespan,M_virus_reduction,M_viruscell,M_virus_contact,M_ifncell,M_ifn_contact]
    if exists("./data_py/data_{}.csv".format(filename)):
        os.remove("./data_py/data_{}.csv".format(filename))
    f = open("./data_py/data_{}.csv".format(filename), "a")
    header=['run','step','infected','protected','dead','exposure','viruscell','viruscontact','ifncell','ifncontact','ifn_prods']
    writer = csv.writer(f)
    writer.writerow(header)
    result=[]
    for istep in range(step):
        M_infected=Mtot[2]
        M_protected=Mtot[3]
        M_dead=Mtot[4]
        M_exposure=Mtot[5]
        M_exposure_timer=Mtot[6]
        M_viruscell = Mtot[9]
        M_virus_contact = Mtot[10]
        M_ifncell = Mtot[11]
        M_ifn_contact = Mtot[12]
        
        infected_count = np.sum(M_infected)
        protected_count = np.sum(M_protected)
        dead_count = np.sum(M_dead)
        exposure_count= np.sum(M_exposure)
        viruscell_count = np.sum(M_viruscell)
        viruscontact_count=np.sum(M_virus_contact)
        ifncell_count = np.sum(M_ifncell)
        ifncontact_count = np.sum(M_ifn_contact)
        data=[irun,istep,infected_count,protected_count,dead_count,exposure_count,viruscell_count,viruscontact_count,ifncell_count,ifncontact_count,ifn_prods]
        writer.writerow(data)
        result.append(data)


        if IsImage == 1:
            
            M_exposure_timer_exceed = 1*(M_exposure_timer>virus_prod_delay)
            M_color = np.stack((M_infected,M_infected-M_exposure_timer_exceed,M_protected), axis=-1)
            M_color = 1*(M_color>0).astype(int)*255   
            fig = plt.figure()
            plt.imshow(M_color)
            plt.title("Red:infected; Green:exposed; Blue:protected")
            plt.text(np.round(grid_size[1]*0.9).astype(int),np.round(grid_size[0]*0.1).astype(int),'{}'.format(istep), color='red',backgroundcolor='aliceblue')
            plt.savefig("./figure_py/infection_{}_{}.png".format(filename,istep))
            #plt.show()
            plt.close()
          
            M_color = np.stack((M_viruscell,M_virus_contact,M_dead), axis=-1)
            M_color = 1*(M_color>0).astype(int)*255
            fig = plt.figure()
            plt.imshow(M_color)
            plt.title("Red:virus generating; Green:virus contacting; Blue:dead")
            plt.text(np.round(grid_size[1]*0.9).astype(int),np.round(grid_size[0]*0.1).astype(int),'{}'.format(istep), color='red',backgroundcolor='aliceblue')
            plt.savefig("./figure_py/virus_{}_{}.png".format(filename,istep))
            #plt.show()
            plt.close()

            M_color = np.stack((M_ifncell,M_ifn_contact,M_dead), axis=-1)
            M_color = 1*(M_color>0).astype(int)*255
            fig = plt.figure()
            plt.imshow(M_color)
            plt.title("Red:IFN generating; Green:IFN contacting; Blue:dead")
            plt.text(np.round(grid_size[1]*0.9).astype(int),np.round(grid_size[0]*0.1).astype(int),'{}'.format(istep), color='red',backgroundcolor='aliceblue')
            plt.savefig("./figure_py/ifn_{}_{}.png".format(filename,istep))
            #plt.show()
            plt.close()
            
        if IsStochastic==1:
            Mtot=matrix_iter_stochastic(Mtot,Par1)
        else:
            Mtot=matrix_iter_homogenous(Mtot,Par1)
 
    df = pd.DataFrame(result, columns =header)
    f.close()
    return df

# Loop Runs
def main():
    args = sys.argv[1:]
    ifn_prods =int(args[0])
    s1=int(args[1])
    s2=int(args[2])
    virus_numb=int(args[3])
    IsStochastic=int(args[4])
    TotalStep = int(args[5])

    TotalRun = 100
    IsImage=0
    n1=100
    n2=100
    grid_size=[n1,n2]
    virus_size=[s1,s2]
    x0=49
    y0=49
    virus_center=[x0,y0]
    a=[0,1]
    virus_prods=2
    virus_prods_n = 0.5
    virus_prods_p = 1-virus_prods/(virus_prods+virus_prods_n)
    virus_diff=1
    
    ifn_diff=5
    ifn_prob = 0.01
    virus_reduction_factor=0.5
    prob_infect = float(0.2)

    virus_prod_delay=12
    ifn_prod_delay=12
    lifespan_mu = 60
    lifespan_sigma = 10
        
    Par0 = [grid_size,virus_size,virus_numb,virus_center,a,
    virus_prods,virus_prods_n,virus_prods_p,virus_diff, 
    ifn_prods, ifn_diff,ifn_prob,
    virus_reduction_factor,prob_infect,virus_prod_delay, 
    ifn_prod_delay, lifespan_mu,lifespan_sigma]
    print(Par0)
    for irun in range(0,TotalRun):
        df=Run(irun,TotalStep,Par0,IsImage,IsStochastic)

if __name__ == "__main__":
    main()
