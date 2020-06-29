# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 14:00:25 2020

@author: yguzmanhernandez
"""
import numpy as np
import matplotlib.pyplot as plt        #import plotting library
from diffusion_funct_defs import diffusion_time_stepping  #import all files in the same folder called 'diffusion func defs'


H = -1000; # depth of water (1 = non-dimensional)
K = 50; # number of grid points

GP = 62.5;  #max growth rate (1 phyto per day) ; 0.001 - 10 "realistic range"
GS = 2.5;
I_0 = 1;   #light intensity (I/I_o); constant (doesn't vary w/season); Watts*m^-2 (light flux)
L_PE = 0.01; #gamma;light dependence =>determines Growth rate sensitivity to light levels
L_PS = 0.001; #unitless (it's in the exponential)
k = 0.5;  #nutrient units (micromolar)
k_PS = 0.05;  #half-saturation parameter; 0.01 - 5 "realistic range"
k_PE = 0.01;  #larger 
m = 0.1;  #mortality rate ("slow") ; << G_max
mPS = 0.1; #"quadratic mortality" is an option


z = np.linspace(0,H,K) # z array
I = I_0*np.exp(K*z) #light profile (light decays exponentially)
kappa = 0.0001 * 86400 # diffusivity value (10^4 m^2 * s^-1); 86400 secs/day 
dz = H/(K-1); # calculate grid spacing
N0 = 20; # bottom boundary condition on nutrient
N = N0*np.ones(z.shape) # initial condition
PE = np.ones(z.shape) #Picoeukaryote vector
PS = np.ones(z.shape) #Synechococcus vector
PEsave = [PE] #saving initial condition
PSsave = [PS]
Nsave = [N] # save initial condition
dt = 0.01 # time step
nt = 10000 # number of time steps
time = np.arange(nt+1)*dt

#G_PE = np.arange(1, 50, 0.1)  #vector for Synechococcus growth rate

#N_PE = (((m/G_PE)*k_PE)/(1-(m/G_PE)))
#N_PS = (((mPS/G_PS)*k_PS)/(1-(m/G_PS)))


for t in np.arange(nt): #for loop: iterate in time, run specified # of times (nt)

    #time-stepping biology
    dN_dt = -GP*(1-np.exp(-L_PE*I))*(N/(k+N))*I*PE + m*PE - GS*(1-np.exp(-L_PS*I))*(N/(k+N))*I*PS + mPS*PS 
    
    dPE_dt= PE*(GP*(1-np.exp(-L_PE*I))*(N/(k+N))*I - m)  #produces time scale to reach equilibrium
    
    dPS_dt= PS*(GS*(1-np.exp(-L_PS*I))*(N/(k+N))*I - mPS)   #produces time scale to reach equilibrium
    
    N = N + dN_dt*dt
    PE = PE + dPE_dt*dt
    PS = PS + dPS_dt*dt
      
#time-stepping diffusion
    N = diffusion_time_stepping(kappa,dz,dt,N0,N) 
    Nsave = np.concatenate((Nsave,[N]),axis = 0) #adding a row each time   
    
    PE = diffusion_time_stepping(kappa, dz, dt, 0, PE)
    PEsave = np.concatenate((PEsave,[PE]),axis = 0) #adding [PE] to the matrix, starting w/initial value 
                                                    # axis tells dimension; zero = rows
    PS = diffusion_time_stepping(kappa, dz, dt, 0, PS)
    PSsave = np.concatenate((PSsave,[PS]),axis=0)


'''fig, ax = plt.subplots()

plt.plot(grid, Nsave[10000,:], label= 'Nutrients (N)')

line1, = ax.plot(time, Nsave[:,0]) #label='N')
line2, = ax.plot(time, Nsave[:,-1], dashes=[6, 2], #label='N')

ax.legend()
plt.xlabel('N')
plt.ylabel('z(m)')

'''
plt.pcolor(time,z,Nsave.T)
plt.colorbar()
plt.xlabel('time in days')
plt.ylabel('depth of water column')

plt.figure()
plt.pcolor(time,z,PEsave.T)
plt.colorbar()
plt.xlabel('time in days')
plt.ylabel('depth of water column')

plt.figure()
plt.pcolor(time,z,PSsave.T)
plt.colorbar()
plt.xlabel('time in days')
plt.ylabel('depth of water column')

'''
#plotting profiles w/depth
plt.plot(Nsave[10000,:], grid)
plt.xlabel('N')
plt.ylabel('z(m)')

#plotting time series @particular depth
plt.plot(time,Nsave[:,0])
plt.plot(time,Nsave[:,-1])
plt.xlabel('time')
plt.ylabel('N')
'''
#plt.plot(time,Nsave);
