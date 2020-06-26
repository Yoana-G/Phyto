# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 14:00:25 2020

@author: yguzmanhernandez
"""
import numpy as np
import matplotlib.pyplot as plt        #import plotting library
from diffusion_funct_defs import diffusion_time_stepping  #import all files in the same folder called 'diffusion func defs'


H = -1; # depth of water (1 = non-dimensional)
K = 20; # number of grid points

GP = 62.5;  #max growth rate (1 phyto per day) ; 0.001 - 10 "realistic range"
GS = 2.5;
I = 1;   #light intensity (I/I_o)
L_PE = 0.01; #gamma;light dependence =>determines Growth rate sensitivity to light levels
L_PS = 0.001;
P0 = 10;  #Phytoplankton intial value 
PS = 10;
PE = P0;
PS = PS;
k = 0.5;  #half-saturation parameter (starting w/half of 10) ; 0.01 - 5 "realistic range"
m = 0.1;  #mortality rate ("slow") ; << G_max
mPS = 0.1; #"quadratic mortality" is an option

grid = np.linspace(0,H,K) # grid array
kappa = 0.01 # diffusivity value
dz = H/(K-1); # calculate grid spacing
N0 = 1; # bottom boundary condition on nutrient
N = np.zeros(grid.shape) # initial condition
Nsave = [N] # save initial condition
dt = 0.01 # time step
nt = 10000 # number of time steps
time = np.arange(nt+1)*dt


for t in np.arange(nt):

    #biological equations 
    dN_dt = -GP*(1-np.exp(-L_PE*I))*(N/(k+N))*I*PE + m*PE - GS*(1-np.exp(-L_PS*I))*(N/(k+N))*I*PS + mPS*PS 
    
    dPE_dt= PE*(GP*(1-np.exp(-L_PE*I))*(N/(k+N))*I - m)  #produces time scale to reach equilibrium
    
    dPS_dt= PS*(GS*(1-np.exp(-L_PS*I))*(N/(k+N))*I - mPS)   #produces time scale to reach equilibrium
    
      
    N = diffusion_time_stepping(kappa,dz,dt,N0,N)
    Nsave = np.concatenate((Nsave,[N]),axis = 0)


plt.pcolor(time,grid,Nsave.T)
plt.colorbar()
plt.xlabel('time in days')
plt.ylabel('depth of water column')

plt.plot(Nsave[10000,:], grid)
plt.xlabel('N')
plt.ylabel('z(m)')

plt.plot(time,Nsave[:,0])
plt.plot(time,Nsave[:,-1])
plt.xlabel('time')
plt.ylabel('N')

#plt.plot(time,Nsave);
