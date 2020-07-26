# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 14:00:25 2020

@author: yguzmanhernandez
"""
import numpy as np
import matplotlib.pyplot as plt        #import plotting library
from diffusion_funct_defs import diffusion_time_stepping  #import all files in the same folder called 'diffusion func defs'
from diffusion_funct_defs import sinking_time_stepping
from diffusion_funct_defs import diffusion2_time_stepping

#all PE --> PS  AND all PS --> PR

# H = -1000
H = 1000; # depth of water (1 = non-dimensional)
K = 1000; # number of grid points (#20)
kpar = 0.05; # light decay 
GS = 5;  #5 #max growth rate (1 phyto per day) ; 0.001 - 10 "realistic range"
GR = 1.5; #1.5
GE = 5;
I_0 = 1;   #light intensity (I/I_o); constant (doesn't vary w/season); Watts*m^-2 (light flux)
L_PS = 5; #gamma;light dependence =>determines Growth rate sensitivity to light levels (0.01)
L_PR = 250; #unitless (it's in the exponential) / (0.001)
L_PE = 5;
#k = 0.5;  #nutrient units (micromolar)
k_PR = 0.25;  #half-saturation parameter; 0.01 - 5 "realistic range" (0.05)
k_PS = 0.5;  #larger (m^2 / s); =0.01 ///// (0.001)
k_PE = 0.5;
mPS = 0.1; #mortality rate ("slow") = day^-1 ; << G_max
mPR = 0.1; #"quadratic mortality" is an option
mPE = 0.1;

z = np.linspace(0,-H,K) # z array
I = I_0*np.exp(kpar*z) #light profile (light decays exponentially)
kappa = 0.001 * 86400 # diffusivity value (m^2 * s^-1); 86400 secs/day; kappa = 0.0001
dz = H/(K-1); # calculate z spacing
N0 = 22; # bottom boundary condition on nutrient
N = (N0-2)*np.ones(z.shape) # initial condition
PS = np.ones(z.shape) #Synechococcus vector
PR = np.ones(z.shape) #Prochlorococcus vector
PE = np.ones(z.shape) #Picoeukaryote vector
PSsave = [PS] #saving initial condition
PRsave = [PR]
PEsave = [PE]
Nsave = [N] # save initial condition
dt = 0.001 # time step
nt = 100000 # number of time steps 
time = np.arange(nt+1)*dt   #won't work bc it saves each time
t_val = np.zeros(z.shape)
tsave = [t_val]
W = 1; #sinking rate
PE0 = 0; # bottom boundary condition for phytoplankton ??


for t in np.arange(nt): #for loop: iterate in time, run specified # of times (nt)

    #time-stepping biology
    dN_dt = -GS*(1-np.exp(-L_PS*I))*(N/(k_PS+N))*PS + mPS*PS - GR*(1-np.exp(-L_PR*I))*(N/(k_PR+N))*PR + mPR*PR 
    
    dPS_dt= PS*(GS*(1-np.exp(-L_PS*I))*(N/(k_PS+N)) - mPS)  #produces time scale to reach equilibrium
    
    dPR_dt= PR*(GR*(1-np.exp(-L_PR*I))*(N/(k_PR+N)) - mPR)   #produces time scale to reach equilibrium
    
    dPE_dz = sinking_time_stepping(dz,PE0,PE)
    #d*PEj-d*PEjm1
        
    dPE_dt = PE*(GE*(1-np.exp(-L_PE*I))*(N/(k_PE+N)) - mPE - W*dPE_dz)
    # sinking term --> W*dPE_dz
    
    N = N + dN_dt*dt      #forward in time (explicit method) - Euler's method
    PS = PS + dPS_dt*dt
    PR = PR + dPR_dt*dt
    PE = PE + dPE_dt*dt
    
      
#time-stepping diffusion
    N = diffusion_time_stepping(kappa,dz,dt,N0,N) 
    PS = diffusion_time_stepping(kappa,dz,dt,0,PS)
    PR = diffusion_time_stepping(kappa,dz,dt,0,PR)
    #PE = diffusion_time_stepping(kappa, dz, dt, 0, PE)
    PE = diffusion2_time_stepping(kappa,dz,dt,PE0,W,PE)

    if np.mod(t,1000) == 0:   #save every 10 time steps 
        Nsave = np.concatenate((Nsave,[N]),axis = 0)
        PSsave = np.concatenate((PSsave,[PS]),axis = 0) #adding [PE] to the matrix, starting w/initial value 
                                                        # axis tells dimension; zero = rows
        PRsave = np.concatenate((PRsave,[PR]),axis=0)
        PEsave = np.concatenate((PEsave,[PE]),axis=0)
        tsave = np.concatenate((tsave,[t*dt*np.ones(z.shape)]),axis=0)
        

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

fig1, ax = plt.subplots()
#plotting profiles w/depth
line1, = ax.plot(Nsave[-1,:], z, dashes=[6, 2], label='Nutrients')
#plt.plot(Nsave[10000,:], z)  #10,000
#plotting profiles w/depth --> PE
line2, = ax.plot(PSsave[-1,:], z, label='Synechococcus (PS)')
#plt.plot(PEsave[10000,:], z)   #10,000
##plotting profiles w/depth --> PS
line3, = ax.plot(PRsave[-1,:], z, label='Prochlorococcus (PR)')
line4, = ax.plot(PEsave[-1,:], z, label='Picoeukaryotes (PE)')

ax.legend()
plt.xlabel('N (uMol/m^3)')  #µ
plt.ylabel('z (meters)')


fig2, ax = plt.subplots()
#plotting time series @particular depth
#plt.plot(time,Nsave[:,0], label='Nutrients')
line1, = ax.plot(tsave[:,0], Nsave[:,0], label='Nutrients')
#line2, = ax.plot(time,Nsave[:,-1], label='Nutrients(at H)')
#plt.plot(time,Nsave[:,-1])
line3, = ax.plot(tsave[:,0], PSsave[:,0], label='Synechococcus (PS)')
#plotting time series @particular depth
#plt.plot(time,PEsave[:,0], label='PE')
#line4, = ax.plot(time,PEsave[:,-1], label='PE (at H)')
#plotting time series @particular depth
line5, = ax.plot(tsave[:,0], PRsave[:,0], label='Prochlorococcus (PR)')
#line6, = ax.plot(time, PSsave[:,-1], label='PS (at H)')
line6, = ax.plot(tsave[:,0], PEsave[:,0], label='Picoeukaryotes (PE)')

ax.legend()
plt.xlabel('time (day^-1)')
plt.ylabel('N (uMol/m^3)')
