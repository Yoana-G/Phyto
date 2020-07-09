# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 10:35:58 2020

@author: yguzmanhernandez
"""
import numpy as np
import matplotlib.pyplot as plt        #import plotting library


H = -1000; # depth of water (1 = non-dimensional)
kpar = 0.05; #20 # number of grid points; @20m half the light; photosynthetically active radiation
#GP = 62.5;  #max growth rate (1 phyto per day) ; 0.001 - 10 "realistic range"
I_0 = 1;   #light intensity (I/I_o); constant (doesn't vary w/season); Watts*m^-2 (light flux)
L_PE = 5; #gamma;light dependence =>determines Growth rate sensitivity to light levels
L_PS = 250; #unitless (it's in the exponential)
k_PS = 0.25;  #half-saturation parameter; 0.01 - 5 "realistic range" // (0.5), (0.05)
k_PE = 0.5;  #larger (m^2 / s); =0.01
m = 0.1; #mortality rate ("slow") = day^-1 ; << G_max
mPS = 0.1; #"quadratic mortality" is an option


z = 50; #50m (Depth) crossover
I = I_0*np.exp(kpar*z) #light profile (light decays exponentially)
kappa = 0.0001 * 86400 # diffusivity value (10^4 m^2 * s^-1); 86400 secs/day; kappa = 0.0001
G_PS0 = np.arange(0, 30)


G_PS = G_PS0*(1 - np.exp(-L_PS*I/I_0))
G_PE0 = ((((G_PS*k_PE*m)-(mPS*m*k_PE)) / (mPS*k_PS) ) + m) / (1 - np.exp(-L_PS*I/I_0))

plt.plot(G_PS0,G_PE0)
plt.xlabel('G_PS0')
plt.ylabel('G_PE0')





