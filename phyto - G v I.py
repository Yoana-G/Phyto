# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 21:14:40 2020

@author: yguzmanhernandez
"""
import numpy as np
import matplotlib.pyplot as plt
 
N0 = 10;     #Nutrient initial value
N = N0;
P0 = 10;     #Phytoplankton intial value
PE = P0;
#GP = 60;     #max growth rate (1 phyto per day) ; 0.001 - 10 "realistic range"
#k = 0.001; #0.05   #half-saturation parameter (starting w/half of 10) ; 0.01 - 5 "realistic range"
kpar = 0.05;     # light decay 
I_0 = 1;  #light intensity (I/I_o)
m = 0.1;  #0.01   #mortality rate ("slow") ; << G_max
mPS = 0.1;   #"quadratic mortality" is an option
z = np.linspace(-250, 0, 15)
#z = np.linspace(-200, 0, 2000)
I = I_0*np.exp(kpar*z) 
# np.arange(10000)
ntime = np.arange(10)   # of timesteps; create time array of 1000 components
dt = 0.01          # time step of (1 day) 
days = dt*np.arange(10)   #outputs number of days; multiplies 'ntime*dt'
#days = dt*np.arange(10000)

k_PS = 0.25;  #half-saturation parameter; 0.01 - 5 "realistic range" (0.05)
k_PE = 0.5;


L_PE = np.linspace(250, 50, 15)
L_PS = L_PE


##### Depth v Phyto Growth #####
fig1, ax = plt.subplots()
for L_PE in [250, 50, 5]:
    GP = 2.8;    #60
    GS = 1.5;    #15
    L_PS = L_PE
    G_PE = GP*(1-np.exp(-L_PE*I/I_0))
    G_PS = GS*(1-np.exp(-L_PS*I/I_0))
    N_PS = ((mPS/G_PS)*k_PS)/(1 - (mPS/G_PS))
    N_PE = ((m/G_PE)*k_PE)/(1 - (m/G_PE))
    line1, = ax.plot(N_PE, z, dashes=[6, 2], label='L_PE=' +str(L_PE))
    line2, = ax.plot(N_PS, z, dashes=[6, 2], label='L_PS=' +str(L_PS))
    

ax.set_xlim((0, 4))    # setting the x axis limits to 10^-4
#ax.set_xlim((0, 0.0005)) 
ax.legend()

plt.xlabel('N* (at equilibrium)')
plt.ylabel('Depth (z)')


fig2, ax = plt.subplots()

line1, = ax.plot(G_PE, z, label='Picoeukaryotes (PE)')
#line1.set_dashes([2, 2, 10, 2])   # 2pt line, 2pt break, 10pt line, 2pt break
line2, = ax.plot(N_PE, z, label='N*_PE')
#line2, = ax.plot(G_PS, z, dashes=[6, 2], label='Synechococcus (PS)')
line3, = ax.plot(N_PS, z, label='N*_PS')

ax.set_xlim((0, 1))
ax.legend()
#plt.plot(G_PE, I)
plt.xlabel('Nutrients and Phytoplankton')
#plt.xlabel('Growth rate of Phytoplankton')
plt.ylabel('Depth(z)')


############## Light intensity (I) v. L_PE ################
'''
L_PE = 200;
G_PE = GP*(1-np.exp(-L_PE*I/I_0))  
fig, ax = plt.subplots()

line1, = ax.plot(G_PE, I, dashes=[6, 2], label='L_PE = 15')
L_PE = 100
G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
line2, = ax.plot(G_PE, I, dashes=[6, 2], label='L_PE = 13')
L_PE = 80
G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
line3, = ax.plot(G_PE, I, dashes=[6, 2], label='L_PE = 11')
L_PE = 60
G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
line4, = ax.plot(G_PE, I, dashes=[6, 2], label='L_PE = 9')
L_PE = 40
G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
line5, = ax.plot(G_PE, I, dashes=[6, 2], label='L_PE = 7')
####### Synechococcus ######
L_PS = 150;
G_PS = GS*(1-np.exp(-L_PS*I/I_0))
line6, = ax.plot(G_PS, I, dashes=[6, 2], label='L_PS = 15')
L_PS = 110
G_PS = GS*(1-np.exp(-L_PS*I/I_0)) 
line7, = ax.plot(G_PS, I, dashes=[6, 2], label='L_PS = 13')
L_PS = 90
G_PS = GS*(1-np.exp(-L_PS*I/I_0)) 
line8, = ax.plot(G_PS, I, dashes=[6, 2], label='L_PS = 11')
L_PS = 70
G_PS = GS*(1-np.exp(-L_PS*I/I_0)) 
line9, = ax.plot(G_PS, I, dashes=[6, 2], label='L_PS = 9')
L_PS = 30
G_PS = GS*(1-np.exp(-L_PS*I/I_0)) 
line10, = ax.plot(G_PS, I, dashes=[6, 2], label='L_PS = 7')
plt.title('Light sensitivity of Phytoplankton')
ax.legend()
plt.xlabel('Growth rate of Phytoplankton')
plt.ylabel('Light intensity(I)')
'''

''' 
L_PE = 15;
#G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
G_PE = PE*(GP*(1-np.exp(-L_PE*I))*(N/(k+N))*I - m) 
fig, ax = plt.subplots()

line1, = ax.plot(G_PE, z, dashes=[6, 2], label='L_PE = 15')
L_PE = 13
G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
line2, = ax.plot(G_PE, z, dashes=[6, 2], label='L_PE = 13')
L_PE = 11
G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
line3, = ax.plot(G_PE, z, dashes=[6, 2], label='L_PE = 11')
L_PE = 9
G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
line4, = ax.plot(G_PE, z, dashes=[6, 2], label='L_PE = 9')
L_PE = 7
G_PE = GP*(1-np.exp(-L_PE*I/I_0)) 
line5, = ax.plot(G_PE, z, dashes=[6, 2], label='L_PE = 7')
'''
####### Synechococcus ######
'''
L_PS = 15;
G_PS = GS*(1-np.exp(-L_PS*I/I_0))
line6, = ax.plot(G_PS, z, dashes=[6, 2], label='L_PS = 15')
L_PS = 13
G_PS = GS*(1-np.exp(-L_PS*I/I_0)) 
line7, = ax.plot(G_PS, z, dashes=[6, 2], label='L_PS = 13')
L_PS = 11
G_PS = GS*(1-np.exp(-L_PS*I/I_0)) 
line8, = ax.plot(G_PS, z, dashes=[6, 2], label='L_PS = 11')
L_PS = 9
G_PS = GS*(1-np.exp(-L_PS*I/I_0)) 
line9, = ax.plot(G_PS, z, dashes=[6, 2], label='L_PS = 9')
L_PS = 7
G_PS = GS*(1-np.exp(-L_PS*I/I_0)) 
line10, = ax.plot(G_PS, z, dashes=[6, 2], label='L_PS = 7')
plt.title('Light sensitivity of Phytoplankton')
ax.legend()
plt.xlabel('Growth rate of Phytoplankton')
plt.ylabel('Depth(z)')
'''