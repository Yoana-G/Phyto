# -*- codinG_max: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
 

N0 = 10;       #Nutrient initial value
P0 = 10;       #Phytoplankton intial value 
PS = 10;
GP = 62.5;      #max growth rate (1 phyto per day) ; 0.001 - 10 "realistic range"
GS = 2.5;
I = 1;          #light intensity (I/I_o)
k = 0.5;          #half-saturation parameter (starting w/half of 10) ; 0.01 - 5 "realistic range"
m = 0.1;       #mortality rate ("slow") ; << G_max
mPS = 0.1;   #"quadratic mortality" is an option
N = N0;
PE = P0;
PS = PS;
L_PE = 0.01;    #gamma;light dependence =>determines Growth rate sensitivity to light levels
L_PS = 0.001;

PSsave = np.empty(0)    #creates new array without initializing entries;random values
Psave = np.empty(0)     #creates new array without initializing entries;random values
Nsave = np.empty(0)     #creates new array without initializing entries
time = np.arange(10000)   #1000 timesteps of 0.01; create time array of 1000 components
dt = 0.01          # time step of (1 day) 
days = dt*np.arange(10000)   #outputs number of days; multiplies 'time*dt'


for t in time:          #iterate in time

    dN_dt = -GP*(1-np.exp(-L_PE*I))*(N/(k+N))*I*PE + m*PE - GS*(1-np.exp(-L_PS*I))*(N/(k+N))*I*PS + mPS*PS 
    
    dPE_dt= PE*(GP*(1-np.exp(-L_PE*I))*(N/(k+N))*I - m)  #produces time scale to reach equilibrium
    
    dPS_dt= PS*(GS*(1-np.exp(-L_PS*I))*(N/(k+N))*I - mPS)   #produces time scale to reach equilibrium
    
    N = N + dN_dt*dt
    PE = PE + dPE_dt*dt
    PS = PS + dPS_dt*dt
    
    Psave = np.append(Psave,PE)      #adds Picoeukaryote (P) to the end of the list
    Nsave = np.append(Nsave,N)      #adds N to the end of the list
    PSsave = np.append(PSsave,PS)   #adds Synechococcus (PS) phyto. to the end of the list

fig, ax = plt.subplots()

plt.plot(days,Nsave, label= 'Nutrients (N)')

line1, = ax.plot(days, Psave, label='Picoeukaryotes (P)')
line1.set_dashes([2, 2, 10, 2])     # 2pt line, 2pt break, 10pt line, 2pt break
line2, = ax.plot(days, PSsave - 0.2, dashes=[6, 2], label='Synechococcus (PS)')

ax.legend()
plt.show()
plt.xlabel('Days')
plt.ylabel('Growth rate (d^-1)')

#print(Nsave[-1])

    
# Does changing mortality rate closer to G_max/growth rate drive phyto. to zero (or non-zero N)?
# Does changing mortality rate of 1 species drive nutrient(s) availability up or down? How does that change
#       other phytoplankton concentration up or down?

#Conclusion 1: Higher sensitivity to light changes (L), bigger mortality rate (m), 
## faster growth rate (G), and k = 2 is dominant in a particular plot. --> high light ecotype(?)

