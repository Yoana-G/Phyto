# -*- codinG_max: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
 

N0 = 10;     #Nutrient initial value
P0 = 10;     #Phytoplankton intial value 
PR = 10;
GS = 5;     #max growth rate (1 phyto per day) ; 0.001 - 10 "realistic range"
GR = 1.5;
I = 1;     #light intensity (I/I_o)
k = 1;     #half-saturation parameter (starting w/half of 10) ; 0.01 - 5 "realistic range"
mPS = 0.1;     #mortality rate ("slow") ; << G_max
mPR = 0.1;   #0.1 #"quadratic mortality" is an option
N = N0;
PS = P0;
PR = PR;
L_PS = 0.5; #gamma;light dependence =>determines Growth rate sensitivity to light levels
L_PR = 250;

PRsave = np.empty(0)    #creates new array without initializing entries;random values
PSsave = np.empty(0)     #creates new array without initializing entries;random values
Nsave = np.empty(0)     #creates new array without initializing entries
ntime = np.arange(100000)   # of timesteps; create time array of 1000 components
dt = 0.01          # time step of (1 day) 
days = dt*np.arange(100000)   #outputs number of days; multiplies 'ntime*dt'


for t in ntime:          #iterate in time

    dN_dt = -GS*(1-np.exp(-L_PS*I))*(N/(k+N))*I*PS + mPS*PS - GR*(1-np.exp(-L_PR*I))*(N/(k+N))*I*PR + mPR*PR 
    
    dPS_dt= PS*(GS*(1-np.exp(-L_PS*I))*(N/(k+N))*I - mPS)   #produces time scale to reach equilibrium
    
    dPR_dt= PR*(GR*(1-np.exp(-L_PR*I))*(N/(k+N))*I - mPR) #produces time scale to reach equilibrium
    
    N = N + dN_dt*dt
    PS = PS + dPS_dt*dt
    PR = PR + dPR_dt*dt


    PSsave = np.append(PSsave,PS)     #adds Picoeukaryote (P) to the end of the list
    Nsave = np.append(Nsave,N)      #adds N to the end of the list
    PRsave = np.append(PRsave,PR)   #adds Synechococcus (PS) phyto. to the end of the list

fig, ax = plt.subplots()

plt.plot(days, Nsave, label= 'Nutrients (N)')

line1, = ax.plot(days, PSsave, label='Synechococcus (PS)')
line1.set_dashes([2, 2, 10, 2])   # 2pt line, 2pt break, 10pt line, 2pt break
line2, = ax.plot(days, PRsave, dashes=[6, 2], label='Prochlorococcus (PR)')

ax.legend()
#plt.show()
plt.xlabel('Time (days)')
plt.ylabel('Growth rate (d^-1)')

#print(Nsave[-1])

    
# Does changing mortality rate closer to G_max/growth rate drive phyto. to zero (or non-zero N)?
# Does changing mortality rate of 1 species drive nutrient(s) availability up or down? How does that change
#       other phytoplankton concentration up or down?

#Conclusion 1: Higher sensitivity to light changes (L), bigger mortality rate (m), 
## faster growth rate (G), and k = 2 is dominant in a particular plot. --> high light ecotype(?)