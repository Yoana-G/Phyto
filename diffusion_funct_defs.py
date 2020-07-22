# -*- coding: utf-8 -*-
import numpy as np


def diffusion_time_stepping(kappa,dz,dt,N0,N):
    '''Solve diffusion equation in 1D using forward Euler time stepping
    
    Inputs:
    kappa:           float: diffusivity
    dz:              float: grid spacing
    dt:              float: time step
    N0:              float: bottom boundary value of nutrient
    N:               float: array of values at current time-step
    
    Output:
    Nj:              float: array of values at next time-step
    '''
    Nj = N
    d = kappa*dt/dz**2  #second derivative
    #double asterisk is an exponent (second derivative in this case)
   
    # Compute arrays for diffusion
    Njp1 = np.append(Nj[1:],N0) # Dirichlet boundary condition (fixed value)
    Njm1 = np.append(Nj[0],Nj[:-1]) # Neumann boundary condition (no-flux)
    
    # Time stepping
        # the minus 2 comes from the stepping process
        # the 1 comes from explicit method (Euler's)
    Nj = d*Njp1+(1-2*d)*Nj+d*Njm1   #second derivative
    return Nj

def sinking_time_stepping(dz,N0,PE):
    PEj = PE
    d = 1 / dz
    
    # Compute arrays for sinking velocity
    #PEjp1 = np.append(PEj[1:],PE0) # Dirichlet boundary condition (fixed value)
    PEjp1 = np.append(PEj[0],PEj[:-1])
    #Njm1 = np.append(Nj[0],Nj[:-1]) # Neumann boundary condition (no-flux)
    
    # Time stepping
    #Nj = d*Njp1+(1-2*d)*Nj+d*Njm1   #second derivative
    PEj = d*PEjp1-d*PEj
    return PEj

