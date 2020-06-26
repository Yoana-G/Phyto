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
    d = kappa*dt/dz**2
    # Compute arrays for diffusion
    Njp1 = np.append(Nj[1:],N0) # Dirichlet boundary condition (fixed value)
    Njm1 = np.append(Nj[0],Nj[:-1]) # Neumann boundary condition (no-flux)
    # Time stepping
    Nj = d*Njp1+(1-2*d)*Nj+d*Njm1
    
    return Nj
