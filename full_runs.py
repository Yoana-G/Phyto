# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 11:47:49 2020

@author: yguzmanhernandez
"""

import xarray as xr # import the package

DS = xr.open_dataset('full_00000.cdf') # open dataset

print(DS.zc[:,1,1]) # print depths of cell centers at one x and y location


