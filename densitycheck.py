#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 15:41:16 2022

@author: leonardobossi1
"""

from Model import *

#%%

syst = systf1

wf = kwant.wave_function(syst, energy=0, params=params_TI) # energy in a band structure

# one index refers to the lead, the other refers to the mode
psi = wf(0)[0]
psi_2 = wf(1)[0]

J_m = kwant.operator.Current(syst)

# Psi* A Psi where A is the current operator

m_current = J_m(psi, params = params_TI)

m2_current = J_m(psi_2, params = params_TI)

#kwant.operator.Current(system, onsite = 1)

D_m = kwant.operator.Density(syst)
m_density = D_m(psi, params = params_TI)
m2_density = D_m(psi_2, params = params_TI)


