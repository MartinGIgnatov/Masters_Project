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



#%%
# Band structure

def plot_bands(syst, momenta, return_bands=False, plot_bands=True, params=None, levels=0):
    bands = kwant.physics.Bands(syst.leads[0], params=params)
    energies = np.array([bands(k) for k in momenta])
    if levels!=0:
        mid = len(energies[0])//2
        a=energies[:,mid-levels:mid]
        b=energies[:,mid:mid+levels]
        energies = np.concatenate((a,b),axis=1)
    if plot_bands:
        fig, ax = plt.subplots()
        ax.set_xlabel('k')
        ax.set_ylabel('Energy')
        ax.plot(momenta, energies)
        #ax.plot(momenta, np.full(len(momenta), 0.005), 'b--')
        #ax.plot(momenta, np.full(len(momenta), 0.015), 'r--')
        plt.show()
    if return_bands:
        return energies
    
flux = '0 - perturbed001'
params_TI["B_x"]     = (5.1)/(W_y*H_z)
params_TI['mu_bulk'] = 0.02
params_TI['mu_lead1'] = 0.042
params_TI['mu_lead2'] = 0.042
momenta = np.linspace(-1, 1, 101)

## flux: 001, 0, 501, 0.5


# Have many levels but only bottom energy levels are relevant
r1 = plot_bands(systf1, momenta, return_bands=True, levels=10, params=params_TI)
plt.title(flux)
plt.ylim(0.04, 0.05)
plt.xlim(-0.25, 0.25)
print(r1)