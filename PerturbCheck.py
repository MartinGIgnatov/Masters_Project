#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 11:57:05 2022

@author: leonardobossi1
"""
from Model import *

# Testing the maximum difference in deivation between the absolute value of the scattering matrices
# to see if the magnitude of the perturbation will affect this value

def probsymm(matrix):
    prob_mat = np.abs(matrix) ** 2
    diff = prob_mat - prob_mat.T
    
    return abs(np.amax(diff))
      
        
params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead

mu_test = np.linspace(-0.05, 0.2, 11)
Smag_test = np.linspace(0, 0.05, 11)
B_dev= np.array([1e-3, 1e-4,  1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10])

sim_data = []

for deviation in B_dev:
    params_TI["B_x"] = deviation/(W_y*H_z)   # half flux is 0.5/(W_y*H_z)  
    print(deviation)
    for idx in range(len(mu_test)):
        print(idx/len(mu_test))
        mu_bulk = mus[get_closest_index(mu_test[idx], mus)] + mu_lead
        params_TI['mu_bulk']  = mu_bulk
        
        for idx2 in range(len(Smag_test)):
        
            smag = Smags[get_closest_index(Smag_test[idx2] , Smags)]
            path = path_generator.generate_path(["Data","BPerturb_Matrices"],f"Scattering_Matrices_{mu_lead}_flux_{deviation}"f"mu_bulk_{mu_bulk}_S_mag_{smag}","txt")
            params_TI['S_mag']= smag
            
            smat_e = scattering_matrix(systf1, params_TI, calibration=None)[0]
            max_diff = probsymm(smat_e)
            
            sim_data.append([deviation, mu_bulk, smag, max_diff])
            
            np.savetxt(fname =path, X = smat_e)
        
