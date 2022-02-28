from Model import *

#%%

    

params_TI['mu_bulk'] = 0.042
params_TI['mu_lead1'] = 0.042
params_TI['mu_lead2'] = 0.042
params_TI["B_x"] = 0.5/(W_y*H_z)    
momenta = np.linspace(-0.2, 0.2, 101)
plot_bands(systf1, momenta, levels=4, params=params_TI)
r1 = get_bands(systf1, [0,1], levels=4, params=params_TI)

print(r1)