from Model import *

#%%
###################################################
#                     Test
###################################################   

mus   =  np.linspace(-0.05, 0.2, 51)
Smags =  np.linspace(0, 0.05, 51)

params_TI['mu_bulk'] = 0
params_TI['mu_lead1'] = 0.042
params_TI['mu_lead2'] = 0.042

params_TI["B_x"] = 0.500712/(W_y*H_z)    
momenta = np.linspace(-0.2, 0.2, 101)
plot_bands(systf1, momenta, levels=10, params=params_TI)

#r1 = get_bands(systf1, [0,1], levels=8, params=params_TI)

#print(r1)

#%%

############################################################################
#               Finding bad points for 0 disorder
############################################################################

mu_lead = 0.042
B = 0.500001

params_TI['S_mag'] = 0
params_TI['mu_bulk'] = 0
params_TI['mu_lead1'] = 0.042
params_TI['mu_lead2'] = 0.042
params_TI["B_x"] = B/(W_y*H_z)  

momenta = np.linspace(-0.2, 0.2, 101)
plot_bands(systf1, momenta, levels=20, params=params_TI)

energy_k_0 = get_bands(systf1, [0], levels = 20, params = params_TI)

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead

mus =  np.linspace(-0.05, 0.2, 51)
"""
for mu in energy_k_0[0,::2]:
    if mu > 0:
        params_TI['mu_bulk'] = mu_lead + mu
        smat_e = scattering_matrix(systf1, params_TI)[0]
        #df = pd.DataFrame(smat_e)
        #display(df)
        print(f"For {mu} the max value in TRS is {np.amax(np.absolute(smat_e + smat_e.T))}")
"""
#print(energy_k_0[:,::2])

#%%

############################################################################
#               Optimization of exact value for half flux
############################################################################

def Reward(val1, val2):
    if abs(val1 - val2) < 10e-6 and val1 != 0:
        return 1000
    else:
        return abs(val1 - val2)
    
number_bands = 3           # depending on lead mu, should be odd
levels = number_bands*2
params_TI['S_mag'] = 0     # shouldnt matter
params_TI['mu_bulk'] = 0   # shouldnt matter
params_TI['mu_lead1'] = 0  # do not touch
params_TI['mu_lead2'] = 0  # same
#params_TI["B_x"] = 0.5/(W_y*H_z)  

"""
energy_k_0 = get_bands(systf1, [0], levels = 5, params = params_TI)
print(energy_k_0)
print(energy_k_0[0,int(len(energy_k_0[0])/2)::2])
momenta = np.linspace(-0.2, 0.2, 51)
plot_bands(systf1, momenta, levels=20, params=params_TI)

"""
offsets = np.linspace(0, 1e-3, 1001)
best_offset = 0
best_reward = 10000000

time_start = time.time()

for j, offset in enumerate(offsets):
    params_TI["B_x"] = (0.5 + offset)/(W_y*H_z)
    energy_k_0 = get_bands(systf1, [0], levels = levels, params = params_TI)[0]
    energy_k_0 = energy_k_0[len(energy_k_0)//2 : : 2]
    energy_k_0 = energy_k_0[:number_bands]
    total_reward = 0
    for i in range(0, number_bands, 2):
        if i == 0:
            total_reward += 2*Reward(0, energy_k_0[i])
        else:
            total_reward += Reward(energy_k_0[i-1], energy_k_0[i])
    if total_reward < best_reward:
        best_reward = total_reward
        best_offset = offset
    
    time_current = time.time()
    percentage = (j + 1)/(len(offsets))
    time_left = (time_current - time_start)/percentage - (time_current - time_start) 
    print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        
if best_reward > 100:
    raise Exception("It dont work")

print("Min reward : ", best_reward, " offset : ", best_offset)
params_TI["B_x"] = (0.5 + best_offset)/(W_y*H_z)
momenta = np.linspace(-0.2, 0.2, 51)
plot_bands(systf1, momenta, levels=levels, params=params_TI)

    
    



