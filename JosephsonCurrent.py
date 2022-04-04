# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 11:15:49 2022

@author: marti
"""

from Model import *


#%%
######################################################################
#########   Generate ABS values and store                 ############   useless
######################################################################

mus   =  np.linspace(-0.05, 0.2, 51)
Smags =  np.linspace(0, 0.05, 51)

mu_lead = 0.042    # 0.02, 0.042
phases = np.linspace(-0.1, 1.1, 121)

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead

params_TI["B_x"]      = 0.500712/(W_y*H_z)   # half flux is 0.5/(W_y*H_z)
flux = "500712"       # quarter, half, zero

#mus = mus[int(len(mus)*0.71):] 

time_start = time.time()

for i, mu in enumerate(mus):
    params_TI['mu_bulk'] = mu_lead + mu
    
    for j, s in enumerate(Smags):
        params_TI['S_mag'] = s
        
        ABS_spectrum = get_ABS_spectrum_calibrated(systf1, params_TI, phases, True, "500712")
        
        path = path_generator.generate_path(["Data","ABS_Spectrums",f"ABS_{params_TI['mu_lead1']}_flux_{flux}_phases_{len(phases)}"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
        np.savetxt(fname = path, X = ABS_spectrum)
        
        time_current = time.time()
        percentage = (i*len(Smags) + j + 1)/(len(Smags)*len(mus))
        time_left = (time_current - time_start)/percentage - (time_current - time_start) 
        print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        



#%%

############################################################################
###################     Generate current data
############################################################################


mus   =  np.linspace(-0.05, 0.2, 51)
Smags =  np.linspace(0, 0.05, 51)
mu_lead = 0.042    # 0.02, 0.042

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
params_TI["B_x"]      = 0.5/(W_y*H_z) 

params_TI['S_mag'] = Smags[get_closest_index( 0.01 , Smags)] 
params_TI['mu_bulk'] = mu_lead + mus[get_closest_index( 0.05, mus)] 

phases = np.linspace(-0.1, 1.1, 121) # make sure pi is in the region

ABS_spectrum = get_ABS_spectrum_calibrated(systf1, params_TI, phases, True, "500712")

only_ABS = ABS_spectrum[1:]

Majorana = ABS_spectrum[0:1]
Majorana[0, 111 : ] = - Majorana[0, 111 : ]

grad_ABS = - np.gradient(only_ABS, axis = 1) 
grad_Majorana = np.gradient(Majorana, axis = 1)

step = (phases [1] - phases[0]) * np.pi

JC1 = (np.sum(grad_ABS, axis = 0) + np.reshape(grad_Majorana, -1)) / step
JC2 = (np.sum(grad_ABS, axis = 0) - np.reshape(grad_Majorana, -1)) / step

JC1 = JC1[(phases >= 0) & (phases <= 1.001)]
JC2 = JC2[(phases >= 0) & (phases <= 1.001)]

ABSC = np.sum(grad_ABS, axis = 0)[(phases >= 0) & (phases <= 1.001)] / step

MajoranaC = ( - np.reshape(grad_Majorana, -1) / step)[(phases >= 0) & (phases <= 1.001)]


#%%

############################################################################
###################     Plot current data
############################################################################

plot_phases = phases[(phases >= 0) & (phases <= 1.001)]

fig, ax = plt.subplots()
ax.set_xlabel(r'$\phi/\pi$')
ax.set_ylabel(r'd($E/\Delta$)/d\phi') 

ax.plot(plot_phases, ABSC, 'C1', label = "ABS Current")
ax.plot(plot_phases, MajoranaC, 'C2', label = "Majorana Current")
ax.plot(plot_phases, JC1, 'C3', label = "Total 1 Current")
ax.plot(plot_phases, JC2, 'C4', label = "Total 2 Current")
ax.legend()
ax.grid()
s_mag = str(params_TI['S_mag'])
mu = str(params_TI['mu_bulk'])
ax.set_title("Josephson Current  mu bulk: " + mu + ", s mag: " + s_mag)

max_index_JC1 = np.argmax(JC1)
max_index_JC2 = np.argmax(JC2)

plt.axvline(x = plot_phases[max_index_JC1], color = 'C3')
plt.axvline(x = plot_phases[max_index_JC2], color = 'C4')

print("Max value for JC1 : ", JC1[max_index_JC1])
print("Max value for JC2 : ", JC2[max_index_JC2])




#%%
sols = np.concatenate((Majorana, ABS_spectrum), axis = 0)

fig, ax = plt.subplots()
ax.set_xlabel(r'$\phi/\pi$')
ax.set_ylabel(r'$E/\Delta$') 

i = 0
for sol in sols:
    ax.plot(phases, sol, 'C'+str(i), label=str(i))
    ax.plot(phases, -1*sol, 'C'+str(i))
    i+=1
ax.legend()


