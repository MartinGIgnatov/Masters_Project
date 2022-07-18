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

params_TI['S_mag'] = Smags[get_closest_index( 0 , Smags)] 
params_TI['mu_bulk'] = mu_lead + mus[get_closest_index( 0.01, mus)] 

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

########################################################
#######    Generate color map plots
########################################################

mus   =  np.linspace(-0.3, 0.2, 101)
mus = np.concatenate((np.linspace(-0.05, -0.3, 51)[::-1], np.linspace(-0.05, 0.2, 51)[1:]))
#mus   =  np.linspace(-0.05, 0.2, 51)

Smags =  np.linspace(0, 0.05, 51)
mu_lead = 0.042    # 0.02, 0.042
flux = "500712" # quarter, half
suffix = ""

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
params_TI["B_x"]      = 0.500712/(W_y*H_z)

phases = np.linspace(-0.1, 1.1, 121) # make sure pi is in the region
step = (phases [1] - phases[0]) * np.pi

time_start = time.time()

Buf1 = []
Buf2 = []
BufM = []

for i, mu in enumerate(mus):
    abs_diff = []
    params_TI['mu_bulk'] = mu_lead + mu
    
    buf1 = []
    buf2 = []
    bufM = []
    
    for j, s in enumerate(Smags):
        params_TI['S_mag'] = s
        
        ABS_spectrum = get_ABS_spectrum_calibrated(systf1, params_TI, phases, True, "500712", False)
        
        only_ABS = ABS_spectrum[1:]
        
        Majorana = ABS_spectrum[0:1]
        Majorana[0, 111 : ] = - Majorana[0, 111 : ]
        
        ABS_C = np.sum(- np.gradient(only_ABS, axis = 1) / step , axis = 0)
        Majorana_C = np.reshape(np.gradient(Majorana, axis = 1) / step , -1)
        
        JC1 = ABS_C + Majorana_C
        JC2 = ABS_C - Majorana_C
        
        max_index_JC1 = np.argmax(JC1)
        max_index_JC2 = np.argmax(JC2)
        
        max_JC1 = JC1[max_index_JC1]
        max_JC2 = JC2[max_index_JC2]
        
        buf1.append(max_JC1)
        buf2.append(max_JC2)
        bufM.append(Majorana_C[np.argmax(Majorana_C)])
        
        time_current = time.time()
        percentage = (i*len(Smags) + j + 1)/(len(Smags)*len(mus))
        time_left = (time_current - time_start)/percentage - (time_current - time_start) 

        print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        
    Buf1.append(buf1)
    Buf2.append(buf2)
    BufM.append(bufM)
        
path = path_generator.generate_path("Data",f"JC1_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
np.savetxt(fname =path, X = np.array(Buf1))
path = path_generator.generate_path("Data",f"JC2_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
np.savetxt(fname =path, X = np.array(Buf2))
path = path_generator.generate_path("Data",f"Maj_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
np.savetxt(fname =path, X = np.array(BufM))

#%%

########################################################
#######    Plot color map plots
########################################################


from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.style.use('styling.mplstyle')

mus   =  np.linspace(-0.3, 0.2, 101)
Smags =  np.linspace(0, 0.05, 51)
mu_lead = 0.042    # 0.02, 0.042
flux = "500712" # quarter, half
suffix = ""
params_TI['mu_bulk'] = mu_lead

path = path_generator.generate_path("Data",f"JC1_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
JC1_max = np.loadtxt(path)
path = path_generator.generate_path("Data",f"JC2_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
JC2_max = np.loadtxt(path)
#path = path_generator.generate_path("Data",f"Maj_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
#Maj_max = np.loadtxt(path)



JC_total_max = np.where(JC1_max > JC2_max, JC1_max, JC2_max)

#plt.plot(mus, JC_total_max[:,0])

v_min = 0.5 # min(np.amin(JC1_max), np.amin(JC2_max))
v_max = 1.3 #max(np.amax(JC1_max), np.amax(JC2_max))



suffix = "_Extended"

fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
#ax.set_title(f'Max Josephson Current mu lead = {mu_lead} flux = {flux}')
im = ax.imshow(JC_total_max,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = v_min, vmax = v_max ,origin='lower', aspect='auto', cmap = 'coolwarm')
cbar = fig.colorbar(im, cax=cax)
cbar.ax.set_ylabel(r'Josephson Current [$2e\Delta/\hbar$]', rotation=270, labelpad=30)
path = path_generator.generate_path("Images",f"JC_Max_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","png")
plt.tight_layout()  
plt.savefig(path, dpi = 600) 
plt.show()




"""
fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title(f'Max Josephson Current mu lead = {mu_lead} flux = {flux}')
im = ax.imshow(JC1_max,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = v_min, vmax = v_max ,origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)
path = path_generator.generate_path("Images",f"JC1_Max_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","png")
plt.savefig(path) 
plt.show()



fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title(f'Max Josephson Current mu lead = {mu_lead} flux = {flux}')
im = ax.imshow(Maj_max,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = v_min, vmax = v_max ,origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)
path = path_generator.generate_path("Images",f"Maj_Max_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","png")
plt.savefig(path) 
plt.show()
"""

#%%
##################################
#     Only mus

mus   =  np.linspace(-0.05, 0.2, 51)
Smags =  np.linspace(0, 0.05, 51)
mu_lead = 0.042    # 0.02, 0.042
flux = "500712" # quarter, half
suffix = ""
params_TI['mu_bulk'] = mu_lead

path = path_generator.generate_path("Data",f"JC1_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
JC1_max = np.loadtxt(path)
path = path_generator.generate_path("Data",f"JC2_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
JC2_max = np.loadtxt(path)

plt.plot(mus, JC1_max[:, 0], label = "Josephson Current 1")
plt.plot(mus, JC2_max[:, 0], label = "Josephson Current 2")
plt.grid()
plt.legend()
plt.ylabel("Josephson Current")
plt.xlabel("Mu mismatch")
plt.show()

#%%

#################################
####   plot the 4pi periodicity

plt.style.use('styling.mplstyle')

mus   =  np.linspace(-0.05, 0.2, 51)
Smags =  np.linspace(0, 0.05, 51)
mu_lead = 0.042    # 0.02, 0.042

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
params_TI["B_x"]      = 0.5/(W_y*H_z) 

params_TI['S_mag'] = Smags[get_closest_index( 0.02 , Smags)] 
params_TI['mu_bulk'] = mu_lead + mus[get_closest_index( 0.03, mus)] 

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

#ABSC = np.sum(grad_ABS, axis = 0)[(phases >= 0) & (phases <= 1.001)] / step

#MajoranaC = ( - np.reshape(grad_Majorana, -1) / step)[(phases >= 0) & (phases <= 1.001)]

plot_phases = np.linspace(0, 4, 401)

fig, ax = plt.subplots()
ax.set_xlabel(r'$\phi/\pi$')
ax.set_ylabel(r'd($E/\Delta$)/d$\phi$') 
ax.set_xlim(0,4)

#total_JC = np.concatenate((  JC1,  -JC1[-2::-1],  JC2[1:],  -JC2[-2::-1]))
total_JC = np.concatenate((  JC2,  -JC1[-2::-1],  JC1[1:],  -JC2[-2::-1]))
ax.plot(plot_phases, total_JC, 'C1', label = "Total 1 Current")

ax.legend()
#ax.grid()
s_mag = str(params_TI['S_mag'])
mu = str(params_TI['mu_bulk'])
ax.set_title("Josephson Current  mu bulk: " + mu + ", s mag: " + s_mag)
plt.tight_layout()  




