from Model import *

#%%

# No disorder find the offset values

mus =  np.linspace(-0.3, 0.2, 10001)
Smag = 0 
mu_lead = 0.042    # 0.02, 0.042

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead

params_TI["B_x"]      = 0.500712/(W_y*H_z)   # half flux is 0.5/(W_y*H_z)
flux = "500712"       # quarter, half, zero

#mus = mus[int(len(mus)*0.71):] 

time_start = time.time()

for i, mu in enumerate(mus):
    params_TI['mu_bulk'] = mu_lead + mu
    params_TI['S_mag'] = Smag

    smat_e = scattering_matrix(systf1, params_TI, calibration=None)[0]
    
    path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}_no_disorder"],f"mu_bulk_{params_TI['mu_bulk']}","txt")
    np.savetxt(fname =path, X = smat_e)
    
    time_current = time.time()
    percentage = (i + 1)/(len(mus))
    time_left = (time_current - time_start)/percentage - (time_current - time_start) 
    print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))

#%%

mus =  np.linspace(-0.3, 0.2, 101)
mu_lead = 0.042 
params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
flux = "500712"
suffix = "_L200"  
params_TI["B_x"]      = 0.500712/(W_y*H_z)
params_TI['S_mag'] = 0.0

r_i = get_closest_index( -0.042, mus)  # reflection index

max_value_1 = []

for i, mu in enumerate(mus):
    params_TI['mu_bulk'] = mu_lead + mu
    path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}_L200"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
    smat_e = np.loadtxt(fname =path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})
    max_value_1.append(np.amax(np.absolute(smat_e + smat_e.T)))


r_mus = []

for i in range(len(mus)):
    if mus[i] >= -0.042:
        r_mus.append(-0.042 - abs(mus[i] -( -0.042)))  
    else:
        r_mus.append(-0.042 + abs(mus[i] - (-0.042 )))
        

plt.rcParams['font.size'] = '12'
fig, ax = plt.subplots()
ax.set_yscale("log")
#plt.grid()
plt.rcParams['font.size'] = '14'
#ax.title("Plot for TRS multichannel mu lead = 0.042 eV")
ax.set_xlabel("Bulk potential mismatch (eV)")
ax.set_ylabel("Magnitude of TRS conservation")
ax.plot(mus, max_value_1, label = "Max values of TRS")
ax.plot(r_mus, max_value_1, label = "Max values of TRS")
plt.show()

    
#%%
mu_lead = 0.042
mus_1 =  np.linspace(-0.05, 0.2, 10001)
max_value_1 = []
flux = "half"  
for i, mu in enumerate(mus_1):
    path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{mu_lead}_flux_{flux}_no_disorder"],f"mu_bulk_{mu_lead + mu}","txt")
    smat_e = np.loadtxt(fname =path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})
    max_value_1.append(np.amax(np.absolute(smat_e + smat_e.T)))
    
mus_2 = np.linspace(-0.05, 0.2, 51)
max_value_2 = []
flux = "500001"  
for i, mu in enumerate(mus_2):
    path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{mu_lead}_flux_{flux}_reorder"],f"mu_bulk_{mu_lead + mu}_S_mag_{0}","txt")
    #smat_e = np.loadtxt(fname =path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})
    #max_value_2.append(np.amax(np.absolute(smat_e + smat_e.T)))
    
plt.rcParams['font.size'] = '12'
fig, ax = plt.subplots()
ax.set_yscale("log")
plt.grid()
plt.rcParams['font.size'] = '14'
#ax.title("Plot for TRS multichannel mu lead = 0.042 eV")
ax.set_xlabel("Bulk potential mismatch (eV)")
ax.set_ylabel("Magnitude of TRS conservation")
ax.plot(mus_1, max_value_1, label = "Max values of TRS")
#plt.plot(mus_2, max_value_2, label = "A few points")
#plt.legend()

path = path_generator.generate_path("Images","Magnitude TRS conservation for no disorder at mu lead = 0.042","png")

plt.savefig(path, dpi = 1000) 

plt.show()