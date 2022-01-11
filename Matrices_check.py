from Model import *


#%%
mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)

mu_lead = 0.02

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
params_TI["B_x"]      = 0   # half flux is 0.5/(W_y*H_z)
flux = "zero"       # quarter, half, zero


time_start = time.time()

for i, mu in enumerate(mus):
    params_TI['mu_bulk'] = mu_lead + mu
    
    for j, s in enumerate(Smags):
        params_TI['S_mag'] = s
        result = scattering_matrix(systf1, params_TI, calibration=None)
        
        smat_e = np.array(pd.DataFrame(result[0]))
        
        path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
        np.savetxt(fname =path, X = smat_e)
        
        time_current = time.time()
        percentage = (i*len(Smags) + j + 1)/(len(Smags)*len(mus))
        time_left = (time_current - time_start)/percentage - (time_current - time_start) 
        print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        
        
#%%
## Generates the data necessary for the plot

mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)
all_max_diff = []
params_TI['mu_lead1'] = 0.042
params_TI['mu_lead2'] = 0.042
flux = "quarter" # quarter, half

computing_error = 1e-12

time_start = time.time()

for i, mu in enumerate(mus):
    abs_diff = []
    params_TI['mu_bulk'] = 0.042+mu
    
    for j, s in enumerate(Smags):
        params_TI['S_mag'] = s
        
        try:
            path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
            smat_e = np.loadtxt(fname = path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})
        
        except:
            print(mu, s)
            
        
        smat_e_T = np.transpose(smat_e)
        
        
        #"""
        abs_smat_e = np.absolute(smat_e)
        abs_smat_e_T = np.transpose(abs_smat_e)
        Average = (abs_smat_e + abs_smat_e_T)/2
        
        zeros = np.zeros(Average.shape)
        norm_diff = np.absolute(smat_e - smat_e_T)/Average
        
        diff = np.absolute(smat_e + smat_e_T)
        
        all_diff = np.where( Average < computing_error, zeros , diff)
        #"""
        
        max_diff = np.amax(all_diff)
        
        abs_diff.append(max_diff)
        
        time_current = time.time()
        percentage = (i*len(Smags) + j + 1)/(len(Smags)*len(mus))
        time_left = (time_current - time_start)/percentage - (time_current - time_start) 

        print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        
        
    all_max_diff.append(abs_diff)
    
    
    
    
"""
df = pd.DataFrame(all_diff)
display(df)
"""

path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}","txt")
np.savetxt(fname =path, X = np.array(all_max_diff))



#%%



from mpl_toolkits.axes_grid1 import make_axes_locatable

mu_lead = 0.042
flux = "quarter"

mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)

path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}","txt")
all_max_diff = np.loadtxt(path)

fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title(f'Abs Diff mu_lead = {mu_lead} flux = {flux}')
im = ax.imshow(all_max_diff,extent = (min(Smags), max(Smags), min(mus), max(mus)), 
                 origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)

path = path_generator.generate_path("Images",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}","png")

plt.savefig(path) 
        