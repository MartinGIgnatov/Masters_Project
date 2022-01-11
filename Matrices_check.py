from Model import *


#%%
##############################################################################
########             Matric generation and save
##############################################################################


mus =  np.linspace(-0.05, 0.2, 51)
Smags =  np.linspace(0, 0.05, 51)

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

        smat_e = scattering_matrix(systf1, params_TI, calibration=None)[0]
        
        path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
        np.savetxt(fname =path, X = smat_e)
        
        time_current = time.time()
        percentage = (i*len(Smags) + j + 1)/(len(Smags)*len(mus))
        time_left = (time_current - time_start)/percentage - (time_current - time_start) 
        print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        
        
#%%

##############################################################################
########             Time reversal check data save
##############################################################################


mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)
all_max_diff = []
params_TI['mu_lead1'] = 0.02
params_TI['mu_lead2'] = 0.02
flux = "half" # quarter, half

computing_error = 1e-12

time_start = time.time()

for i, mu in enumerate(mus):
    abs_diff = []
    params_TI['mu_bulk'] = 0.02 + mu
    
    for j, s in enumerate(Smags):
        params_TI['S_mag'] = s
        
        try:
            path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
            smat_e = np.loadtxt(fname = path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})
        except:
            print(mu, s)
            raise Exception("ff")
            
        
        smat_e_T = np.transpose(smat_e)
        
        
        #"""
        abs_smat_e = np.absolute(smat_e)
        abs_smat_e_T = np.transpose(abs_smat_e)
        Average = (abs_smat_e + abs_smat_e_T)/2
        
        zeros = np.zeros(Average.shape)
        norm_diff = np.absolute(smat_e - smat_e_T)/Average
        

        mask_non_diag = np.ones(zeros.shape) - np.identity(zeros.shape[0])
        
        diff = np.absolute(smat_e + smat_e_T) * mask_non_diag

        
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
##############################################################################
########             Time reversal check data display
##############################################################################

from mpl_toolkits.axes_grid1 import make_axes_locatable

mu_lead = 0.02
flux = "half"

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

#%%
##############################################################################
########             Check tables
##############################################################################
        
mu_lead = 0.042
flux = "half" # quarter, half

mus = np.linspace(-0.05, 0.2, 51) + 0.042
Smags = np.linspace(0, 0.05, 51)

mu = mus[-1]
smag = Smags[-1]

path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{mu_lead}_flux_{flux}"],f"mu_bulk_{mu}_S_mag_{smag}","txt")
smat_e = np.loadtxt(fname = path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})

df = pd.DataFrame(smat_e)
display(df)


#%%

##############################################################################
########             Time reversal check data save, transmission and reflection
##############################################################################

mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)

all_max_diff              = []
all_max_diff_diagonal    = []
all_max_diff_reflection   = []
all_max_diff_transmission = []

mu_lead = 0.02        #0.02 0.042
params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
flux = "half" # quarter, half, zero

computing_error = 1e-12

time_start = time.time()

for i, mu in enumerate(mus):
    abs_diff              = []
    abs_diff_diagonal     = []
    abs_diff_reflection   = []
    abs_diff_transmission = []
    
    params_TI['mu_bulk'] = mu_lead + mu
    
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
        
        side = zeros.shape[0]
        indeces = np.indices(Average.shape)
        
        non_diag_mask     = np.ones(zeros.shape) - np.identity(zeros.shape[0])
        top_left_mask     = np.ones(zeros.shape)*((indeces[0] < side/2) & (indeces[1] < side/2)) 
        bottom_right_mask = np.ones(zeros.shape)*((indeces[0] >= side/2) & (indeces[1] >= side/2)) 
        reflection_mask   = top_left_mask + bottom_right_mask
        transmission_mask = np.ones(zeros.shape) - reflection_mask
        
        diff                  = np.absolute(smat_e + smat_e_T) * non_diag_mask
        all_diff              = np.where( Average < computing_error, zeros , diff)
        all_diff_diagonal     = np.absolute(smat_e + smat_e_T) * np.identity(zeros.shape[0])
        all_diff_reflection   = all_diff * reflection_mask
        all_diff_transmission = all_diff * transmission_mask
        
        max_diff              = np.amax(all_diff)
        max_diff_diagonal     = np.amax(all_diff_diagonal)
        max_diff_reflection   = np.amax(all_diff_reflection)
        max_diff_transmission = np.amax(all_diff_transmission)
        
        abs_diff.append(max_diff)
        abs_diff_diagonal.append(max_diff_diagonal)
        abs_diff_reflection.append(max_diff_reflection)
        abs_diff_transmission.append(max_diff_transmission)
        
        time_current = time.time()
        percentage = (i*len(Smags) + j + 1)/(len(Smags)*len(mus))
        time_left = (time_current - time_start)/percentage - (time_current - time_start) 

        print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        
        
    all_max_diff.append(abs_diff)
    all_max_diff_diagonal.append(abs_diff_diagonal)
    all_max_diff_reflection.append(abs_diff_reflection)
    all_max_diff_transmission.append(abs_diff_transmission)
    
    

path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}","txt")
np.savetxt(fname =path, X = np.array(all_max_diff))
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}_diagonal","txt")
np.savetxt(fname =path, X = np.array(all_max_diff_diagonal))
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}_reflection","txt")
np.savetxt(fname =path, X = np.array(all_max_diff_reflection))
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}_transmission","txt")
np.savetxt(fname =path, X = np.array(all_max_diff_transmission))

#%%
##############################################################################
########             Time reversal check data display
##############################################################################

from mpl_toolkits.axes_grid1 import make_axes_locatable

mu_lead = 0.02
flux = "half"

mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)

path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}","txt")
all_max_diff = np.loadtxt(path)
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}_diagonal","txt")
all_max_diff_diagonal = np.loadtxt(path)
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}_reflection","txt")
all_max_diff_reflection = np.loadtxt(path)
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}_transmission","txt")
all_max_diff_transmission = np.loadtxt(path)

fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title(f'Abs Diff mu_lead = {mu_lead} flux = {flux}')
im = ax.imshow(all_max_diff,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = 0, vmax = np.max(all_max_diff) ,origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)
path = path_generator.generate_path("Images",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}","png")
plt.savefig(path) 
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title(f'Abs Diff mu_lead = {mu_lead} flux = {flux} diagonal')
im = ax.imshow(all_max_diff_diagonal,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = 0, vmax = np.max(all_max_diff_diagonal) ,origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)
path = path_generator.generate_path("Images",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}_diagonal","png")
plt.savefig(path) 
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title(f'Abs Diff mu_lead = {mu_lead} flux = {flux} reflection non diagonal')
im = ax.imshow(all_max_diff_reflection,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = 0, vmax = np.max(all_max_diff),origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)
path = path_generator.generate_path("Images",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}_reflection","png")
plt.savefig(path) 
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title(f'Abs Diff mu_lead = {mu_lead} flux = {flux} transmission')
im = ax.imshow(all_max_diff_transmission,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = 0, vmax = np.max(all_max_diff),origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)
path = path_generator.generate_path("Images",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}_transmission","png")
plt.savefig(path) 
plt.show()


