from Model import *


#%%
##############################################################################
########             Matric generation and save
##############################################################################


mus   =  np.linspace(-0.05, 0.2, 51)
Smags =  np.linspace(0, 0.05, 51)

mu_lead = 0.042    # 0.02, 0.042

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead

params_TI["B_x"]      = 0.000001/(W_y*H_z)   # half flux is 0.5/(W_y*H_z)
flux = "000001"       # quarter, half, zero

#mus = mus[int(len(mus)*0.71):] 

time_start = time.time()

for i, mu in enumerate(mus):
    params_TI['mu_bulk'] = mu_lead + mu
    
    for j, s in enumerate(Smags):
        params_TI['S_mag'] = s

        smat_e = scattering_matrix(systf1, params_TI, calibration=None)[0]
        
        path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}_reorder"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
        np.savetxt(fname =path, X = smat_e)
        
        time_current = time.time()
        percentage = (i*len(Smags) + j + 1)/(len(Smags)*len(mus))
        time_left = (time_current - time_start)/percentage - (time_current - time_start) 
        print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        

#%%
# Plotting specific points from spectrum that are anomalous

mu_lead = 0.042

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
params_TI["B_x"]      = (0.5)/(W_y*H_z)   # half flux is 0.5/(W_y*H_z)
flux = "half"        # quarter, half, zero
channel = "single" if mu_lead == 0.02 else "multi"
plot = "diag"

params_TI['S_mag'] = 0.002
mismatch = 0.05

params_TI['mu_bulk'] = mu_lead + mismatch


#smat_e = scattering_matrix(systf1, params_TI, calibration=None)[0]

phases = np.linspace(0,2,45)
plot_ABS_spectrum_calibrated(systf1, params_TI, phases = phases)

plt.title(f"{channel}_flux={flux}_{plot}_S_mag={params_TI['S_mag']}_mismatch={mismatch}")

path = path_generator.generate_path("Images",f"SpecificPoints_{channel}_flux={flux}_{plot}_S_mag={params_TI['S_mag']}_mismatch={mismatch}","png")

plt.savefig(path) 


        
#%%

##############################################################################
########             Time reversal check data save
##############################################################################


mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)
all_max_diff = []
mu_lead = 0.042
params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
flux = "500712" # quarter, half
suffix = "_new_disorder"

computing_error = 1e-12

time_start = time.time()

for i, mu in enumerate(mus):
    abs_diff = []
    params_TI['mu_bulk'] = mu_lead + mu
    
    for j, s in enumerate(Smags):
        params_TI['S_mag'] = s
        
        try:
            path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}{suffix}"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
            smat_e = np.loadtxt(fname = path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})
        except:
            print(mu, s)
            raise Exception("ff")
        
        all_diff = np.absolute(smat_e + smat_e.T)
        max_diff = np.amax(all_diff)
        abs_diff.append(max_diff)
        
        time_current = time.time()
        percentage = (i*len(Smags) + j + 1)/(len(Smags)*len(mus))
        time_left = (time_current - time_start)/percentage - (time_current - time_start) 

        print("Time elapsed : ", "{:.1f}".format(time_current - time_start), "Time remaining : ", "{:.1f}".format(time_left), "Percentage : ", "{:.3f}".format(percentage*100))
        
    all_max_diff.append(abs_diff)

path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}{suffix}","txt")
np.savetxt(fname =path, X = np.array(all_max_diff))



#%%
##############################################################################
########             Time reversal check data display
##############################################################################

from mpl_toolkits.axes_grid1 import make_axes_locatable

mu_lead = 0.042
flux = "500001"
suffix = "_negative_disorder" #_negative_disorder, _new_disorder

mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)

path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}{suffix}","txt")
all_max_diff = np.log10(np.loadtxt(path))

plt.rcParams['font.size'] = '10'
fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)


cax = div.append_axes('right', '5%', '5%')
#ax.set_title(f'Log10 of Abs Diff mu_lead = {mu_lead} flux = {flux}')
im = ax.imshow(all_max_diff,extent = (min(Smags), max(Smags), min(mus), max(mus)), 
                 origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)

path = path_generator.generate_path("Images",f"TRS_log_mu_lead_{mu_lead}_flux_{flux}{suffix}","png")

plt.savefig(path, dpi = 1000) 


#%%
##############################################################################
########             Check tables
##############################################################################

def get_closest_index(val, list_data):
    index = 0
    diff = 1000
    for i in range(len(list_data)):
        if abs(val - list_data[i]) < diff:
            diff = abs(val - list_data[i])
            index = i
    return index

def check_probabilities(matrix):
    buf = np.abs(matrix)**2
    df = pd.DataFrame(buf)
    display(df)


mu_lead = 0.042
flux = "000001" # quarter, half, zero

mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)

# 0.155 0, 0.1  0.05

mu_bulk = 0 + mu_lead # 1.19137784e-01 + mu_lead # mus[-5] + mu_lead # mus[get_closest_index( 0.170, mus)] + mu_lead
smag = 0 # Smags[27] #   Smags[get_closest_index( 0.001 , Smags)]
print(mu_bulk - mu_lead, smag)

params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
params_TI["B_x"]      = 0#.500001/(W_y*H_z)   # half flux is 0.5/(W_y*H_z)     
params_TI['mu_bulk']  = mu_bulk
params_TI['S_mag']    = smag

smat_e_new = scattering_matrix(systf1, params_TI, calibration = np.identity(4), reorder=False)[0]
#smat_e_old = scattering_matrix(systf1, params_TI, calibration = np.identity(4))[0]
path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{mu_lead}_flux_{flux}_reorder"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
#smat_e = np.loadtxt(fname = path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})

print("Original")
df = pd.DataFrame(smat_e_new)
display(df)
print("Added to T\n", abs(smat_e_new + smat_e_new.T))
"""
print("Original ABS")
df = pd.DataFrame(np.absolute(smat_e_new))
display(df)
print("angle old")
df = pd.DataFrame(np.angle(smat_e_old)/np.pi)
display(df)
print("angle new")
df = pd.DataFrame(np.angle(smat_e_new)/np.pi)
display(df)

print("Added to T\n", abs(smat_e_new + smat_e_new.T))
print("Saved matrix Added to T\n", abs(smat_e + smat_e.T))
"""
"""
print("Added to T")
df = pd.DataFrame(smat_e_new + smat_e_new.T)
display(df)

print("Subtracted from T")
df = pd.DataFrame(smat_e_new - smat_e_new.T)
#display(df)
print("ABSolutes, Subtracted from T")
df = pd.DataFrame(np.abs(smat_e_new) - np.abs(smat_e_new.T))
#display(df)

#"""

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

mu_lead = 0.042        #0.02 0.042
params_TI['mu_lead1'] = mu_lead
params_TI['mu_lead2'] = mu_lead
flux = "500001" # quarter, half, zero

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
            path = path_generator.generate_path(["Data","Antisymmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}_new_disorder"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
            #path = path_generator.generate_path(["Data","Symmetric_Scattering_Matrices",f"Scattering_Matrices_{params_TI['mu_lead1']}_flux_{flux}"],f"mu_bulk_{params_TI['mu_bulk']}_S_mag_{params_TI['S_mag']}","txt")
            smat_e = np.loadtxt(fname = path,dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})
        except:
            print(mu, s)
            
        
        smat_e_T = np.transpose(smat_e)
        
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
    
    
#######   check if the adress is asym or sym
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

mu_lead = 0.042
flux = "500001"

mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)

path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}","txt")
all_max_diff = np.log(np.loadtxt(path)) / np.log(10)
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}_diagonal","txt")
all_max_diff_diagonal = np.log(np.loadtxt(path)) / np.log(10)
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}_reflection","txt")
all_max_diff_reflection = np.log(np.loadtxt(path)) / np.log(10)
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{params_TI['mu_lead1']}_flux_{flux}_transmission","txt")
all_max_diff_transmission = np.log(np.loadtxt(path)) / np.log(10)  

v_min = min(np.min(all_max_diff), np.min(all_max_diff_diagonal))
v_max = max(np.max(all_max_diff), np.max(all_max_diff_diagonal))

fig, ax = plt.subplots()
ax.set_xlabel('Correlated disorder strength (eV)')
ax.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title(f'Abs Diff mu_lead = {mu_lead} flux = {flux}')
im = ax.imshow(all_max_diff,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = v_min, vmax = v_max ,origin='lower', aspect='auto')
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
im = ax.imshow(all_max_diff_diagonal,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = v_min, vmax = v_max ,origin='lower', aspect='auto')
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
im = ax.imshow(all_max_diff_reflection,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = v_min, vmax = v_max, origin='lower', aspect='auto')
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
im = ax.imshow(all_max_diff_transmission,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = v_min, vmax = v_max ,origin='lower', aspect='auto')
fig.colorbar(im, cax=cax)
path = path_generator.generate_path("Images",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}_transmission","png")
plt.savefig(path) 
plt.show()

#%%

####################################################
#         Present the data for poster
####################################################


from mpl_toolkits.axes_grid1 import make_axes_locatable

mu_lead = 0.042
flux = "500001"
#flux_old = "half"

mus = np.linspace(-0.05, 0.2, 51)
Smags = np.linspace(0, 0.05, 51)

path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}","txt")
all_max_diff_new = np.log(np.loadtxt(path)) / np.log(10)
path = path_generator.generate_path("Data",f"Diff_asym_scattering_mu_lead_{mu_lead}_flux_{flux}_diagonal","txt")
all_max_diff_diagonal_new = np.log(np.loadtxt(path)) / np.log(10)
all_max_diff_new = np.where(all_max_diff_new > all_max_diff_diagonal_new, all_max_diff_new, all_max_diff_diagonal_new)
"""
path = path_generator.generate_path("Data",f"Diff_sym_scattering_mu_lead_{mu_lead}_flux_{flux_old}","txt")
all_max_diff_old = np.log(np.loadtxt(path)) / np.log(10)
path = path_generator.generate_path("Data",f"Diff_sym_scattering_mu_lead_{mu_lead}_flux_{flux_old}_diagonal","txt")
all_max_diff_diagonal_old = np.log(np.loadtxt(path)) / np.log(10)
all_max_diff_old = np.where(all_max_diff_old > all_max_diff_diagonal_old, all_max_diff_old, all_max_diff_diagonal_old)
"""

v_min = min(np.min(all_max_diff_new), np.min(all_max_diff_new))
v_max = max(np.max(all_max_diff_new), np.max(all_max_diff_new))

plt.rcParams['font.size'] = '18'

fig, ax1 = plt.subplots()
fig.set_size_inches(10, 7.5)
ax1.set_xlabel('Correlated disorder strength (eV)')
ax1.set_ylabel('Bulk potential mismatch (eV)')
div = make_axes_locatable(ax1)
cax = div.append_axes('right', '5%', '5%')
#ax1.set_title('Maximum TRS difference for half flux with callibration')
im1 = ax1.imshow(all_max_diff_new ,extent = (min(Smags), max(Smags), min(mus), max(mus)),vmin = v_min, vmax = v_max ,origin='lower', aspect='auto')
fig.colorbar(im1, cax=cax)
path = path_generator.generate_path("Images",f"Poster_comparison","png")
plt.savefig(path, dpi = 1600) 

plt.show()

