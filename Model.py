import kwant
print(kwant.version.version)

import matplotlib

import matplotlib.pyplot as plt
#plt.style.use('genstyle')

import pandas as pd

import scipy.sparse.linalg as sla

import numpy as np
import kwant
import kwant.continuum
import peierls as peierls

import adaptive
import holoviews as hv
from holoviews import opts
#adaptive.notebook_extension()     only for jupyter
from concurrent.futures import ProcessPoolExecutor
from operator import itemgetter

import  sympy
from sympy.physics.matrices import msigma, Matrix
from sympy import eye
from sympy.physics.quantum import TensorProduct

from sympy.utilities.exceptions import SymPyDeprecationWarning
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)


import scipy.signal
from scipy.stats import multivariate_normal


import time
import os
import path_generator as path_generator

##############################




pd.options.display.float_format = "{:,.2f}".format




##############################
#      SYSTEM CREATION



# shape of system
def get_shape(L, W, H):
    def shape(site):
        (x, y, z) = site.pos
        return (0 <= x <= L and 0 <= y <= W and 0 <= z <= H)
    return shape, np.array([0, 0, 0])

# shape of leads
def get_shape_lead_0(L, W, H):
    def shape(site):
        (x, y, z) = site.pos
        return (0 <= y <= W and 0 <= z <= H)
    return shape, np.array([0, 0, 0])

def get_shape_lead_1(L, W, H):
    def shape(site):
        (x, y, z) = site.pos
        return (0 <= y <= W and 0 <= z <= H)
    return shape, np.array([L, 0, 0])

L_x = 100
W_y = 100
H_z = 100




############################
#       3D HAMILTONIAN



import kwant.continuum


norbs = 8 # 8 orbitals (2 for particle-hole x 2 for spin up-down x 2 for orbitals A/B)
ham_TI = ("- mu_bulk * kron(sigma_z, sigma_0, sigma_0) + "
          "{epsilon} * kron(sigma_z, sigma_0, sigma_0) + "
          "{M} * kron(sigma_z, sigma_0, sigma_z) - "
          "A_perp * k_x * kron(sigma_z, sigma_y, sigma_x) + "
          "A_perp * k_y * kron(sigma_0, sigma_x, sigma_x) + "
          "A_z * k_z * kron(sigma_z, sigma_0, sigma_y) + "
          "m_z * kron(sigma_z, sigma_z, sigma_0) + "
          "S_mag * S_imp(site) * kron(sigma_z, sigma_0, sigma_0)")

ham_TI_lead1 = ("- mu_lead1 * kron(sigma_z, sigma_0, sigma_0) + "
          "{epsilon} * kron(sigma_z, sigma_0, sigma_0) + "
          "{M} * kron(sigma_z, sigma_0, sigma_z) - "
          "A_perp * k_x * kron(sigma_z, sigma_y, sigma_x) + "
          "A_perp * k_y * kron(sigma_0, sigma_x, sigma_x) + "
          "A_z * k_z * kron(sigma_z, sigma_0, sigma_y) + "
          "m_z * kron(sigma_z, sigma_z, sigma_0) ")

ham_TI_lead2 = ("- mu_lead2 * kron(sigma_z, sigma_0, sigma_0) + "
          "{epsilon} * kron(sigma_z, sigma_0, sigma_0) + "
          "{M} * kron(sigma_z, sigma_0, sigma_z) - "
          "A_perp * k_x * kron(sigma_z, sigma_y, sigma_x) + "
          "A_perp * k_y * kron(sigma_0, sigma_x, sigma_x) + "
          "A_z * k_z * kron(sigma_z, sigma_0, sigma_y) + "
          "m_z * kron(sigma_z, sigma_z, sigma_0) ")


epsilon = "(C_0 - C_perp * (k_x**2 + k_y**2) - C_z * k_z**2)"
M = "(M_0 - M_perp * (k_x**2 + k_y**2) - M_z * k_z**2)"

ham_TI = ham_TI.format(epsilon=epsilon, M=M, C_0="C_0")
ham_discr, coords = kwant.continuum.discretize_symbolic(ham_TI)


ham_TI_lead1 = ham_TI_lead1.format(epsilon=epsilon, M=M, C_0="C_0")
ham_discr_lead1, coords = kwant.continuum.discretize_symbolic(ham_TI_lead1)

ham_TI_lead2 = ham_TI_lead2.format(epsilon=epsilon, M=M, C_0="C_0")
ham_discr_lead2, coords = kwant.continuum.discretize_symbolic(ham_TI_lead2)


signs = [-1, -1, -1, -1, 1, 1, 1, 1]
vector_potential='[-B_z * y, -B_x * z, -B_y * x]'
ham_discr = peierls.apply(ham_discr, coords, A=vector_potential, signs=signs)
ham_discr_lead1 = peierls.apply(ham_discr_lead1, coords, A=vector_potential, signs=signs)
ham_discr_lead2 = peierls.apply(ham_discr_lead2, coords, A=vector_potential, signs=signs)





#################################




a = 10
ti_syst = kwant.continuum.build_discretized(ham_discr, coords, grid=a)
ti_lead_0 = kwant.continuum.build_discretized(ham_discr_lead1, coords, grid=a)
ti_lead_1 = kwant.continuum.build_discretized(ham_discr_lead2, coords, grid=a)

###########


syst1 = kwant.Builder()
_ = syst1.fill(ti_syst, *get_shape(L_x, W_y, H_z))

lat = kwant.lattice.cubic(a, norbs=norbs)


###########



sigma_0 = np.identity(2)
sigma_x = np.array([[0, 1], [1, 0]])
sigma_z = np.array([[1, 0], [0, -1]])
conservation_law = -np.kron(sigma_z, np.kron(sigma_0, sigma_0))
particle_hole = np.kron(sigma_x, np.kron(sigma_0, sigma_0))
sigma_TR = np.array([[0, -1], [1, 0]])
# i*sigma_y
time_reversal = np.kron(sigma_z, np.kron(sigma_TR, sigma_0))

nleads = 2

lead_0 = kwant.Builder(kwant.TranslationalSymmetry((-a, 0, 0)),
                       conservation_law=conservation_law,
                       particle_hole=particle_hole
                      )
lead_1 = kwant.Builder(kwant.TranslationalSymmetry((a, 0, 0)),
                       conservation_law=conservation_law,
                       particle_hole=particle_hole
                      )

lead_0.fill(ti_lead_0, *get_shape_lead_0(L_x, W_y, H_z))
lead_1.fill(ti_lead_1, *get_shape_lead_1(L_x, W_y, H_z))

syst1.attach_lead(lead_0)
syst1.attach_lead(lead_1)


##########



systf1 = syst1.finalized()



########################################
#       DEFINE DISORDER


# Define a random number (following gaussian distribution) table for disorder
disorder_3D = np.array([[[kwant.digest.gauss('('+str(ind_x)+ ',' +str(ind_y) + ',' +str(ind_z)+')')
                       for ind_z in range(H_z//a+1)]
                      for ind_x in range(L_x//a+1)]
                     for ind_y in range(W_y//a+1)])

disorder_x = np.array([[[kwant.digest.gauss(str(ind_x))
                       for ind_z in range(H_z//a+1)]
                      for ind_x in range(L_x//a+1)]
                     for ind_y in range(W_y//a+1)])

disorder_y = np.array([[[kwant.digest.gauss(str(ind_y))
                       for ind_z in range(H_z//a+1)]
                      for ind_x in range(L_x//a+1)]
                     for ind_y in range(W_y//a+1)])

disorder_z = np.array([[[kwant.digest.gauss(str(ind_z))
                       for ind_z in range(H_z//a+1)]
                      for ind_x in range(L_x//a+1)]
                     for ind_y in range(W_y//a+1)])

def get_S_imp_3D():
    def S_imp(site):
        ind_x = int(site.pos[0]/a)
        ind_y = int(site.pos[1]/a)
        ind_z = int(site.pos[2]/a)
        return disorder_3D[ind_y, ind_x, ind_z]
    return S_imp

def get_S_imp_x():
    def S_imp(site):
        ind_x = int(site.pos[0]/a)
        ind_y = int(site.pos[1]/a)
        ind_z = int(site.pos[2]/a)
        return disorder_x[ind_y, ind_x, ind_z]
    return S_imp

def get_S_imp_y():
    def S_imp(site):
        ind_x = int(site.pos[0]/a)
        ind_y = int(site.pos[1]/a)
        ind_z = int(site.pos[2]/a)
        return disorder_y[ind_y, ind_x, ind_z]
    return S_imp

def get_S_imp_z():
    def S_imp(site):
        ind_x = int(site.pos[0]/a)
        ind_y = int(site.pos[1]/a)
        ind_z = int(site.pos[2]/a)
        return disorder_z[ind_y, ind_x, ind_z]
    return S_imp

################################


correlation_length = 10# in units of a (lattice constant)
#Smag_imp = 0.05

x, y, z = np.mgrid[-3*correlation_length:3*correlation_length+1:1, -3*correlation_length:3*correlation_length+1:1, -3*correlation_length:3*correlation_length+1:1]
pos = np.stack((x, y, z), axis=3)
rv = multivariate_normal([0, 0, 0], [[correlation_length**2, 0, 0], [0, correlation_length**2, 0], [0, 0, correlation_length**2]])
filter_kernel = rv.pdf(pos)

disorder = np.array([[[kwant.digest.gauss('('+str(ind_x)+ ',' +str(ind_y) + ',' +str(ind_z)+')')
                       for ind_z in range(H_z//a+1)]
                      for ind_x in range(L_x//a+1)]
                     for ind_y in range(W_y//a+1)])

disorder_fil = scipy.signal.fftconvolve(disorder, filter_kernel, mode='same')
disorder_fil = disorder_fil/np.std(disorder_fil.flatten())

def get_S_imp_correlated():
    def S_imp(site):
        ind_x = int(site.pos[0]/a)
        ind_y = int(site.pos[1]/a)
        ind_z = int(site.pos[2]/a)
        #return disorder[ind_y, ind_x, ind_z]
        return disorder_fil[ind_y, ind_x, ind_z]
    return S_imp




#####################################
#      PARAMS DICT



#params_toy = dict(C_0=0.0,
#                  C_2=3.0,
#                  mu=0.012
#                 )


params_toy = dict(C_0=0.0,
                  C_2=3.0,
                  mu=0.0205,
                  mu_lead1=0.0205,
                  mu_lead2=0.021,
                  S_imp = get_S_imp_x(),
                  S_mag = 0.000
                 )

params_toy_sim = dict(C_0=0.0,
                  C_2=3.0,
                  mu=0.0,
                 )

#params_TI = dict(A_perp=3.0,
#                 A_z=3.0,
#                 M_0=0.3,
#                 M_perp=15.0,
#                 M_z=15.0,
#                 C_0=0.0,
#                 C_perp=0.0,
#                 C_z=0.0,
#                 m_z=0.0,
#                 mu=0.1
#                )

params_TI = dict(A_perp=3.0,
                 A_z=3.0,
                 M_0=0.3,
                 M_perp=15.0,
                 M_z=15.0,
                 C_0=0.0,
                 C_perp=0.0,
                 C_z=0.0,
                 m_z=0.0,
                 mu_bulk=0.042,
                 mu_lead1=0.042,
                 mu_lead2=0.042,
                 S_imp = get_S_imp_correlated(),
                 S_mag = 0.0,
                 B_x = 0.5/(W_y*H_z),
                 B_y = 0,
                 B_z = 0,
                 phi_0=1.0,
                 exp=np.exp,
                 a = 10
                 )

params_TI_sim = dict(A_perp=3.0,
                 A_z=3.0,
                 M_0=0.3,
                 M_perp=15.0,
                 M_z=15.0,
                 C_0=0.0,
                 C_perp=0.0,
                 C_z=0.0,
                 m_z=0.0,
                 mu=0
                )



###################################
#      ANDREEV SPECTRUM, FUNCTIONS


def scattering_matrix(syst, p, calibration=None):
    smat = kwant.smatrix(syst, params=p)
    #print(smat.submatrix(0, 0))
    #print(smat.submatrix(1, 1))
    size_L = smat.submatrix((0, 0), (0, 0)).shape[0]
    size_R = smat.submatrix((1, 0), (1, 0)).shape[0]
    mask_e = np.zeros((2*size_L+2*size_R, 2*size_L+2*size_R), dtype=bool)
    mask_h = np.zeros((2*size_L+2*size_R, 2*size_L+2*size_R), dtype=bool)
    

    mask_e[:size_L, :size_L] = True
    mask_h[size_L:2*size_L, size_L:2*size_L] = True
    mask_e[:size_L, 2*size_L:2*size_L+size_R] = True
    mask_e[2*size_L:2*size_L+size_R,:size_L] = True
    mask_h[size_L:2*size_L, 2*size_L+size_R:2*size_L+2*size_R] = True
    mask_h[2*size_L+size_R:2*size_L+2*size_R,size_L:2*size_L] = True
    mask_e[2*size_L:2*size_L+size_R,2*size_L:2*size_L+size_R] = True
    mask_h[2*size_L+size_R:2*size_L+2*size_R,2*size_L+size_R:2*size_L+2*size_R] = True
    

    smat_e = smat.data[mask_e].reshape((size_L+size_R,size_L+size_R))
    smat_h = smat.data[mask_h].reshape((size_L+size_R,size_L+size_R))
    #smat_h = reverse(smat_h)

    
    smat_h_copy = np.array(list(smat_h))
    
    i = 0
    for a in smat_h:
        j = 0
        for b in smat_h[0]:
            if i < size_L:
                new_i = size_L - 1 - i
            else:
                new_i = 2*size_L + size_R - 1 - i
            if j < size_L: 
                new_j = size_L - 1 - j
            else:
                new_j = 2*size_L + size_R - 1 - j
            smat_h[new_i][new_j] = smat_h_copy[i][j]
            j+=1
        i+=1

    if calibration is None:
        calibration_e = np.identity(size_L+size_R, dtype=complex)

    ## Assuming the left lead and the right lead have the same dimensions (same number of propagating modes)
    ## Calibrate the phase shift        
        for i in range(0, size_L):
            shift_e = np.angle(smat_e[i][i+size_L]) - np.angle(smat_e[i+size_L][i]) - np.pi
            calibration_e[i+size_L][i+size_L] = np.exp(1j*shift_e/2)
    else:
        print('same calibration!')
        calibration_e = calibration
    
    smat_e = calibration_e@smat_e@calibration_e.conj()
    return [smat_e, size_L, size_R]



def energies_over_delta_calibrated(smat_e, size_L, size_R, phases=np.zeros(nleads)):
## Add in the user-defined phase difference between the two leads
    mat_phase_1 = np.identity(size_L)*np.exp(1j*phases[0])
    mat_phase_2 = np.identity(size_R)*np.exp(1j*phases[1])
    mat_phase = np.block([
        [mat_phase_1, np.zeros((size_L, size_R))],
        [np.zeros((size_R, size_L)), mat_phase_2]
    ])
    #    Use particle-hole symmetry    
    smat_prod = smat_e.T.conj() @ mat_phase @ smat_e.T @ mat_phase.conj()
        #print(smat_prod)
    
    operator = 0.5 * np.eye(smat_prod.shape[0]) + 0.25 * (smat_prod + smat_prod.T.conj())
    return np.sqrt(np.linalg.eigvalsh(operator))


def plot_ABS_spectrum_calibrated(syst, p, phases, r=False, calibration=None):

    params = p.copy()
    fig, ax = plt.subplots()
    ax.set_xlabel(r'$\phi/\pi$')
    ax.set_ylabel(r'$E/\Delta$')

    sol_list = []
    smat_e, size_L, size_R = scattering_matrix(syst, p, calibration)
    for p in (phases)*np.pi:
        phase = [0, p]
        sol_list.append(energies_over_delta_calibrated(smat_e, size_L, size_R, phases=phase))
    
    sol_list = np.array(sol_list).T
    sol_list2 = sol_list[::2]
    i = 0
    for sol in sol_list2:
        ax.plot(phases, sol, 'C'+str(i), label=str(i))
        ax.plot(phases, -1*sol, 'C'+str(i))
        i+=1
    ax.legend()
    if r:
        return sol_list2
    
def find_gap(syst, p, n=0, calibration=None):
    smat_e, size_L, size_R = scattering_matrix(syst, p, calibration)
    phases = [0, np.pi]
    return 2*energies_over_delta_calibrated(smat_e, size_L, size_R, phases=phases)[n]

def find_SB(syst, p, n=0, calibration=None):
    smat_e, size_L, size_R = scattering_matrix(syst, p, calibration)
    SB_matrix = np.abs(smat_e - smat_e.T)
    return np.max(SB_matrix)





