import os
import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
from liblibra_core import *
from libra_py import units as units
from libra_py import influence_spectrum as infsp
from libra_py import data_visualize
from libra_py import data_conv
from libra_py import data_read
from libra_py import data_stat
import libra_py.workflows.nbra.mapping as mapping
import libra_py.workflows.nbra.step4 as step4
import libra_py.workflows.nbra.step3 as step3
import libra_py.workflows.nbra.step2_many_body as step2_many_body
import libra_py.workflows.nbra.step2_analysis as step2_analysis
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

#######################################################################################
# Helper functions
def make_T_matricies_energy_ordered( reindex_nsteps, coeffs, keep_identity_order=False ):
    """
    This function ordered the rows (Slater determinant (SD) states) of the T matrix 
    T matrix = matrix to transform between SD and Many-body bases
    
    Inputs:
        reindex_nsteps (list of lists of ints) - energy ordering of the SD states at each step 
                                                 in terms of their original index
        coeffs (list of lists) - mapping of the SD states to the MB states
        Ex)

            Consider we have a SD basis:
 
            phi 0 = ground state
            phi 1 = first excited state through promotion of alpha elctron

            # Psi 1  phi 0
            [ [1.0], [0] ],

            # Psi 2  phi 1, -phi 1
            [ [norm, -norm], [1, 1] ],

            Where norm = 1.0/sqrt(2)

    Returns
        T matrix for each step, where rows are the SD (energy ordered) and the columns are the MB states
    """

    nsteps = len(reindex_nsteps)
    nSDs   = len(reindex_nsteps[0])
    nMBs   = len(coeffs)

    T_identity_ordered = CMATRIX( nSDs, nMBs )
    for i in range( nMBs ):
        count = 0
        for j in coeffs[i][1]:
            T_identity_ordered.set(j,i,coeffs[i][0][count])
            count += 1

    T_energy_ordered = []
    for step in range( nsteps ):
        T_energy_ordered.append( CMATRIX( nSDs, nMBs ) )
        if len(reindex_nsteps[step]) != nSDs:
            print("\n len(reindex_nsteps[step]) != nstates_sd")
            print("Exiting now")
            sys.exit(0)
        else:
            for i in range( nSDs ):
                for j in range( nMBs ):
                    # To sort by energy
                    T_energy_ordered[step].set( i, j, T_identity_ordered.get( int( reindex_nsteps[step][i] ), j ) )

    return T_energy_ordered


def compute_overlaps_in_parallel( step ):
    """
    Function used for the making the computation of overlaps parallel via python multiprocessing 
    """
    s_sd  = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step], S_ks[0][step],  use_minimal=False )
    s_sd  = data_conv.MATRIX2nparray(s_sd)
    return s_sd

def compute_toverlaps_in_parallel( step ):
    """
    Function used for the making the computation of time-overlaps parallel via python multiprocessing 
    """
    st_sd = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step+1], St_ks[0][step], use_minimal=False )
    st_sd = data_conv.MATRIX2nparray(st_sd)
    return st_sd



#######################################################################################
"""
 1. Read the files that have the energies, overlaps, and time-overlap matricies in the 
    Kohn-Sham basis: E_ks, S_ks, and St_ks.
"""
#######################################################################################

os.system("rm -rf res_mb_sp")
res_dir = "res_mb_sp"
os.mkdir(res_dir)

path = os.getcwd()
res_dir_mb   = path+"/../../../step2/res/"
data_dim     = 62 # rows in E_ks
active_space = range(0,int(data_dim/2)) # alpha channel only here #range(data_dim)
istep = 150 # initial step
fstep = 154 # final step + 1  
dt    = 1.0*units.fs2au

params = { "data_set_paths" : [res_dir_mb], "data_dim":data_dim, "active_space":active_space, "isnap":istep,  "fsnap":fstep }
# Fetching E_ks
params.update({ "data_re_prefix" : "E_ks_",  "data_re_suffix" : "_re", "data_im_prefix" : "E_ks_",  "data_im_suffix" : "_im"  } )
E_ks = data_read.get_data_sets(params)
E_ks[0][-1].show_matrix()

# Fetching S_ks
params.update({ "data_re_prefix" : "S_ks_", "data_re_suffix" : "_re", "data_im_prefix" : "S_ks_", "data_im_suffix" : "_im"  } )
S_ks = data_read.get_data_sets(params)

# Fetching St_ks
params.update({ "data_re_prefix" : "St_ks_", "data_re_suffix" : "_re", "data_im_prefix" : "St_ks_", "data_im_suffix" : "_im"  } )
St_ks = data_read.get_data_sets(params)
St_ks[0][-2].show_matrix()
#sys.exit(0)

step3.apply_orthonormalization_general( S_ks[0], St_ks[0] )
step3.apply_phase_correction_general( St_ks[0] )
#sys.exit(0)

ks_homo_index = 28
min_band_ks   = 28 
max_band_ks   = 29
ks_orbital_indicies = list( range( min_band_ks, max_band_ks + 1 ) )
sd_basis_states_unique = [
                            # Phi 0 
                            [1, -1],
                            # Phi 1 
                            [2, -1],
                            # Phi 2 
                            [-2, 1],
                                                     ]    

# 2.3. Reindex the single-particle excitations into a format expected by Libra
#sd_states_reindexed = step2_many_body.reindex_cp2k_sd_states( ks_homo_index, ks_orbital_indicies, sd_basis_states_unique, sd_format=2 )

# 2.4. Order / sort the single-particle excitations at each timestep by energy or identity
sd_states_reindexed = sd_basis_states_unique
E_sd, sd_states_unique_sorted, sd_states_reindexed_sorted, reindex_nsteps = step3.sort_unique_SD_basis( E_ks[0], sd_basis_states_unique, sd_states_reindexed, istep, fstep, sorting_type="energy" )
#print( reindex_nsteps )
#print(sd_states_reindexed_sorted)
#sys.exit(0)


#######################################################################################
"""
 3. Compute overlaps and time-overlaps between the single-particle excitations 
    3.1. Compute overlaps in the basis of single-particle excitations
    3.2. Compute time-overlaps in the basis of single-particle excitations
    3.3. Convert the data types to Libra's CMATRIX
"""
#######################################################################################

# 3.1. Overlaps
pool = mp.Pool( 24 )
tmp_s_sd = pool.map( compute_overlaps_in_parallel, list( range( fstep - istep ) ) )
pool.close()
pool.join()

# 3.2. Time-overlaps
pool = mp.Pool( 24 )
tmp_st_sd = pool.map( compute_toverlaps_in_parallel, list( range( fstep - istep - 1 ) ) )
pool.close()
pool.join()

# 3.3. Convert the data types to Libra's CMATRIX
S_sd, St_sd = [], []
for step in range( fstep - istep ):
    S_sd.append(  data_conv.nparray2CMATRIX( tmp_s_sd[step]  ) )
for step in range( fstep - istep - 1 ):
    St_sd.append( data_conv.nparray2CMATRIX( tmp_st_sd[step] ) )
#sys.exit(0)



#############################################
"""
 3. Make T matrix and convert to SAC basis
"""

# Define a transformaiton matrix to a SAC basis
norm = 1.0/math.sqrt(2)
coeffs = [
            # Psi 1  phi 0
            [ [1.0], [0] ],
            # Psi 2  phi 1, -phi1
            [ [norm, -norm], [1, 2] ],
         ]

T_energy_ordered = make_T_matricies_energy_ordered( reindex_nsteps, coeffs )  
T_energy_ordered[0].show_matrix()
#sys.exit(0)

E_sac_step, E_sac  = [], []
St_sac = []
for step in range( fstep - istep ):
    E_sac_step.append( T_energy_ordered[step].H() * E_sd[step]  * T_energy_ordered[step] )
for step in range( fstep - istep - 1 ):
    E_sac.append( 0.5 * ( E_sac_step[step] + E_sac_step[step+1] ) )
    St_sac.append( T_energy_ordered[step].H() * St_sd[step] * T_energy_ordered[step+1] )
print("E_sac[0]")
E_sac[0].show_matrix()
print("St_sac[0]")
St_sac[0].show_matrix()
#sys.exit(0)

####################
# Apply the state reordering followed by a correction to the phases.
params = {}
params["do_state_reordering"]    = 2
params["state_reordering_alpha"] = 0.0
step3.apply_state_reordering_general( St_sac, E_sac, params)
#step3.apply_phase_correction_general( St_sac )
#sys.exit(0)

# Out put sac basis data
for step in range( fstep - istep - 1 ):
    energy_sac = E_sac[step]
    nacs_sac   = 0.5j/dt * ( St_sac[step] - St_sac[step].H() )
    hvib_sac   = energy_sac - nacs_sac
    hvib_sac.real().show_matrix("%s/Hvib_sac_%d_re" % (res_dir,   istep+step))
    hvib_sac.imag().show_matrix("%s/Hvib_sac_%d_im" % (res_dir,   istep+step))
    St_sac[step].real().show_matrix("%s/St_sac_%d_re" % (res_dir, istep+step))
#END

