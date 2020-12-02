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

step3.apply_normalization( S_ks[0], St_ks[0] )
step3.apply_phase_correction( St_ks[0] )
#sys.exit(0)


ks_homo_index = 28
min_band_ks   = 28 
max_band_ks   = 40
ks_orbital_indicies = list( range( min_band_ks, max_band_ks + 1 ) )
sd_basis_states_unique = [
                            # Phi 2  
                            [ [28,30], "alp" ],
                            # Phi 3  
                            [ [28,31], "alp" ],
                            # Phi 4  
                            [ [28,32], "alp" ],
                            # Phi 5  
                            [ [28,33], "alp" ],
                            # Phi 6  
                            [ [28,34], "alp" ],
                            # Phi 7  
                            [ [28,35], "alp" ],
                            # Phi 1  
                            [ [28,29], "alp" ],
                         ]    


# 2.3. Reindex the single-particle excitations into a format expected by Libra
sd_states_reindexed = step2_many_body.reindex_cp2k_sd_states( ks_homo_index, ks_orbital_indicies, sd_basis_states_unique, sd_format=2 )

# 2.4. Order / sort the single-particle excitations at each timestep by energy or identity
E_sd, sd_states_unique_sorted, sd_states_reindexed_sorted, reindex_nsteps = step3.sort_unique_SD_basis( E_ks[0], sd_basis_states_unique, sd_states_reindexed, istep, fstep, sorting_type="energy" )
#E_sd, sd_states_unique_sorted, sd_states_reindexed_sorted, reindex_nsteps = step3.sort_unique_SD_basis( E_ks[0], sd_basis_states_unique, sd_states_reindexed, istep, fstep, sorting_type="identity" )

#print(reindex_nsteps)
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


#######################################################################################
"""
 4. Make and output the Hvib in the basis of single-particle excitations
    5.1. Compute single-particle excitation energies at the midpoints
    5.2. Apply orthonormalization to the many-body basis, apply state reordering, 
    and apply phase corrections
    5.3. Make the hvib in the basis of single-particle excitations
"""
#######################################################################################

# 5.1. Compute single-particle excitation energies at the midpoints
sd_midpoint_energies = []
for step in range( fstep - istep - 1 ):
    sd_midpoint_energy = 0.5 * ( E_sd[step] + E_sd[step+1] )
    sd_midpoint_energies.append( sd_midpoint_energy )
print("\nNormalize SD basis before output")

# 5.2. Apply orthonormalization to the many-body basis, apply state reordering, and apply phase corrections
step3.apply_orthonormalization_general( S_sd, St_sd )
params2 = {"do_state_reordering":2, "state_reordering_alpha":0}
step3.apply_state_reordering_general( St_sd, sd_midpoint_energies, params2 )

# 5.3. Make the Hvib in the basis of single-particle excitations
for step in range( fstep - istep - 1 ):
    sd_nacs = (  0.5j / dt ) * CMATRIX ( ( St_sd[step] - St_sd[step].H() ).real() )
    sd_hvib = sd_midpoint_energies[step] - sd_nacs
    sd_hvib.real().show_matrix("%s/Hvib_sd_%d_re" % (res_dir, int( istep+step )))
    sd_hvib.imag().show_matrix("%s/Hvib_sd_%d_im" % (res_dir, int( istep+step )))

for step in range( fstep - istep ):
    S_sd[step].real().show_matrix("%s/S_sd_%d_re" % (res_dir, int( istep+step )))
for step in range( fstep - istep - 1 ):
    St_sd[step].real().show_matrix("%s/St_sd_%d_re" % (res_dir, int( istep+step )))
# END


