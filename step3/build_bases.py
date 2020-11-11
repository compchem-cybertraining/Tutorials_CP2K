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
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp




#######################################################################################
"""
 1. Read the files that have the energies, overlaps, and time-overlap matricies in the 
 Kohn-Sham basis: E_ks, S_ks, and St_ks.
"""

path = os.getcwd()
res_dir_mb   = path+"/../step2/res/"
data_dim     = 42 # rows in E_ks
active_space = range(0,int(data_dim/2)) # alpha channel only here #range(data_dim)
start_time   = 90  # initial step
finish_time  = 100 # final step + 1  
dt           = 1.0*units.fs2au

params  = { "data_set_paths" : [res_dir_mb], "data_dim":data_dim, "active_space":active_space, "isnap":start_time,  "fsnap":finish_time }

# Fetching E_ks
params.update({ "data_re_prefix" : "E_ks_",  "data_re_suffix" : "_re", "data_im_prefix" : "E_ks_",  "data_im_suffix" : "_im"  } )
E_ks_job = data_read.get_data_sets(params)
print ("\n")
E_ks_job[0][-1].show_matrix()
#sys.exit(0)

# Fetching S_ks
params.update({ "data_re_prefix" : "S_ks_", "data_re_suffix" : "_re", "data_im_prefix" : "S_ks_", "data_im_suffix" : "_im"  } )
S_ks_job = data_read.get_data_sets(params)
print ("\n")
S_ks_job[0][-1].show_matrix()
#sys.exit(0)

# Fetching St_ks
params.update({ "data_re_prefix" : "St_ks_", "data_re_suffix" : "_re", "data_im_prefix" : "St_ks_", "data_im_suffix" : "_im"  } )
St_ks_job = data_read.get_data_sets(params)
print ("\n")
St_ks_job[0][-2].show_matrix()
#sys.exit(0)

S_ks  = S_ks_job[0]
St_ks = St_ks_job[0]
step3.apply_orthonormalization_general( S_ks, St_ks )
step3.apply_phase_correction_general( St_ks )
ks_nacs = 0.5j/dt * ( St_ks[0] - St_ks[0].H() )
print("\nKS data")
S_ks[0].show_matrix()
S_ks[1].show_matrix()
St_ks[0].show_matrix()
ks_nacs.show_matrix()
#sys.exit(0)




#######################################################################################
#######################################################################################
"""
2.0. We begin by obtaining the following items for all timesteps:
    2.1. Excitation energies
    2.2. Normalized "CI" coefficients
    2.3. Update a list of unique SD states that comprise the CI-like excitations
"""
     
res_dir = "res_mb_sp"
os.system("rm -rf res_mb_sp")
os.mkdir(res_dir)
curr_step = start_time
min_band=15
max_band=35
ks_orbital_homo_index=28

# Define ks_orbital_indicies based on the min_band and max_band    
ks_orbital_indicies = list( range( min_band, max_band+1 ) )

# Update params with ks_orbital_indicies
params.update({"ks_orbital_indicies":ks_orbital_indicies})
params["logfile_directory"] = "../excitation_analysis/all_logfiles"
params["es_software"]       = "cp2k"
params["curr_step"]         = curr_step
params["isUKS"]             = 0
params["number_of_states"]  = 10
params["tolerance"]         = 0.0
params["logfile_name"]      = "step_"+str(params["curr_step"])+".out"

# The unique basis of Slater determinants basis
S_sd_job  = []
St_sd_job = []
sd_basis_states_unique = []
# The overlap matrices for excited states for each time step in this job
S_ci_job  = []
St_ci_job = []
# The configuration interaction coefficients, energies, and basis states with their spin components for each time step in this job
ci_coefficients_job   = []
ci_basis_states_job   = []
ci_energies_job       = []
spin_components_job   = []
excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = step2_many_body.get_excitation_analysis_output( params )
ci_coefficients_raw_norm = step2_many_body.normalize_ci_coefficients(ci_coefficients_raw_unnorm)
print( excitation_energies )
print( ci_basis_raw )
print( ci_coefficients_raw_unnorm )
print( spin_components )
print( ci_coefficients_raw_norm )
# Now append the extracted excitation analysis output in _job variables
ci_basis_states_job.append( ci_basis_raw )
ci_coefficients_job.append(   ci_coefficients_raw_norm )
ci_energies_job.append( excitation_energies )
spin_components_job.append( spin_components )
# Extract the uniquie SD basis states from the ci basis states
for ci_basis_state_index in range( len( ci_basis_raw ) ):
    for sd_basis_state_index in range( len( ci_basis_raw[ ci_basis_state_index ] ) ):
         sd_basis_state_and_spin = [ ci_basis_raw[ ci_basis_state_index ][ sd_basis_state_index ] , spin_components[ ci_basis_state_index ][ sd_basis_state_index ] ]
         if sd_basis_state_and_spin not in sd_basis_states_unique:
             sd_basis_states_unique.append( sd_basis_state_and_spin )
#print( "Slater determinant basis states = ", sd_basis_states_unique )
curr_step += 1

# All other steps after initial step for this job
for step in range( finish_time-start_time-1 ):
    params.update({ "curr_step":curr_step })
    excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = step2_many_body.get_excitation_analysis_output( params )
    ci_coefficients_raw_norm = step2_many_body.normalize_ci_coefficients(ci_coefficients_raw_unnorm)
    # Extract the uniquie SD basis states from the ci basis states
    for ci_basis_state_index in range( len( ci_basis_raw ) ):
        for sd_basis_state_index in range( len( ci_basis_raw[ ci_basis_state_index ] ) ):
            sd_basis_state_and_spin = [ ci_basis_raw[ ci_basis_state_index ][ sd_basis_state_index ] , spin_components[ ci_basis_state_index ][ sd_basis_state_index ] ]
        if sd_basis_state_and_spin not in sd_basis_states_unique:
            sd_basis_states_unique.append( sd_basis_state_and_spin )
    #print( "Slater determinant basis states = ", sd_basis_states_unique )
    # Now append the extracted excitation analysis output in _job variables
    ci_basis_states_job.append( ci_basis_raw )
    ci_coefficients_job.append(   ci_coefficients_raw_norm )
    ci_energies_job.append( excitation_energies )
    spin_components_job.append( spin_components )
    curr_step += 1




#######################################################################################
#######################################################################################
"""
    3.0. We need to process the unique SD bases, which involves the following

        1. Reindexing the SD bases from the native ES format to what Libra expects
        2. Sorting the SD bases at each timestep by their energies

"""
# 3.1. Reindex the SDs into the format expected by Libra
sd_states_reindexed = step2_many_body.reindex_cp2k_sd_states( ks_orbital_homo_index, ks_orbital_indicies, sd_basis_states_unique, sd_format=2 )
#print("The reindexed Slater determinants = \n", sd_states_reindexed)
#sys.exit(0)

# 3.2. Sort the SDs at each timestep
E_sd_job = []
sd_states_reindexed_sorted = []
sd_states_unique_sorted = []
SD_energy_corr = [0.0]*len(sd_states_reindexed)
for step in range( finish_time-start_time ):
    # At this step, compute the energy of the SD
    E_this_sd  = mapping.energy_mat_arb( sd_states_reindexed, E_ks_job[0][step], SD_energy_corr )
    nstates_sd = len(sd_states_reindexed)
    # Make an array of zeros, these will be overwritten with the energy of each SD
    e = np.zeros( nstates_sd )
    for state in range(nstates_sd):
        e[state] =  E_this_sd.get(state,state).real
    # Obtain the indexing fo the SDs by their energies
    reindex = np.argsort(e)
    # We need to write the energies into matrix form 
    E_sd_job.append(  CMATRIX(nstates_sd,nstates_sd) )
    sd_states_reindexed_sorted.append( [] )
    # For each SD basis, make the energy matrix and reindex the list of basis according to their energies
    for i in range(len(reindex)):
        # This is making the energy matrix
        E_sd_job[step].set(  i,i, E_this_sd.get(  int(reindex[i]), int(reindex[i])) )
        # This is reindexing the list of SD bases at this time step according to their energies 
        sd_states_reindexed_sorted[step].append( sd_states_reindexed[ int(reindex[i]) ] )
    # Make another list but omit the ground state, we will manually make this later. the below variable goes to make the T matrix
    sd_states_unique_sorted.append( [] )
    for i in range(1,len(reindex)):
        sd_states_unique_sorted[step].append( sd_basis_states_unique[ int(reindex[i])-1 ] )

# Before we make the Sd overlaps, apply a normalization to the S_ks and St_ks files, in case they for some reason are not in an orthonormal basis
step3.apply_normalization( S_ks_job[0], St_ks_job[0] )
step3.apply_phase_correction( St_ks_job[0] )

def myfunc_s( step ):
    s_sd  = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step],   S_ks_job[0][step],  use_minimal=False )
    #print(step, "s_sd as CMATRIX  <1|1> = ", s_sd.get(1,1) )
    s_sd  = data_conv.MATRIX2nparray(s_sd)
    #print(step, "s_sd as np.array <1|1> = ", s_sd[1,1])
    return s_sd

def myfunc_st( step ):
    st_sd = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step+1], St_ks_job[0][step], use_minimal=False )
    #print(step, "st_sd as CMATRIX  <1|1> = ", st_sd.get(1,1) )
    st_sd = data_conv.MATRIX2nparray(st_sd)
    #print(step, "st_sd as np.array <1|1> = ", st_sd[1,1])
    return st_sd

pool = mp.Pool( 24 )
tmp_st_sd = pool.map( myfunc_st, list(range(finish_time-start_time-1)) )
pool.close()
pool.join()

pool = mp.Pool( 24 )
tmp_s_sd = pool.map( myfunc_s, list(range(finish_time-start_time)) )
pool.close()
pool.join()

for step in range( finish_time-start_time ):
    S_sd_job.append(  data_conv.nparray2CMATRIX( tmp_s_sd[step]  ) )
for step in range( finish_time-start_time-1 ):
    St_sd_job.append( data_conv.nparray2CMATRIX( tmp_st_sd[step] ) )



#######################################################################################
#######################################################################################
'''
    4.0. We need to make the linear transformation matrices. This involes:

        1. Make the list of ci_coefficients for each step in the way Libra accepts
        2. Making the T matrix from this list keeping in mind the SD bases

'''
number_of_states = params["number_of_states"]
ci_coefficients = []
# Add one to the number of CI states because the ground state is not included yet
nSDs = len( sd_basis_states_unique ) + 1
nCIs = number_of_states + 1
SD2CI = []
for step in range( finish_time-start_time ):

    # Make the list of ci_coefficients for each step in the way Libra accepts
    ci_coefficients.append( [] )
    # Start with the ground state. This is not explicitly given by electronic strcture calculations
    ci_coefficients[step].insert( 0, [0.0] * nSDs )
    ci_coefficients[step][0][0] = 1.0

    # For each ci state for this step
    for i in range( len( ci_coefficients_job[step] ) ):
        count = 0
        # The ci wavefunction is a linear combination of SD states. Make a list of zeros the size of the number of unique
        # SD states + 1 for the ground state
        ci_coefficients[step].append( [0.0] * nSDs )
        # Exclude ground state here in the index, that info is not explicitly contained 
        # in the ci_coefficients_dynamics list from electronic structure calculations
        tmp_ci_basis_state_and_spin = []
        # For each ci_coefficient in this ci state for this step, get the ci coefficients and spin (alp or bet)
        for k in range(len(ci_coefficients_job[step][i])):
            tmp_ci_basis_state_and_spin.append( [ci_basis_states_job[step][i][k] , spin_components_job[step][i][k]] )
        # Now, loop over the SDs (excluding the ground state) to assign the coefficients
        for j in range( nSDs-1 ):
            # Check to see if one of the SDs from the list of unique SDs comprises this ci state
            if sd_states_unique_sorted[step][j] in tmp_ci_basis_state_and_spin:   
                # ok, it has found a match, now what is the index of the SD in the list of unique SDs?
                item_index = tmp_ci_basis_state_and_spin.index(sd_states_unique_sorted[step][j])
                ci_coefficients[step][i+1][j+1] = float(ci_coefficients_job[step][i][item_index])
    SD2CI.append( CMATRIX( nSDs, nCIs ) )
    for i in range( nSDs ):
        for j in range( nCIs ):
            SD2CI[step].set( i, j, ci_coefficients[step][j][i] * (1.0+0.0j) )
    # Output the transformation matrix. This is how you can double check that it worked ( it does ... :) )
    SD2CI[step].show_matrix( "%s/T_%s.txt" % (res_dir, str(step)) )




#######################################################################################
#######################################################################################
'''
    5.0. Use the SD to CI transformation matrices to convert from St_sd -> St_ci
         We also need to make the CI energy matrix from the excitation energies
         And finally make the Hvib in the CI basis
'''
# Now, compute the CI energy matrix at each-point and the mid-points
# For each step
#print("Computing the CI energy matrices....")
ci_energies_job_cmatrix = []
for step in range( finish_time-start_time ):
    ci_energies_job_cmatrix.append( CMATRIX( number_of_states + 1, number_of_states + 1 ) )
    for state in range( number_of_states + 1 ):
        if state == 0:
            ci_energies_job_cmatrix[step].set( state, state, 0.0 )
        else:
            ci_energies_job_cmatrix[step].set( state, state, ( ci_energies_job[step][state-1]  * units.ev2Ha )  )
# At the midpoints
ci_midpoint_energies = []
for step in range( finish_time-start_time-1 ):
    total_energy_mid_point = 0.0 #0.5 * ( total_energies_job[step] + total_energies_job[step+1] )
    ci_midpoint_energies.append( CMATRIX( number_of_states + 1, number_of_states + 1 ) )
    for state in range( number_of_states + 1 ):
        if state == 0:
            ci_midpoint_energies[step].set( state, state, total_energy_mid_point )
        else:
            midpoint_energy = 0.5 * ( ci_energies_job[step][state-1] + ci_energies_job[step+1][state-1] )
            ci_midpoint_energies[step].set( state, state, total_energy_mid_point + ( midpoint_energy  * units.ev2Ha )  )






# For each step make St_ci
print("Making St_ci matrices....")
for step in range( finish_time-start_time ):
    s_ci  = SD2CI[step].H() * S_sd_job[step]  * SD2CI[step]
    S_ci_job.append(  s_ci  )

for step in range( finish_time-start_time-1 ):
    st_ci = SD2CI[step].H() * St_sd_job[step] * SD2CI[step+1]
    St_ci_job.append( st_ci )

tmp_ci_nacs = (  0.5j / dt ) * CMATRIX ( ( St_ci_job[0] - St_ci_job[0].H() ).real() )
print("\nCI data before phase corrections")
S_ci_job[0].show_matrix()
S_ci_job[1].show_matrix()
St_ci_job[0].show_matrix()
tmp_ci_nacs.show_matrix()

print("\nNormalize CI basis before output")
step3.apply_orthonormalization_general( S_ci_job, St_ci_job )
params2 = {"do_state_reordering":2, "state_reordering_alpha":0}
step3.apply_state_reordering_general( St_ci_job, ci_midpoint_energies, params2 )
step3.apply_phase_correction_general( St_ci_job )

# Output CI data to res directory
#print("Outputting the CI data to the res directory..." )
for step in range( finish_time-start_time ):
    S_ci_job[step].real().show_matrix("%s/S_ci_%d_re" % (res_dir, int(start_time+step)))
    ci_energies_job_cmatrix[step].real().show_matrix("%s/E_ci_%d_re"   % (res_dir, int(start_time+step)))
for step in range( finish_time-start_time-1 ):
    St_ci_job[step].real().show_matrix("%s/St_ci_%d_re" % (res_dir, int(start_time+step)))

# Now, compute the CI NACs and compute the CI Hvib
#print("Computing and outputting the CI NACs...")
for step in range( finish_time-start_time-1 ): 
    ci_nacs = (  0.5j / dt ) * CMATRIX ( ( St_ci_job[step] - St_ci_job[step].H() ).real() )    
    ci_hvib = ci_midpoint_energies[step] - ci_nacs
    ci_hvib.real().show_matrix("%s/Hvib_ci_%d_re" % (res_dir, int( start_time+step )))
    ci_hvib.imag().show_matrix("%s/Hvib_ci_%d_im" % (res_dir, int( start_time+step )))
print("All ci steps were done successfully for this job!")
print("\nCI data")
S_ci_job[0].show_matrix()
S_ci_job[1].show_matrix()
St_ci_job[0].show_matrix()
ci_nacs.show_matrix()
#sys.exit(0)


# At the very end output the SD data
# At the midpoints
sd_midpoint_energies = []
for step in range( finish_time-start_time-1 ):
    sd_midpoint_energy = 0.5 * ( E_sd_job[step] + E_sd_job[step+1] )
    sd_midpoint_energies.append( sd_midpoint_energy )

tmp_sd_nacs = (  0.5j / dt ) * CMATRIX ( ( St_sd_job[0] - St_sd_job[0].H() ).real() )
print("\nSD data before phase corrections")
S_sd_job[0].show_matrix()
S_sd_job[1].show_matrix()
St_sd_job[0].show_matrix()
tmp_sd_nacs.show_matrix()

print("\nNormalize SD basis before output")
step3.apply_orthonormalization_general( S_sd_job, St_sd_job )
params2 = {"do_state_reordering":2, "state_reordering_alpha":0}
step3.apply_state_reordering_general( St_sd_job, sd_midpoint_energies, params2 )
step3.apply_phase_correction_general( St_sd_job )

# Now, compute the CI NACs and compute the CI Hvib
#print("Computing and outputting the SD NACs...")
for step in range( finish_time-start_time-1 ):
    sd_nacs = (  0.5j / dt ) * CMATRIX ( ( St_sd_job[step] - St_sd_job[step].H() ).real() )
    sd_hvib = sd_midpoint_energies[step] - sd_nacs
    sd_hvib.real().show_matrix("%s/Hvib_sd_%d_re" % (res_dir, int( start_time+step )))
    sd_hvib.imag().show_matrix("%s/Hvib_sd_%d_im" % (res_dir, int( start_time+step )))
print("All sd steps were done successfully for this job!")

print("\nSD data")
S_sd_job[0].show_matrix()
S_sd_job[1].show_matrix()
St_sd_job[0].show_matrix()
sd_nacs.show_matrix()
#sys.exit(0)


# END
