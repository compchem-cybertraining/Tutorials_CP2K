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
# Helper functions
def get_step2_mb_sp_properties( params ):
    """
        This function extracts information from es output files. It currently works for
        cp2k, dftb+, gaussian

        params ( dictionary )

            params["ks_orbital_indicies"] (list)   - ks orbitals that were used in step2. Ex) [15,16,17, ... , 32,33,34]
            params["logfile_directory"]   (string) - "../excitation_analysis/all_logfiles"
            params["es_software"]         (string) - "cp2k"
            params["isUKS"]               (int)    - 0 for spin-unpolarized, 1 for spin-polarized
            params["number_of_states"]    (int)    - number of ci states to extract
            params["tolerance"]           (float)  - cutoff for SD contribution
            params["start_time"]          (int)    - file index to start reading from
            params["finish_time"]         (int)    - file index to stop  reading from
   
        returns S_sd_job, St_sd_job, sd_basis_states_unique, ci_basis_states_job, ci_coefficients_job, ci_energies_job, spin_components_job
         
        S_sd_job  (list) - overlaps of SD at each step
        St_sd_job (list) - time-overlaps of SD at each step
        sd_basis_states_unique (list) - 1 of each of the SP transitions (and its spin) that made up the considered CI states
        ci_basis_states_job (list) - similar to sd_basis_states_unique, but all SP transitions encountered
        ci_coefficients_job (list) - all ci coefficients encountered
        ci_energies_job (list) - excitation energies
        spin_components_job (list) - just the spin components (alpha or beta excitaiton?)
    """


    start_time  = params["start_time"]
    finish_time = params["finish_time"]

    curr_step = start_time

    # Define ks_orbital_indicies based on the min_band and max_band    

    # Update params with ks_orbital_indicies
    params["curr_step"]           = curr_step
    params["logfile_name"]        = "step_"+str(params["curr_step"])+".out"

    S_sd_job, St_sd_job = [], []
    ci_energies_job, ci_basis_states_job, ci_coefficients_job = [], [], []
    spin_components_job = []
    sd_basis_states_unique = []

    excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = step2_many_body.get_excitation_analysis_output( params )
    ci_coefficients_raw_norm = step2_many_body.normalize_ci_coefficients(ci_coefficients_raw_unnorm)

    print( excitation_energies )
    print( ci_basis_raw )
    print( ci_coefficients_raw_unnorm )
    print( spin_components )
    print( ci_coefficients_raw_norm )

    # Now append the extracted excitation analysis output to the respective lists
    ci_basis_states_job.append( ci_basis_raw )
    ci_coefficients_job.append( ci_coefficients_raw_norm )
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

    return S_sd_job, St_sd_job, sd_basis_states_unique, ci_basis_states_job, ci_coefficients_job, ci_energies_job, spin_components_job





def sort_unique_SD_basis( E_ks, sd_states_unique, sd_states_reindexed, sorting_type ):
    """
        This function computes the energies of the SP transitions (according to the sum of 1 electron terms) - no J or K
        It then may sort the order of the sd_states either based on their energy at each timestep
       
        E_ks (list of CMATRIX) - KS orbital energies
        sd_states_unique (list of lists) - all SP transitions and which spin it was
                                           Ex) [ [ ['28 29'], ['alp'] ]. [ ['28 30'], ['alp'] ] ]
        sd_states_reindexed (list of lists) - sd_states_unique but in internal  Libra notation 
                                           Ex) [ [1,-1,3,-2], [3,-1,2,-2] ]
        sorting_type ( (string) - "energy"   - sort by energy
                                  "identity" - sort by identity

        Returns E_sd, sd_states_unique_sorted, sd_states_reindexed_sorted
        
        E_sd (list of CMATRIX) - SD energies at each timestep
        sd_states_unique_sorted - all SP transitions and which spin it is, but now sorted somehow ^
        sd_states_reindexed_sorted - sd_states_unique_sorted, but in Libra's notation
    """

    E_sd = []
    sd_states_reindexed_sorted = []
    sd_states_unique_sorted = []
    SD_energy_corr = [0.0]*len(sd_states_reindexed)
    nstates_sd = len(sd_states_reindexed)

    for step in range( finish_time-start_time ):

        # Append a CMATRIX of dimension nstates_sd x nstates_sd
        E_sd.append( CMATRIX( nstates_sd, nstates_sd) )

        # At this step, compute the energy of the SD
        E_this_sd = mapping.energy_mat_arb( sd_states_reindexed, E_ks[step], SD_energy_corr )

        # Make a list for the final ordering of the sd_states_unique.
        # This will not contain the ground state, which we will manually add later. 
        sd_states_unique_sorted.append( [] )
        sd_states_reindexed_sorted.append( [] )

        if sorting_type == "identity":

            for state in range(nstates_sd):

                # This is making the energy matrix. And the "sorted" sds (written in Libra's input format) are just appened
                # to sd_states_reindexed_sorted[step]. For identity ordering - no energy sorting is done!
                E_sd[step].set(  state, state, E_this_sd.get( state, state ) )
                sd_states_reindexed_sorted[step].append( sd_states_reindexed[ state ] )
                print( sd_states_reindexed_sorted[step][ state ], ( E_sd[step].get( state, state ) - E_sd[step].get( 0, 0 ) ).real * units.au2ev )

                # This is reindexing the list of SD bases at this time step according to their energies 
                # We are adding the ground state SD later, so skip it for now. In this list sd_states_unique,
                # the ground state is not there - this list is the single-particle transitions (and spin) given by the
                # ES software. So, for example, if nstates_sd = 4, we take only the first 3, because the ground state is not
                # in sd_states_unique
                # Ex) sd_states_unique = [  [ ['28,29'], ['alp'] ] , [ ['27,29'], ['alp'] ] , [ ['26,29'], ['alp'] ] ]    
                # 28 = homo
                if state < nstates_sd-1:
                    sd_states_unique_sorted[step].append( sd_states_unique[ state ] )

        elif sorting_type == "energy":

            # Make an array of zeros, these will be overwritten with the energy of each SD
            e = np.zeros( nstates_sd )
            for state in range(nstates_sd):
                e[state] =  E_this_sd.get(state,state).real
                # Obtain the indexing fo the SDs by their energies
            reindex = np.argsort(e)

            # For each SD basis, make the energy matrix and reindex the list of basis according to their energies
            for i in range(len(reindex)):
                # This is making the energy matrix
                E_sd[step].set( i, i, E_this_sd.get(  int(reindex[i]), int(reindex[i])) )
                # This is reindexing the list of SD bases at this time step according to their energies 
                sd_states_reindexed_sorted[step].append( sd_states_reindexed[ int(reindex[i]) ] )
                print( sd_states_reindexed_sorted[step][i], ( E_sd[step].get( i, i ) - E_sd[step].get( 0, 0 ) ).real * units.au2ev )

            for i in range(1,len(reindex)):
                sd_states_unique_sorted[step].append( sd_states_unique[ int(reindex[i])-1 ] )


    return E_sd, sd_states_unique_sorted, sd_states_reindexed_sorted





def make_T_matricies( ci_coefficients, ci_basis_states, spin_components, sd_states_unique_sorted, params, verbose=1):

    number_of_states = params["number_of_states"]
    ci_coefficients_libra = []
    nSDs = len( sd_states_unique_sorted[0] ) + 1
    # Add one to the number of CI states because the ground state is not included yet
    nCIs = number_of_states + 1
    SD2CI = []

    for step in range( finish_time-start_time ):

        # Make the list of ci_coefficients for each step in the way Libra accepts
        ci_coefficients_libra.append( [] )
        # Start with the ground state. This is not explicitly given by electronic strcture calculations
        ci_coefficients_libra[step].insert( 0, [0.0] * nSDs )
        ci_coefficients_libra[step][0][0] = 1.0

        # For each ci state for this step
        for i in range( len( ci_coefficients[step] ) ):
            count = 0
            # The ci wavefunction is a linear combination of SD states. Make a list of zeros the size of the number of unique
            # SD states + 1 for the ground state
            ci_coefficients_libra[step].append( [0.0] * nSDs )
            # Exclude ground state here in the index, that info is not explicitly contained 
            # in the ci_coefficients_dynamics list from electronic structure calculations
            tmp_ci_basis_state_and_spin = []
            # For each ci_coefficient in this ci state for this step, get the ci coefficients and spin (alp or bet)
            for k in range(len(ci_coefficients[step][i])):
                tmp_ci_basis_state_and_spin.append( [ci_basis_states[step][i][k] , spin_components[step][i][k]] )
            # Now, loop over the SDs (excluding the ground state) to assign the coefficients
            for j in range( nSDs-1 ):
                # Check to see if one of the SDs from the list of unique SDs comprises this ci state
                if sd_states_unique_sorted[step][j] in tmp_ci_basis_state_and_spin:   
                    # ok, it has found a match, now what is the index of the SD in the list of unique SDs?
                    item_index = tmp_ci_basis_state_and_spin.index(sd_states_unique_sorted[step][j])
                    ci_coefficients_libra[step][i+1][j+1] = float(ci_coefficients[step][i][item_index])


        # Sanity check. Make sure sum of squared elements of columns == 1:
        for i in range( nCIs ):
            check_norm = 0
            for j in range( nSDs ):
                check_norm += ci_coefficients_libra[step][i][j]**2
            if verbose == 1:
                print("Step", step, "state", i, "check_norm", check_norm)
            if check_norm < 0.99 or check_norm > 1.01:
                print("Warning: Step, ", step) 
                print("Column ", i, "in SD2Ci (T) matrix has norm either < 0.99 or > 1.01")
                print("Exiting now")
                sys.exit(0)

        SD2CI.append( CMATRIX( nSDs, nCIs ) )
        for i in range( nSDs ):
            for j in range( nCIs ):
                SD2CI[step].set( i, j, ci_coefficients_libra[step][j][i] * (1.0+0.0j) )
 
        # Output the transformation matrix. This is how you can double check that it worked ( it does ... :) )
        SD2CI[step].show_matrix( "%s/T_%s.txt" % (res_dir, str(step)) )

    return SD2CI




def compute_ci_energies_midpoint( ci_energies, params ):
    """
    Function compute the excitation energies energies at the midpoint from a list of excitation energies at each step. 
    At each step, there are many electronic states. This function takes a list as an input, and is meant to be used 
    in the NBRA workflow calculatiosn where lists may be more convenient than matricies. 

    Energies are assumed to be energies from TDDFT calculatons. This function adds zero as the ground state total energy

    ci_energies - list of list of energies 
    params ( dictionary )
        params["number_of_states"] - number of ci states

    Returns ci_midpoint_energies
    ci_midpoint_energies (list of CMATRIX) - energies in Ha. Ground state energy is set to zero
    """

    nstates     = params["number_of_states"]
    start_time  = params["start_time"]
    finish_time = params["finish_time"]

    # Now, compute the CI energy matrix at each-point and the mid-points
    # For each step
    #print("Computing the CI energy matrices....")
    ci_energies_cmatrix = []
    for step in range( finish_time-start_time ):
        ci_energies_cmatrix.append( CMATRIX( params["number_of_states"] + 1, params["number_of_states"] + 1 ) )
        for state in range( params["number_of_states"] + 1 ):
            if state == 0:
                ci_energies_cmatrix[step].set( state, state, 0.0 )
            else:
                ci_energies_cmatrix[step].set( state, state, ( ci_energies[step][state-1]  * units.ev2Ha )  )

    # At the midpoints
    ci_midpoint_energies = []
    for step in range( finish_time-start_time-1 ):
        total_energy_mid_point = 0.0 #0.5 * ( total_energies[step] + total_energies[step+1] )
        ci_midpoint_energies.append( CMATRIX( params["number_of_states"] + 1, params["number_of_states"] + 1 ) )
        for state in range( params["number_of_states"] + 1 ):
            if state == 0:
                ci_midpoint_energies[step].set( state, state, total_energy_mid_point )
            else:
                midpoint_energy = 0.5 * ( ci_energies[step][state-1] + ci_energies[step+1][state-1] )
                ci_midpoint_energies[step].set( state, state, total_energy_mid_point + ( midpoint_energy  * units.ev2Ha )  )
    
    return ci_midpoint_energies




def myfunc_s( step ):
    """
    Function used for the making the computation of overlaps parallel via python multiprocessing 
    """
    s_sd  = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step], S_ks[0][step],  use_minimal=False )
    #print(step, "s_sd as CMATRIX  <1|1> = ", s_sd.get(1,1) )
    s_sd  = data_conv.MATRIX2nparray(s_sd)
    #print(step, "s_sd as np.array <1|1> = ", s_sd[1,1])
    return s_sd

def myfunc_st( step ):
    """
    Function used for the making the computation of overlaps parallel via python multiprocessing 
    """
    st_sd = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step+1], St_ks[0][step], use_minimal=False )
    #print(step, "st_sd as CMATRIX  <1|1> = ", st_sd.get(1,1) )
    st_sd = data_conv.MATRIX2nparray(st_sd)
    #print(step, "st_sd as np.array <1|1> = ", st_sd[1,1])
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
res_dir_mb   = path+"/../step2/res/"
data_dim     = 62 # rows in E_ks
active_space = range(0,int(data_dim/2)) # alpha channel only here #range(data_dim)
start_time   = 150 # initial step
finish_time  = 201 # final step + 1  
dt           = 1.0*units.fs2au

params = { "data_set_paths" : [res_dir_mb], "data_dim":data_dim, "active_space":active_space, "isnap":start_time,  "fsnap":finish_time }
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








#######################################################################################
"""
 2. Read the TDDFT output files and get all the needed information 
    2.1. Set parameters
    2.2. Get the information from the TDDFT calculations
    2.3. Reindex the single-particle excitations into a format expected by Libra
    2.4. Order / sort the single-particle excitations at each timestep by energy or identity
"""
#######################################################################################

# 2.1. Set parameters
min_band_ks   = 10
max_band_ks   = 40
ks_homo_index = 28
ks_orbital_indicies = list( range( min_band_ks, max_band_ks + 1 ) )
params["ks_orbital_indicies"] = ks_orbital_indicies
params["logfile_directory"]   = "../excitation_analysis/all_logfiles"
params["es_software"]         = "cp2k"
params["isUKS"]               = 0
params["number_of_states"]    = 10 
params["tolerance"]           = 0.0
params["start_time"]  = start_time
params["finish_time"] = finish_time

# 2.2. Get the information from the TDDFT calculations
S_sd, St_sd, sd_basis_states_unique, ci_basis_states, ci_coefficients, ci_energies, spin_components = get_step2_mb_sp_properties( params )

# 2.3. Reindex the single-particle excitations into a format expected by Libra
sd_states_reindexed = step2_many_body.reindex_cp2k_sd_states( ks_homo_index, ks_orbital_indicies, sd_basis_states_unique, sd_format=2 )

# 2.4. Order / sort the single-particle excitations at each timestep by energy or identity
E_sd, sd_states_unique_sorted, sd_states_reindexed_sorted = sort_unique_SD_basis( E_ks[0], sd_basis_states_unique, sd_states_reindexed, sorting_type="energy" )








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
tmp_s_sd = pool.map( myfunc_s, list(range(finish_time-start_time)) )
pool.close()
pool.join()

# 3.2. Time-overlaps
pool = mp.Pool( 24 )
tmp_st_sd = pool.map( myfunc_st, list(range(finish_time-start_time-1)) )
pool.close()
pool.join()

# 3.3. Convert the data types to Libra's CMATRIX
for step in range( finish_time-start_time ):
    S_sd.append(  data_conv.nparray2CMATRIX( tmp_s_sd[step]  ) )
for step in range( finish_time-start_time-1 ):
    St_sd.append( data_conv.nparray2CMATRIX( tmp_st_sd[step] ) )








#######################################################################################
"""
 4. 
    4.1. Take the list of excitation energies at each timestep, and compute the midpoints
    4.2. Make the transformation matrix from the single-particle to many-body basis at each 
         time step. 
    4.3. Transform from single-particle to many-body (CI) basis
    4.4. Apply orthonormalization to the many-body basis, apply state reordering, 
    and apply phase corrections
    4.5. Output many-body basis overlaps and time-overlaps to the res directory
    4.6. Make the Hvib in the many-body basis
"""
#######################################################################################


# 4.1. Take the list of excitation energies at each timestep, and compute the midpoints 
ci_midpoint_energies = compute_ci_energies_midpoint( ci_energies, params )

# 4.2. Make the transformation matrix from the single-particle to many-body basis at each time step
SD2CI = make_T_matricies( ci_coefficients, ci_basis_states, spin_components, sd_states_unique_sorted, params )

# 4.3. Transform from single-particle to many-body (CI) basis
S_ci, St_ci  = [], []
for step in range( finish_time-start_time ):
    s_ci  = SD2CI[step].H() * S_sd[step]  * SD2CI[step]
    S_ci.append(  s_ci  )
for step in range( finish_time-start_time-1 ):
    st_ci = SD2CI[step].H() * St_sd[step] * SD2CI[step+1]
    St_ci.append( st_ci )

# 4.4. Apply orthonormalization to the many-body basis, apply state reordering, apply phase corrections
step3.apply_orthonormalization_general( S_ci, St_ci )
params2 = {"do_state_reordering":2, "state_reordering_alpha":0}
step3.apply_state_reordering_general( St_ci, ci_midpoint_energies, params2 )
step3.apply_phase_correction_general( St_ci )

# 4.5. Output many-body basis overlaps and time-overlaps to the res directory
print("Outputting the CI data to the res directory..." )
for step in range( finish_time-start_time ):
    S_ci[step].real().show_matrix("%s/S_ci_%d_re" % (res_dir, int(start_time+step)))
for step in range( finish_time-start_time-1 ):
    St_ci[step].real().show_matrix("%s/St_ci_%d_re" % (res_dir, int(start_time+step)))

# 4.6. Make the Hvib in the many-body basis
for step in range( finish_time-start_time-1 ): 
    ci_nacs = (  0.5j / dt ) * CMATRIX ( ( St_ci[step] - St_ci[step].H() ).real() )    
    ci_hvib = ci_midpoint_energies[step] - ci_nacs
    ci_hvib.real().show_matrix("%s/Hvib_ci_%d_re" % (res_dir, int( start_time+step )))
    ci_hvib.imag().show_matrix("%s/Hvib_ci_%d_im" % (res_dir, int( start_time+step )))








#######################################################################################
"""
 5. Make and output the Hvib in the basis of single-particle excitations
    5.1. Compute single-particle excitation energies at the midpoints
    5.2. Apply orthonormalization to the many-body basis, apply state reordering, 
    and apply phase corrections
    5.3. Make the hvib in the basis of single-particle excitations
"""
#######################################################################################

# 5.1. Compute single-particle excitation energies at the midpoints
sd_midpoint_energies = []
for step in range( finish_time-start_time-1 ):
    sd_midpoint_energy = 0.5 * ( E_sd[step] + E_sd[step+1] )
    sd_midpoint_energies.append( sd_midpoint_energy )
print("\nNormalize SD basis before output")

# 5.2. Apply orthonormalization to the many-body basis, apply state reordering, and apply phase corrections
step3.apply_orthonormalization_general( S_sd, St_sd )
params2 = {"do_state_reordering":2, "state_reordering_alpha":0}
step3.apply_state_reordering_general( St_sd, sd_midpoint_energies, params2 )
step3.apply_phase_correction_general( St_sd )

# 5.3. Make the Hvib in the basis of single-particle excitations
for step in range( finish_time-start_time-1 ):
    sd_nacs = (  0.5j / dt ) * CMATRIX ( ( St_sd[step] - St_sd[step].H() ).real() )
    sd_hvib = sd_midpoint_energies[step] - sd_nacs
    sd_hvib.real().show_matrix("%s/Hvib_sd_%d_re" % (res_dir, int( start_time+step )))
    sd_hvib.imag().show_matrix("%s/Hvib_sd_%d_im" % (res_dir, int( start_time+step )))
# END


