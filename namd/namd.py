import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
from liblibra_core import *
from libra_py import units as units
from libra_py import influence_spectrum as infsp
from libra_py import data_visualize
from libra_py import data_conv
from libra_py import data_stat
import libra_py.workflows.nbra.step4 as step4
import libra_py.workflows.nbra.lz as lz
import libra_py.workflows.nbra.decoherence_times as decoherence_times
import numpy as np
import matplotlib.pyplot as plt
import os
import multiprocessing as mp









####################
# 1. Get the Hvibs from step2 
print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../step3/res_mb_sp/")
params["Hvib_re_prefix"] = "Hvib_ci_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_ci_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 9
params["nstates"]        = 11 # total number of electronic states
params["init_times"]     = [90]
params["active_space"]   = list(range(params["nstates"])) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_mb = step4.get_Hvib2(params)
hvib_mb[0][-1].show_matrix()
print ("Length of hvib_mb is: ", len(hvib_mb[0]))
#sys.exit(0)

print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../step3/res_mb_sp/")
params["Hvib_re_prefix"] = "Hvib_sd_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_sd_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 9
params["nstates"]        = 23 # total number of electronic states
params["init_times"]     = [90]
#params["active_space"]   = list(range(params["nstates"])) # indexing is from 0!
params["active_space"]   = list(range(0,11)) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_mixed_sd = step4.get_Hvib2(params)
hvib_mixed_sd[0][-1].show_matrix()
print ("Length of hvib_mixed_sd is: ", len(hvib_mixed_sd[0]))
#sys.exit(0)






#### Helper functions

def compute_state_energies_vs_time( hvib ):
    """
    Computes the states energies vs time for a given hvib.
        hvib ( list of list of matrices ): The vibronic hamiltonian for all timesteps
    returns: a list of energies vs time for all states in hvib
    """
    nsteps   = len(hvib)
    nstates  = hvib[0].num_of_rows
    energies = []
    for state in range( nstates ):
        energies.append( [] )
        for step in range( nsteps ):
            energies[ state ].append( hvib[ step ].get( state, state ).real - hvib[ step ].get( 0, 0 ).real )
    return np.array( energies )


def run_nbra_namd_wrapper( hvib, istate, outfile_name, params ):
    """
    A small wrapper function to run the nbra namd dynamics. 
        hvib ( list of list of matrices ): The vibronic hamiltonian for all timesteps
        istate ( int ): Index of the initial state. This index is from 0
        outfile ( string ): The name of the output file
        params ( dict ): A dictionary of important dynamical parameters
    returns: A list of energy values that is the electronic energy vs time. Returned energy is in eV
    """

    print("decoherence_method = ", params["decoherence_method"])
    params["outfile"]  = outfile_name
    params["nstates"]  = hvib[0].num_of_rows
    params["istate"]   = istate
    res_nbra_namd = step4.run( [ hvib ], params )
    nbra_namd     = data_conv.MATRIX2nparray( res_nbra_namd )
    # For SH energies
    energy_nbra_namd = ( nbra_namd[:,3*params["nstates"]+2] - nbra_namd[:,1] )*units.au2ev
    return energy_nbra_namd


def run_bllz_wrapper( hvib, istate, outfile_name, params ):
    """
    A small wrapper function to run the bllz dynamics. 
        hvib ( list of list of matrices ): The vibronic hamiltonian for all timesteps
        istate ( int ): Index of the initial state. This index is from 0
        outfile ( string ): The name of the output file
        params ( dict ): A dictionary of important dynamical parameters
    returns: A list of energy values that is the electronic energy vs time. Returned energy is in eV
    """

    params["outfile"]  = outfile_name
    params["nstates"]  = hvib[0].num_of_rows
    params["istate"]   = istate
    res_bllz, P = lz.run( [hvib], params )
    bllz_namd   = data_conv.MATRIX2nparray( res_bllz )
    # Recall that for bllz, we use the schrodiger not the surface hopping energy
    energy_bllz = ( bllz_namd[:,3*params["nstates"]+1] - bllz_namd[:,1] )*units.au2ev
    # for bllz surface hopping energy
    #energy_bllz = ( bllz_namd[:,3*params["nstates"]+2] - bllz_namd[:,1] )*units.au2ev
    return energy_bllz


def get_initial_states( energies, num_istates_to_generate, prefactor, average_excitation_energy, ref_state ):
    """
    This function gets initial states for a given electronic energies vs. time
    """
   
    nstates = len(energies)    
    # Obtain energy of first step for each electronic states
    initial_excitation_energies = []
    for state in range( nstates ):
        ref_energy = energies[ref_state][0]
        initial_excitation_energy = ( energies[state][0] - ref_energy ) * units.au2ev
        initial_excitation_energies.append( initial_excitation_energy  )
        print( "State", state, initial_excitation_energy )

    rnd = Random()
    excitation_energies = [] 
    # Compute energy window for a given E_avg
    for i in range( num_istates_to_generate ):
        ksi = prefactor * rnd.normal()
        excitation_energies.append( average_excitation_energy + ksi )
    print(excitation_energies)
   
    istates = []
    for i in range( num_istates_to_generate ):
        index = min( range(len(initial_excitation_energies)), key=lambda j: abs(initial_excitation_energies[j]-excitation_energies[i]))
        istates.append( index )

    return istates




#####################
# 2. Divide up into many nuclear trajectories. These are to be consdiered our independent nuclear trajectories
nuclear_traj_parser = [
                        [0, 0], #[0, 199], [0, 399], [0, 599], [0, 799]
                      ]

num_nuclear_trajs = len( nuclear_traj_parser )
nuclear_traj_len  = params["nfiles"]
params["dt"] = 1.0 * units.fs2au

hvib_mb_trajs = []
hvib_mixed_sd_trajs = []

nstates_mb = hvib_mb[0][0].num_of_rows
nstates_sp = hvib_mixed_sd[0][0].num_of_rows

for i in range( num_nuclear_trajs ):
    hvib_mb_trajs.append( [] )
    hvib_mixed_sd_trajs.append( [] )

# Extract the trajectories
for traj in range( num_nuclear_trajs ):
    traj_id    = nuclear_traj_parser[ traj ][ 0 ]
    start_time = nuclear_traj_parser[ traj ][ 1 ] 
    for t in range( start_time, start_time + nuclear_traj_len ):
        hvib_mb_trajs[ traj ].append( hvib_mb[ traj_id ][ t ] )
        hvib_mixed_sd_trajs[ traj ].append( hvib_mixed_sd[ traj_id ][ t ] )

#####################
# 3. For each nuclear trajectory, we need to get the state energies vs. time, which is needed to get the istates
istates_mb_all_trajs  = []
energies_mb_all_trajs = []
istates_mixed_sd_all_trajs  = []
energies_mixed_sd_all_trajs = []

num_istates_to_generate = 100
for traj in range( num_nuclear_trajs ):
    
    energies_mb_all_trajs.append( compute_state_energies_vs_time( hvib_mb_trajs[ traj ] ) )
    istates_mb_all_trajs.append( get_initial_states( energies_mb_all_trajs[ traj ], num_istates_to_generate, 0.1, 9, 0 ) )

    energies_mixed_sd_all_trajs.append( compute_state_energies_vs_time( hvib_mixed_sd_trajs[ traj ] ) )
    istates_mixed_sd_all_trajs.append( get_initial_states( energies_mixed_sd_all_trajs[ traj ], num_istates_to_generate, 0.1, 9, 0 ) )

print(istates_mb_all_trajs)
print(istates_mixed_sd_all_trajs)
#sys.exit(0)


#========================================================
# Setting up now for nbra namd calculations
params["T"]          = 300.0
params["sh_method"]  = 1 
params["Boltz_opt"]  = 1 
params["nsteps"]     = nuclear_traj_len
params["init_times"] = [0]
params["ntraj"] = 1

#========================================================
# Setting up now for BLLZ calculations
# Looking on the "SE" populations - Markov chain approach
params["target_space"]       = 1
params["gap_min_exception"]  = 0
params["Boltz_opt_BL"]       = 1     # Option to incorporate hte frustrated hops into BL probabilities
params["evolve_Markov"]      = True  # Rely on the Markov approach
params["evolve_TSH"]         = False # don't care about TSH
params["extend_md"]          = False
params["extend_md_time"]     = 0




os.system("rm -rf fssh")
os.system("rm -rf ida")
os.system("rm -rf msdm")
os.system("rm -rf bllz")

os.system("mkdir fssh")
os.system("mkdir ida")
os.system("mkdir msdm")
os.system("mkdir bllz")

def myfunc( istate_index ):

    # At this point, we have the nuclear trajectories for all bases, here, we show just for the mb
    # This is a specific istate, so, here each nuclear trajectory, we need to run the dynamics at this istate
    # For each nuclear trajectory, run dynamics for this istate

    md_time = [ i for i in range( nuclear_traj_len ) ]
    md_time = np.array( md_time ) * params["dt"] * units.au2fs

    for nuc_traj_index in range( num_nuclear_trajs ):

        # Get nuc_traj_index hvbis for each basis
        mb_hvib     = hvib_mb_trajs[ nuc_traj_index ]    
        istate_mb   = istates_mb_all_trajs[ nuc_traj_index ][ istate_index  ]
        energies_mb = energies_mb_all_trajs[ nuc_traj_index ]

        mixed_sd_hvib     = hvib_mixed_sd_trajs[ nuc_traj_index ]
        istate_mixed_sd   = istates_mixed_sd_all_trajs[  nuc_traj_index ][ istate_index  ]
        energies_mixed_sd = energies_mixed_sd_all_trajs[ nuc_traj_index ]

        #"""
        ### FSSH
        params["decoherence_method"] = 0
        outfile = "fssh/_mb_fssh_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"
        mb_nbra_fssh = run_nbra_namd_wrapper( mb_hvib, istate_mb, outfile, params )
        outfile = "fssh/_mixed_sd_fssh_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"
        mixed_sd_nbra_fssh = run_nbra_namd_wrapper( mixed_sd_hvib, istate_mixed_sd, outfile, params )
        #"""
        
        #"""
        ### ID-A
        params["decoherence_method"] = 1
        outfile = "ida/_mb_ida_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"
        mb_nbra_ida = run_nbra_namd_wrapper( mb_hvib, istate_mb, outfile, params )
        outfile = "ida/_mixed_sd_ida_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"
        mixed_sd_nbra_ida = run_nbra_namd_wrapper( mixed_sd_hvib, istate_mixed_sd, outfile, params )
        #"""

        #"""
        ### mSDM
        params["decoherence_method"] = 2
        outfile = "msdm/_mb_msdm_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"
        mb_nbra_msdm = run_nbra_namd_wrapper( mb_hvib, istate_mb, outfile, params )
        outfile = "msdm/_mixed_sd_msdm_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"
        mixed_sd_nbra_msdm = run_nbra_namd_wrapper( mixed_sd_hvib, istate_mixed_sd, outfile, params )
        #"""

        #"""
        ### bllz
        outfile = "bllz/_mb_bllz_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"
        mb_nbra_bllz = run_bllz_wrapper( mb_hvib, istate_mb, outfile, params )
        outfile = "bllz/_mixed_sd_bllz_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"
        mixed_sd_nbra_bllz = run_bllz_wrapper( mixed_sd_hvib, istate_mixed_sd, outfile, params )
        #"""





pool = mp.Pool( 32 )
pool.map( myfunc, list( range( num_istates_to_generate ) ) )
pool.close()
pool.join()

