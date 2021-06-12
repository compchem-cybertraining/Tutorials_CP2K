from libra_py import data_stat
from libra_py import CP2K_methods
from libra_py.workflows.nbra import step2_many_body 
import os
import sys
import multiprocessing as mp
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.colors
import math
import time




###############################################################
"""
  1. Number of states, tolerance, curr_step
"""
params = { }
params["number_of_states"] = 10
params["tolerance"] = 0.0
thermal = False #True






###############################################################
"""
  2. Extract CI-like coefficients
"""

if thermal == True:
    logfiles = glob.glob('all_logfiles/*.log')

elif thermal == False:
    logfiles = glob.glob('../tddft/*.log')

ci_coeffs = []
for logfile in logfiles:
    params.update({"logfile_name": logfile})
    excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = CP2K_methods.read_cp2k_tddfpt_log_file( params ) 
    ci_coefficients_raw_norm = step2_many_body.normalize_ci_coefficients(ci_coefficients_raw_unnorm)
    for j in range(len(ci_coefficients_raw_norm)):
        for k in range(len(ci_coefficients_raw_norm[j])):
            ci_coefficients_raw_norm[j][k] = ci_coefficients_raw_norm[j][k]**2
    ci_coeffs.append(ci_coefficients_raw_norm)






###############################################################
"""
  3. Post process the Ci-like coefficients, plot
"""

nsteps = len(ci_coeffs)
nstates = params["number_of_states"]
nsds = 5
coeffs = []
coeffs_avg   = []
coeffs_error = []

plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
plt.subplot(1,1,1)


if thermal == True:
    plt.title("Adamantane, 300 K", fontsize=10)
else:
    plt.title("Adamantane, 0 K", fontsize=10)


plt.ylim(0,1)
plt.xlabel('State Index', fontsize=10)
plt.ylabel('< c$_{i}^2$ >',   fontsize=10)
plt.xticks([0,2,4,6,8,10])

for state in range(nstates):

    coeffs.append( [] )
    coeffs_avg.append( [] )
    coeffs_error.append( [] )

    for sd in range( nsds ):

        coeffs[state].append( [] )
        coeffs_avg[state].append( [] )
        coeffs_error[state].append( [] )

        for step in range( nsteps ):
            if len( ci_coeffs[step][state] ) < nsds and sd > len( ci_coeffs[step][state] )-1:
                coeffs[state][sd].append( 0.0 )
            else:
                coeffs[state][sd].append( ci_coeffs[step][state][sd] )
     
        mb_coeff_avg, mb_coeff_std = data_stat.scalar_stat( coeffs[state][sd] )
        coeffs_avg[state][sd].append( mb_coeff_avg )
        coeffs_error[state][sd].append( 1.96 * mb_coeff_std / math.sqrt(nsteps) )

        if sd == 0:
            print("std = ", mb_coeff_std)     

        if sd == 0:
            plt.plot(     state+1, mb_coeff_avg, color="black", marker='s', markerfacecolor='green', markeredgewidth=0.4, markersize=5)
            plt.errorbar( state+1, mb_coeff_avg, yerr=coeffs_error[state][sd], linestyle="None", color='black')

        elif sd == 1:
            plt.plot(     state+1, mb_coeff_avg, color="black", marker='s', markerfacecolor='blue',  markeredgewidth=0.4, markersize=5)
            plt.errorbar( state+1, mb_coeff_avg, yerr=coeffs_error[state][sd], linestyle="None", color='black')

        elif sd == 2:
            plt.plot(     state+1, mb_coeff_avg, color="black", marker='s', markerfacecolor='red',  markeredgewidth=0.4, markersize=5)
            plt.errorbar( state+1, mb_coeff_avg, yerr=coeffs_error[state][sd], linestyle="None", color='black')

        elif sd == 3:
            plt.plot(     state+1, mb_coeff_avg, color="black", marker='s', markerfacecolor='purple',  markeredgewidth=0.4, markersize=5)
            plt.errorbar( state+1, mb_coeff_avg, yerr=coeffs_error[state][sd], linestyle="None", color='black')

        else:
            plt.plot(     state+1, mb_coeff_avg, color="black", marker='s', markerfacecolor='cyan',  markeredgewidth=0.4, markersize=5)
            plt.errorbar( state+1, mb_coeff_avg, yerr=coeffs_error[state][sd], linestyle="None", color='black')
 
plt.tight_layout()
if thermal == True:
    plt.savefig('Adamantane_300K.png', dpi=300)
else:
    plt.savefig('Adamantane_0K.png', dpi=300)

