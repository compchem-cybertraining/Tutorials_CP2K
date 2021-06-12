#***********************************************************
# * Copyright (C) 2020 Brendan Smith, Mohammad Shakiba, Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/
import os
import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import data_stat
from libra_py import units
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.optimize import curve_fit


### Helpfer functions

# From Story
def distrib(data, xmin, xmax, dx):
    # take real components
    for i in range(len(data)):
        data[i] = data[i].real            
    data = DATA(data)    
    #============= Build the bin support ===========
    bin_support = []
    x = xmin
    while x <= xmax:        
        bin_support.append(x)
        x = x + dx
    dens, cum = data.Calculate_Distribution(bin_support)
    return bin_support, dens, cum


# From Story
def plot_pop_dens(data, avg_data, dt, dx, emin, emax, title, basis, printout=0):

    t_points = len(data[0])
    ntraj = len(data)
    
    dens_list = [] # list of lists of densities at each time interval
    # calculating probility densities at each state and time interval
    for j in range(t_points):
        pop = []
        for k in range(ntraj):
            pop.append(data[k][j])
        if printout:
            print(pop)
        # calculate density at this time point
        bin_supp, dens, cum = distrib(pop, emin, emax, dx)
        if printout:
            print(dens)
        dens_list.append(dens)
 
    #converting densities into surface as np array       
    surface = np.array(dens_list)
    # originally each row is a time step, but we want time steps as columns
    surface = np.transpose(surface)
    # the matrix here is upside down, so has to be flipped vertically 
    surface = np.flipud(surface)
    
    # normalized
    max_value = np.max(surface)
    if max_value != 0:
        surface /=  max_value

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)   
    #plt.title("Cubic, "+title, fontsize=12)

    if basis == "mb":
        plt.title(title+", Many-body", fontsize=12)
    elif basis == "mixed_sd":
        plt.title(title+", Single-particle", fontsize=12)

    plt.xlabel('Time, fs', fontsize=10)
    plt.ylabel('Excess energy, eV', fontsize=10)
    plt.ylim(0,2.0)
    plt.imshow(surface, cmap="hot", aspect = 'auto', extent=[0,((t_points-1)*dt),emin,emax], interpolation='gaussian')
    plt.plot(xdata, avg_data, linewidth=2, color="white", label="")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("Energy_decay_"+title+"_"+basis+".png")







####################
# 1. Define the needed helper functions. This includes the fitting function and a wrapper funciton to call it
def stretch_compressed(t, tau, beta, E0):
    return E0*np.exp( -(t/tau)**beta )

def single_exponential(t, tau, E0):
    return E0*np.exp( -(t/tau) )

def gaussian_exponential(t, a, tau1, tau2, E0):
    return a*np.exp( -( t/tau1 ) ) + (E0-a)*np.exp( -( t/tau2 )**2 )


def fit_data( xdata, ydata, function_option ):

    ydata_fit = []

    # stretch-compressed 
    if function_option == 0:
        popt, pcov = curve_fit( stretch_compressed, xdata, ydata, bounds=([0.0, 0.0, ydata[0]-0.001], [np.inf, np.inf, ydata[0]+0.001]))
        tau, beta, E0 = popt
        for t in range( len(xdata) ):
            ydata_fit.append( stretch_compressed( xdata[t], *popt ) )
        residuals  = ydata - stretch_compressed(xdata, *popt)

    # single-exponential
    elif function_option == 1:
        popt, pcov = curve_fit( single_exponential, xdata, ydata, bounds=([0.0, ydata[0]-0.001], [np.inf, ydata[0]+0.001]))
        tau, E0 = popt
        for t in range( len(xdata) ):
            ydata_fit.append( single_exponential( xdata[t], *popt ) )   
        residuals  = ydata - single_exponential(xdata, *popt)

    # gaussian-exponential
    elif function_option == 3:
        popt, pcov = curve_fit( gaussian_exponential, xdata, ydata, bounds=([0.0, 0.0, 0.0, ydata[0]-0.001], [np.inf, np.inf, np.inf, ydata[0]+0.001]))
        a, tau1, tau2, E0 = popt
        for t in range( len(xdata) ):
            ydata_fit.append( gaussian_exponential( xdata[t], *popt ) )
        residuals  = ydata - gaussian_exponential(xdata, *popt)


    ss_res, ss_tot = np.sum(residuals**2), np.sum((ydata - np.mean(ydata))**2)
    r_squared  = 1.0 - (ss_res / ss_tot)

    """
    # Uncomment the below section to plot the fit. This can be done to check a special case of intrest, if need be 
    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    plt.title("TEST", fontsize=8) 
    plt.xlabel('Time, fs',   fontsize=8)
    plt.plot(xdata, ydata,     linewidth=1.5, color="black", label="data")
    plt.plot(xdata, ydata_fit, linewidth=1,   color="green", label="fit") 
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig("tmp.png")
    """

    # stretch compressed
    if function_option == 0:
        return ydata_fit, tau, beta, E0, r_squared

    # single exponential
    elif function_option == 1:
        return ydata_fit, tau, E0, r_squared

    # gaussian exponential
    elif function_option == 3:
        return ydata_fit, a, tau1, tau2, E0, r_squared





####################
# 2. Set paramters for reading / sorting the data
namd_length = 50
dt = 1.0 #fs
decoherence_options       = ["fssh","ida","msdm","bllz"]
decoherence_options_names = ["FSSH","ID-A","mSDM","BLLZ"]
fit_option = 3
# 0 = "Hot_energy decay"
dynamics_option   = 0 
num_istates       = 100
num_nuclear_trajs = 1
taus  = []
betas = []
dyns  = []
decoherence_option_count = 0
# Hot energy
if dynamics_option == 0:
    basis_options = ["mb","mixed_sd"]




#####################
# 3. Enter the "ultimate" fitting loop. For each decoherence option, we will perofrm dynamics using both the SD and MB bases
for decoherence_option in decoherence_options:

    taus.append(  [] )
    betas.append( [] )
    dyns.append(  [] )

    # 3.1. For each basis option choice (SD and MB) we have a number of initial conditions, we will obtain a tau for each
    basis_option_count = 0
    for basis_option in basis_options:
        
        taus[decoherence_option_count].append(  [] )
        betas[decoherence_option_count].append( [] )
        dyns[decoherence_option_count].append(  [] )     

        # 3.2. Now, the initial conditions are for 4 nuclear sub-trajectories, each having 3 initial "istate" conditions
        for nuc_traj_index in range( num_nuclear_trajs ):

            # 3.3.. Again, each istate is also an initial condition, here. So for each initial condition, we collect a tau
            for istate_index in range(num_istates):

                xdata     = []
                ydata     = []
                ydata_fit = []

                filename = "../"+decoherence_option+"/_"+basis_option+"_"+decoherence_option+"_"+str( nuc_traj_index )+"_"+str( istate_index )+".txt"

                f = open(filename)
                A = f.readlines()
                sz = len(A)
                f.close()

                for i in range( namd_length ):

                    namd_data_line = A[i].strip().split()

                    # Hot Energy Decay
                    if dynamics_option == 0:

                        if decoherence_options[decoherence_option_count] == "bllz":

                            # mb dynamics
                            if basis_option_count == 0:
                                y = float( namd_data_line[ 34 ] ) - float( namd_data_line[4] )
                                if y < 0:
                                    y = 0.00

                            # mixed_sd dynamics
                            elif basis_option_count == 1:
                                y = float( namd_data_line[ 34 ] ) - float( namd_data_line[4] )
                                if y < 0:
                                    y = 0.00

                        else:             

                            # mb dynamics
                            if basis_option_count == 0:
                                y = float( namd_data_line[ 35 ] ) - float( namd_data_line[4] )
                                if y < 0:
                                    y = 0.00
 
                            # mixed_sd dynamics
                            elif basis_option_count == 1:
                                #y = float( namd_data_line[ 71 ] ) - float( namd_data_line[4] )
                                y = float( namd_data_line[ 35 ] ) - float( namd_data_line[4] )
                                if y < 0:
                                    y = 0.00

                        y = y * units.au2ev

                    xdata.append(i*dt)
                    ydata.append(y)

                dyns[decoherence_option_count][basis_option_count].append( np.array( ydata ) )

        dyns[decoherence_option_count][basis_option_count] = np.array( dyns[decoherence_option_count][basis_option_count] )

        basis_option_count += 1

    dyns[decoherence_option_count] = np.array( dyns[decoherence_option_count] )

    ####################
    # 4.0 We have compute the dynamics for all bases (mb, mixed_sd, elec_sd, hole_sd) for a given decoherence scheme. Now, let's make our plots
    ydata_avg_mb       = sum(dyns[decoherence_option_count][0]) / len(dyns[decoherence_option_count][0])
    ydata_avg_mixed_sd = sum(dyns[decoherence_option_count][1]) / len(dyns[decoherence_option_count][1])

    dx   = 0.01
    emin = 0.0
    emax = 2.0
    plot_pop_dens( dyns[decoherence_option_count][0], ydata_avg_mb, dt, dx, emin, emax, decoherence_options_names[decoherence_option_count], "mb", printout=0)
    plot_pop_dens( dyns[decoherence_option_count][1], ydata_avg_mixed_sd, dt, dx, emin, emax, decoherence_options_names[decoherence_option_count], "mixed_sd", printout=0)

    if fit_option == 0:
        ydata_avg_mb_fit, tau_avg_mb, beta_mb, E0, r_squared_mb                         = fit_data( xdata, ydata_avg_mb, fit_option )
        ydata_avg_mixed_sd_fit, tau_avg_mixed_sd, beta_mixed_sd, E0, r_squared_mixed_sd = fit_data( xdata, ydata_avg_mixed_sd, fit_option )
        print ("\nFinished fitting average dynamics, stretch compressed", decoherence_options[decoherence_option_count], [basis_option_count])
        print ("tau_avg_mb", decoherence_options[decoherence_option_count], tau_avg_mb, " fs")
        print ("beta mb   ", beta_mb)
        print ("R2 mb     ", r_squared_mb)
        print ("\n")
        print ("tau_avg_sd", decoherence_options[decoherence_option_count], tau_avg_mixed_sd, " fs")
        print ("beta sd   ", beta_mixed_sd)
        print ("R2 sd     ", r_squared_mixed_sd)
         
    elif fit_option == 1:
        ydata_avg_mb_fit, tau_avg_mb, E0, r_squared_mb                   = fit_data( xdata, ydata_avg_mb, fit_option )
        ydata_avg_mixed_sd_fit, tau_avg_mixed_sd, E0, r_squared_mixed_sd = fit_data( xdata, ydata_avg_mixed_sd, fit_option )
        print ("\nFinished fitting average dynamics, single exponential", decoherence_options[decoherence_option_count], [basis_option_count])
        print ("tau_avg_mb", decoherence_options[decoherence_option_count], tau_avg_mb, " fs")
        print ("R2 mb     ", r_squared_mb)
        print ("\n")
        print ("tau_avg_sd", decoherence_options[decoherence_option_count], tau_avg_mixed_sd, " fs")
        print ("R2 sd     ", r_squared_mixed_sd)

    if fit_option == 3:
        ydata_avg_mb_fit, a_mb, tau1_avg_mb, tau2_avg_mb, E0_mb, r_squared_mb                               = fit_data( xdata, ydata_avg_mb, fit_option )
        ydata_avg_mixed_sd_fit, a_mixed_sd, tau1_avg_mixed_sd, tau2_avg_mixed_sd, E0_mixed_sd, r_squared_mixed_sd = fit_data( xdata, ydata_avg_mixed_sd, fit_option )
        print ("\nFinished fitting average dynamics, gaussian exponential", decoherence_options[decoherence_option_count], [basis_option_count])
        print ("tau1_avg_mb ", decoherence_options[decoherence_option_count], tau1_avg_mb, " fs")
        print ("tau2_avg_mb ", decoherence_options[decoherence_option_count], tau2_avg_mb, " fs")
        print ("a mb        ", a_mb)
        print ("E0 mb       ", E0_mb)
        print ("R2 mb       ", r_squared_mb)
        timescale_mb = ( a_mb / E0_mb ) * tau1_avg_mb + ( ( E0_mb - a_mb ) / E0_mb ) * tau2_avg_mb
        print ("computed timescale mb ", decoherence_options[decoherence_option_count], timescale_mb )
        print ("computed rate mb", decoherence_options[decoherence_option_count], 1.0 / timescale_mb )

        print ("\n")
        print ("tau1_avg_mixed_sd ", decoherence_options[decoherence_option_count], tau1_avg_mixed_sd, " fs")
        print ("tau2_avg_mixed_sd ", decoherence_options[decoherence_option_count], tau2_avg_mixed_sd, " fs")
        print ("a mixed_sd        ", a_mixed_sd)
        print ("E0 mixed_sd       ", E0_mixed_sd)
        print ("R2 mixed_sd       ", r_squared_mixed_sd)
        timescale_mixed_sd = ( a_mixed_sd / E0_mixed_sd ) * tau1_avg_mixed_sd + ( ( E0_mixed_sd - a_mixed_sd ) / E0_mixed_sd ) * tau2_avg_mixed_sd
        print ("computed timescale mixed_sd ", decoherence_options[decoherence_option_count], timescale_mixed_sd )
        print ("computed rate mixed_sd", decoherence_options[decoherence_option_count], 1.0 / timescale_mixed_sd )

    decoherence_option_count += 1

