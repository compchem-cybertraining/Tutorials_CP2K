import os
import sys
import time
import numpy as np
import multiprocessing as mp
import glob
import matplotlib.pyplot as plt



def convolve_pdos(cp2k_log_file: str, time_step: int, sigma: float, coef: float, npoints: int, energy_conversion: float, angular_momentum_cols: list):
    """
    This function reads the pdos file produced by CP2K and extract the pdos at each time step and 
    then convolve them with Gaussian functions.
    
    Args:
    
        cp2k_log_file (str): The CP2K output file.
        
        time_step (int): The time step of molecular dynamics.
        
        sigma (float): The standard deviation in Gaussian function.
        
        coef (float): The coefficient multiplied in Gaussian function.
        
        npoints (int): The number of points used in convolution.
        
        energy_conversion (float): The energy conversion unit from Hartree. For example 27.211386 is
                                   for unit conversion from Hartree to eV.
				   
        angular_momentum_cols (list): The angular momentum columns in the *.pdos files produced by CP2K.
	
    Returns:
	
        energy_grid (numpy array): The energy grid points vector.
		
        convolved_pdos (numpy array): The convolved pDOS vector.
		
        homo_energy (float): The average HOMO energy.
		
    """

    # Opening the file
    file = open(cp2k_log_file,'r')
    lines = file.readlines()
    file.close()
    
    # Lines with 'DOS'
    lines_with_dos = []
    
    # Finding the lines with 'DOS'
    for i in range(0,len(lines)):
        if 'DOS'.lower() in lines[i].lower().split():
            lines_with_dos.append(i)
    
    # Finding the first and last index of PDOS for each time step
    if len(lines_with_dos)==1:
        # First index
        first_index = 2
        # Last index
        last_index = int(lines[len(lines)-1].split()[0])
    elif len(lines_with_dos)>1:
        # First index
        first_index = 2
        # Last index
        last_index = int(lines_with_dos[1]-1)
    
    # Find the number of columns in the PDOS file showing the number 
    # of orbital components, energy, and occupation column.
    num_cols = len(lines[first_index].split())
    
    # Number of energy levels considered for PDOS
    num_levels = last_index - first_index + 1

    
    # Finding the homo and lumo energy level by appending the 
    # pdos numerical values of unoccupied states only
    pdos_unocc = []
    # Energy levels
    energy_levels = []
    for i in range(first_index, last_index + 1):
        energy_levels.append(float(lines[i].split()[1])*energy_conversion)
        if float(lines[i].split()[2])==0:
            pdos_unocc.append(i)
    # LUMO energy level
    lumo_level = int(lines[min(pdos_unocc)].split()[0])
    # HOMO energy level
    homo_level = lumo_level-1
    # HOMO energy
    homo_energy = float(lines[homo_level].split()[1])*energy_conversion
    # Minimum energy level
    min_energy = float(lines[first_index].split()[1])*energy_conversion
    # Maximum energy level
    max_energy = float(lines[last_index].split()[1])*energy_conversion
    
    
    # Now we make an equispaced energy vector from min_energy ad max_energy with npoints.
    energy_grid = np.linspace( min_energy-2, max_energy+2, npoints )
    energy_grid = np.array(energy_grid)
    
    
    # Appending the energy lines with their component densities of states
    energy_lines = []
    for i in range( time_step * ( num_levels + 2 ) + 2, ( time_step + 1 ) * ( num_levels + 2 ) ):
        # Appending the energy lines into enrgy_lines
        energy_lines.append( lines[i].split() )

    for i in range(0, len(energy_lines)):
        
        for j in range(0,len(energy_lines[0])):
            
            energy_lines[i][j] = float(energy_lines[i][j])
            
    energy_lines = np.array(energy_lines)

    # Now we sum the PDOSs defined in angular_momentum_cols by user
    pdos_sum = []
    for k in range(0, len(energy_lines)):
        
        # A temporary vector for summation of the PDOS
        tmp_vec = []
        tmp_vec.append(energy_lines[k][1])
        
        for i in range(0,len(angular_momentum_cols)):
            # Initializing a new sum variable
            # print("angular_momentum_cols[i]",angular_momentum_cols[i])
            tmp_sum = 0
            for j in angular_momentum_cols[i]:
                
                # If j is less than the number of columns 
                # then sum the PDOS
                if j<=num_cols:
                    tmp_sum += energy_lines[k][j]
            # Appending tmp_sum into tmp_vec
            tmp_vec.append(tmp_sum)
        
        # Now append tmp_vec into pdos_sum, we will
        # then use this pdos_sum for convolution
        pdos_sum.append(tmp_vec)
    
    convolved_pdos = []
    t1 = time.time()
    for j in range(1,len(angular_momentum_cols)+1):
        # Initialize a vector of zeros summing the weighted PDOS
        tmp_weighted_pdos = np.zeros(energy_grid.shape)

        for i in range(0,num_levels):
            # The Guassian function
            gaussian_fun = (coef/(sigma*np.sqrt(2.0*np.pi)))*(np.exp(-0.5*np.power(((energy_grid-float(pdos_sum[i][0])*energy_conversion)/sigma),2)))
            
            tmp_weighted_pdos = tmp_weighted_pdos + gaussian_fun * float( pdos_sum[i][j] )
        convolved_pdos.append(tmp_weighted_pdos)
    print('Elapsed time for convolving ',cp2k_log_file,': ',time.time()-t1,' seconds')
    convolved_pdos = np.array(convolved_pdos)
    
    energy_grid = energy_grid
    
    return energy_grid, convolved_pdos, homo_energy









####################################################################################################################################################
#============================= Main part starts from here and we use the function above combined with multiprocessing =============================#

#pdos_type = "orbital_resolved"
pdos_type = "atom_resolved"
#thermal=True
thermal=False

if pdos_type == "orbital_resolved":

    #===================== Orbital resolved columns
    #                               C, total            C, s               C, p                C, d
    angular_momentum_cols = [ [ list(range(3,12)), list(range(3,4)), list(range(4,7)), list(range(7,12)) ], \
    #                               H, total            H, s               H, p        
                              [ list(range(3,7)), list(range(3,4)), list(range(4,7)) ] ]
    # This is for orbital resolved
    labels = ['C, s','C, p','C, d','H, s','H, p']
    # Colors, the color orders are based on the labels, Here are chosen for elements
    colors = ['red','blue','green','orange', 'purple']
    outname = 'orbital'

elif pdos_type == "atom_resolved":

    #===================== Atom resolved columns
    #                           C, total          C, total
    angular_momentum_cols = [ [ list(range(3,12)), list(range(3,12))],\
    #                           H, total          H, total
                              [ list(range(3,7)), list(range(3,7)) ]  ]
    # This is for atom resolved
    labels = ['C','H']
    # This row is for atom resolved
    colors = ['blue','red']
    outname = 'atom'





#===================== Other inputs
# The time step is MD, For static calculations we only use 0
time_step = 0
sigma = 0.1
coef = 1
npoints = 2000 # number of points for the grid mesh vector
energy_conversion = 27.211386 # Hartree to eV
nprocs = 24 # number of processors

t1 = time.time()
# create the pool of processors
pool = mp.Pool(nprocs)

# This variable will contain the convolved PDOS for each element and each angular momentum list
Total_results = []
homos = []
#========================== Convolution
# Here we change the order,
# In fact the k1 is C, k2 is H
for k in [1,2]:
    # Create an empty list for summation of the convolved PDOS
    # This is used for Total DOS
    dos_summation_angular = []
    energy_grid_ave = []
    
    # Define a zero vector to sum the convolved PDOS for each angular momentum list
    # It is based on the number of points
    zero_vec = np.zeros((npoints))
    # Now for each angular momentum append the zero vector
    # The same for enegy grid average since we want the energy grid average
    for i in range(len(angular_momentum_cols[k-1])):
        dos_summation_angular.append(zero_vec)
        energy_grid_ave.append(zero_vec)

    # Create the variables for pool.starmap
    vars_for_pool = []
    # Find all the pdos files of an element in all_pdosfiles for the first trajectory

    if thermal == True:
        DOS_files1 = glob.glob('all_pdosfiles/*k%d-1.pdos'%k)
    else:
        DOS_files1 = glob.glob('../tddft/*k%d-1.pdos'%k)

    # For each of the pdos files we create the list of variables
    for DOS_file in DOS_files1:
        vars_for_pool.append((DOS_file,time_step,sigma,coef,npoints,energy_conversion,list(angular_momentum_cols[k-1])))       
        
    # The results of convolve_pdos function for an element
    results_for_k = pool.starmap(convolve_pdos,vars_for_pool)
    # We initialize all the homos average: homos_ave
    # This variable is the same for each element so it will be repeated but doesn't 
    # change the results
    homos_ave = 0
    # Now we need to take the average of them (The same is for average HOMO eergy level
    # and the energy grid as well)
    # Energy grid is the 0th element, convolved DOS is the 1st element, HOMO energy levels is the 2nd one
    # 1. We add the convolved ones to dos_summation_angular
    # Here we start to take the averages by summing the results.
    for i in range(len(DOS_files1)):
        energy_grid_ave       += results_for_k[i][0] # First element is the energy grid output by the function
        dos_summation_angular += results_for_k[i][1] # Second element is the convolved DOS
        homos_ave             += results_for_k[i][2] # Third element is the HOMO energy level
        
        
    # 2. We take the average for that by dividing by len(DOS_files)
    # 3. We append it to Total_results
    Total_results.append(dos_summation_angular/len(DOS_files1))
    # Uncomment only if you use two trajectories
    #Total_results.append(dos_summation_angular/(len(DOS_files1)+len(DOS_files2)))
    energy_grid_ave /= len(DOS_files1)
    # Uncomment only if you use two trajectories
    #energy_grid_ave /= (len(DOS_files1)+len(DOS_files2))
    homos_ave /= len(DOS_files1)
    # Uncomment only if you use two trajectories
    #homos_ave /= (len(DOS_files1)+len(DOS_files2))
    

# Close the pool
pool.close()
pool.join()
#================================== End of convolution



#================================== Total density
# We first compute the total density of states through the first computed angular momentum 
# for [3,12] as is shown above - So we only choose the 0 index
# Make a zero vector of npoints
total_density = np.zeros((npoints))
# i here is each element
for i in range(len(Total_results)):
    # Sum the total by Total_results of zero index angular momentum column for each element
    total_density += Total_results[i][0]

#============================== Plotting 
figure = plt.figure(num=None, figsize=(3.21, 2.41), dpi=1200, edgecolor='black', frameon=True)


# Plot the total density by black color
plt.plot(energy_grid_ave[0]-homos_ave,total_density,label='Total',color='black', linewidth=2.0)


# set up a counter for labels
c = 0
for i in range(len(Total_results)):
    for j in range(1,len(Total_results[i])):
        plt.plot(energy_grid_ave[0]-homos_ave,Total_results[i][j],label=labels[c],color=colors[c], linewidth=2)
        c += 1



plt.xlim(-7,15)
#plt.ylim(0,60)
plt.legend(fontsize=6.75, ncol=1, loc='upper center')
plt.xlabel('Energy, eV',fontsize=12)
plt.ylabel('DOS, 1/eV',fontsize=12)

if thermal == True:
    plt.title("Adamantane, 300 K",fontsize=12)
else:
    plt.title("Adamantane, 0 K",fontsize=12)

plt.tight_layout()
plt.show()

if thermal == True:
    if pdos_type == "orbital_resolved":
        plt.savefig('orbital_DOS_300K_'+outname+'_average.png', dpi=300)
    elif pdos_type == "atom_resolved":
        plt.savefig('atom_DOS_300K_'+outname+'_average.png', dpi=300)

elif thermal == False:
    if pdos_type == "orbital_resolved":
        plt.savefig('orbital_DOS_0K_'+outname+'.png', dpi=300)
    elif pdos_type == "atom_resolved":
        plt.savefig('atom_DOS_0K_'+outname+'.png', dpi=300)

print("Total Elapsed time:",time.time()-t1)

