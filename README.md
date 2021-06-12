# CP2K tutorials for electronic structure calculations

This repository contains input examples for calculations associated with CP2k package, including 
the NonAdiabatic Molecular Dynamics (NA-MD) calculations via CP2K/Libra interface. 


## 1. General instructions

This repository summarizes inputs for different types of electronic structure calculations. 
The same files are also available in other places as well 
([project_cp2k_libra](https://github.com/AkimovLab/Project_Libra_CP2K), 
[project_perovskite_crystal_symmetry](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP)) for use in
nonadiabatic dynamics. However, here we provide detailed instructions on how to use CP2K and how the 
functionality and timings change with different inputs.

Here we use CP2K v6.1 compiled with Intel parallel studio 2019. For TD-DFT with hybrid functionals, 
we will use CP2K v7.1 which is compiled with GCC-8.3 compiler. This is because 
in lower versions of CP2K the TD-DFT calculation does not converge for hybrid functional due to 
a [problem](https://groups.google.com/g/cp2k/c/SEglKzKlVLQ/m/MyTavEqYBQAJ) with 
converging ADMM calculations. It is also worth noting that the computed results are done using Intel(R) Xeon(R) E5-2699 v4 @ 2.20GHz CPUs. 

For CP2K installation we recommend the use of `./install_too_chain.sh` in the `tools/toolchain` folder of CP2K. 
For compilation with Intel parallel studio you can use 
the instructions given in [XConfigure website](https://xconfigure.readthedocs.io/en/latest/cp2k/).


Here, we use a 2D perovskite of (BA)2PbI4 as a test and which its [cif](http://crystallography.net/cod/2102937.cif) 
file is available from the [crystallography website](http://crystallography.net/). We will use the unit 
cell with 156 atoms and a 2x2x2 supercell with 1248 atoms and check the performance of CP2K calculations functionality. 


## 2. Types of calculations 

* geometry preparation
* energy (single-point) calculations
* convergence with respect to basis set parameters
* geometry and cell optimization
* single-point calculations at the TD-DFT level
* single-point calculations with hybrid density functionals
* single-point calculations of huge systems
* molecular dynamics
* time-overlap calculations using cube files
* nonadiabatic coupling calculations 

The `legacy` folder contains a number of recent tutorials by Mohammad Shakiba and Brendan Smith, but 
they are yet to be organized and revised. These tutorials contain detailed and very useful instructions
which will eventually be migrated to the main folder here.



* [3. Extra resources]
 
  Here, we provide links to the repositories with the codes used in the following publications: 

  * [Nonadiabatic Dynamics in Si and CdSe Nanoclusters: Many-Body vs Single-Particle Treatment of Excited States](https://pubs.acs.org/doi/10.1021/acs.jctc.0c01009), [Technical details](https://github.com/AkimovLab/Project_Libra_CP2K)

  * [Crystal Symmetry and Static Electron Correlation Greatly Accelerate Nonradiative Dynamics in Lead Halide Perovskites](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.0c03799), [Technical details](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP)


We highly welcome improving the functionality of the input files in this repository. 
So, please feel free to share your inputs and timings with us if you used these inputs.


