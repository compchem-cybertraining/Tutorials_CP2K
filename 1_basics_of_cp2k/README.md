# CP2K manual for electronic structure calculations

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

The procedure adopted here is as follows:

* 2.1. [Obtaining the `xyz` coordinates from a `cif` file](1_structure_preparation/README.md)
* 2.2. [Convergence analysis](2_convergence_analysis/README.md)
* 2.3. [Geometry optimization](3_geometry_optimization/README.md)
* 2.4. [Cell optimization](4_cell_optimization/README.md)
* 2.5. [TD-DFT calculations](5_tddft/README.md)
* 2.6. [Molecular dynamics](6_molecular_dynamics/README.md)
* 2.7. [Calculations with hybrid functionals](7_hybrid_functionals/README.md)
* The effect of CP2K parameters on the cube file (These files are used in [Libra](https://github.com/Quantum-Dynamics-Hub/libra-code) package.)
    * Cube file sizes
    * Molecular orbitals overlap 
    * Molecular orbitals time-overlap
        * The effect of time step 

We highly welcome improving the functionality of the input files in this repository. 
So, please feel free to share your inputs and timings with us if you used these inputs.


