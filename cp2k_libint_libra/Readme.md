# Computing the MO overlaps using Libra/Libint interface for DFT calculations

The example in this repository, shows the computation of the MO overlaps for a Li@C60 molecular system. The number of electrons for this system is odd and therefore we need
to consider the unrestricted spin calculations in the CP2K input files, both for MD and single-point calculations. This is done by adding the `UKS .TRUE.` to CP2K input files.

Here are the instructions to run the step2 for computing the MO overlaps for this system. First, you need to modify the `run_template.py` file. Here is the description for each 
variable in this file and the ones that you need to change for run:

`path`: The full path to your `step2` folder. 

`params['nprocs']`: Number of processors 

`params['istep']`: The initial step of the job. You do not need to fill this and Libra will fill that based on the number of jobs and number of steps.

`params['fstep']`: The final step of the job. Again, leave this and do not fill that.

`params['init_state']`: The minimum state number to be considered (`min_band` as in previous Libra/CP2K workflows using cube files).

`params['final_state']`: The maximum state number to be considered (`max_band` as in previous Libra/CP2K workflows using cube files).

`params['isxTB']`: The flag for extended tight-bingind (xTB) calculations. Note that if you set this flag to `True`, you will need to provide the inputs for both 
orbital transformation (OT) and diagonalization methods for xTB in `params['cp2k_ot_input_template']` and `params['cp2k_diag_input_template']`. 

`params['isUKS']`: The flag for unrestricted spin calculations. If set to `True`, then it will consider the spin-polarized case (with alpha and beta spins). Make sure you have already added the `UKS .TRUE.` to your CP2K input template file in `step2_compute_overlap` folder (Here is `pbe.inp`).

`params['is_spherical']`: The MO coefficients in spherical or Cartesian basis.

`params['remove_molden']`: This example uses only the `molden` files for computation of the MO overlaps and these can become quite large for very large systems. If this flag
is set to `True`, then it will remove this file after each step is done. We also recommned to set this flag to `True`.

`params['res_dir']`: The full path to save the MO overlaps and energies.

`params['all_pdosfiles']`: The full path to save all the `.pdos` files, produced by CP2K, in each step.

`params['all_logfiles']`: The full path to sav all the `.log` files of CP2K run.

`params['cp2k_exe']`: The full path to the CP2K executable or if you load it through `module load` or added that using `export PATH`, you can just set the executable name
such as `cp2k.psmp`.

`params['cp2k_ot_input_template']`: This is used for xTB calculations and it is just the CP2K input that computes the electronic structures using the OT method. If you do not 
use the xTB, please set it to an empty string like `''`.

`params['cp2k_diag_input_template']`: This is the CP2K input file that uses the diagonalization method to compute the electronic structure calculations. This is for both 
DFT and xTB calculations. 

`params['trajectory_xyz_filename']`: The full path to trajectory `.xyz` file.





