# Computing the MO overlaps using Libra/Libint interface

## Running the calculations

The example in this repository, shows the computation of the MO overlaps for a Li@C60 molecular system in DFT calculations. 
The number of electrons for this system is odd and therefore we need
to consider the unrestricted spin calculations in the CP2K input files, both for MD and single-point calculations. This is done by adding the `UKS .TRUE.` to CP2K input files.

Here are the instructions to run the step2 for computing the MO overlaps for this system. 

**NOTE:** You can use this file for extended tight-binding calculations as well.

First, you need to modify the `run_template.py` file. Here is the description for each 
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


After you applied the changes to the `run_template.py`, you need to modify the `distribute_jobs.py` file. The variables in this file are as follows:


`run_slurm`: The flag to submit the job using slurm environment by `sbatch`. 

`submit_template`: The submit template file name. You can also add the full path to this file as well.

`run_python_file`: The python file to run the calculations.

`istep`: The first MD step in the trajectory `xyz` file.

`fstep`: The final MD step in the trajectory `xyz` file.

`njobs`: The number of jobs to be submitted if you are using `sbatch`. The code will distribute the jobs for you using the `CP2K_methods.distribute_cp2k_libint_jobs` function.

Now, the only thing that is left is to run the calculations by `python distribute_jobs.py`. This will start computing the MO overlaps. 



## Check the results

The energies and MO overlaps are stored in sparse format to reduce the disk space usage. We use the `scipy.sparse` library to save and load the sparse files which 
have the `npz` format. As an example, you can see the results in the `res` directory using the following commands in environments such as `python`, `ipython`, or Jupyter 
files:

```
In [1]:

import numpy as np
# import the scipy.sparse library
import scipy.sparse as sp
# load the sparse matrix
S_sparse = sp.load_npz('res/S_ks_1.npz')
# transform to dense format (numpy ndarray)
S_dense = S_sparse.todense()
# print the diagonal elements of the S matrix (the MO overlaps of the same molecular configuration)
# to check the orthonormality of the MOs
print('The diagonal elements of the S matrix are:\n')
print(np.diag(S_dense))

Out [1]:

The diagonal elements of the S matrix are:

[0.9999995 +0.j 0.99999952+0.j 0.99999923+0.j 0.99999944+0.j
 0.99999939+0.j 0.99999971+0.j 0.99999929+0.j 0.99999933+0.j
 1.00000141+0.j 0.99999871+0.j 0.99999868+0.j 0.99999885+0.j
 0.99999882+0.j 0.99999887+0.j 0.99999914+0.j 1.00000061+0.j
 1.00000095+0.j 0.99999843+0.j 0.99999951+0.j 0.99999926+0.j
 0.99999923+0.j 0.99999943+0.j 1.00000224+0.j 1.00000019+0.j
 1.00000017+0.j 1.00000145+0.j 0.99999868+0.j 0.99999871+0.j
 0.99999867+0.j 0.99999882+0.j 0.99999885+0.j 1.00000092+0.j
 0.9999993 +0.j 1.0000007 +0.j 1.00000099+0.j 0.99999839+0.j]

```

This can be used for all types of `npz` files, for example the energies.



