# Step2 - Compute the overlap matrices and energies in the Kohn-Sham basis and perform TD-DFT calculations along the precomputed nuclear trajectories

## 1. CP2K input for electronic structure calculations

`cp2k_input_template.inp`

This file is a template of a cp2k input file and is used to compute the electronic structure of the system along the precomputed nuclear trajectory. SCF calculations are called with the parameter `RUN_TYPE ENERGY`. 

In the section FORCE_EVAL < PROPERTIES, is the input needed to compute TD-DFT calculations. If one wishes to compute only KS properties (overlaps, energies) one can delete this section. One would also need to set completion_level = 0 in the file submit_template.slm. One can set completion_level = 0 and also perform TDDFT calculations as well. In this case, the computation of the MB basis would be done in a post-process style. This is what we do here. 

If the user needs to change the cp2k input, the changes should be done to this input file. One must make sure that the cube files can be produced via WRITE_CUBE .TRUE. 

Other required files for running the CP2K input file are basis set and pseudopotential files or any other files required to run the calculations, such as `dftd3.dat`. The full path to these files in the `cp2k_input_template.inp` file shoud be specified.

[CP2K paper](https://aip.scitation.org/doi/pdf/10.1063/5.0007045)

[Difference between TD-DFT and TD-DFPT](https://groups.google.com/g/cp2k/c/xj8udnSyeEI)

The keyword `RESTART` increases the speed of calculations, both for SCF and TD-DFPT calculations. Therefore, the `RESTART` is required to be set to `.TRUE.` in `TDDFPT` section. Also, the `WFN_RESTART_FILE_NAME` in this section should exist with a random `tdwfn` file name. The same is also needed in the `FORCE_EVAL` section. `WFN_RESTRAT_FILE_NAME` should exist in the input with an random `wfn` file name. 

In the `&MO_CUBES` section the number of occupied and unoccupied orbitals must be specified. This is dependent on the TD-DFPT calculations and the user have to make a good guess to make sure that the cube files of all the states in the excitation analysis of TD-DFPT calculations exist. This guess can be obtained from running the calculations for 5-10 steps.

## 2. Bash file for running the calculations for one job

The standard sample bash file for submitting the calculations and running the Python code through `slurm` and `sbatch` is the `submit_template.slm`. First one should load all needed modules including the modules required for loading CP2K. This file contains the input variables required for calculations of the overlap matrices and NACs. We list the variables as follows:

`nprocs`: The number of processors used for calculations. Note that the same number of processors should be specified in the `#SBATCH --ntasks-per-node` above.

`cp2k_exe`: The executable CP2K or the full path to executable CP2K folder.

`res`: The directory for storing the overlap matrices and NACs.

`min_band`: The minimum KS orbital index to be considered.

`max_band`: The maximum KS orbital index to be considered.

`ks_orbital_homo_index`: The HOMO index for KS orbitals.

`MO_images_directory`: The directory where the molecular orbital isosurfaces images are stored. Not used in the current computations

`path_to_tcl_file`: The path to the `tcl` file for plotting the cube files. This file is fed to VMD for plotting the cube files. The file can be different according to the user need. Here we have uploaded a sample file for our purpose as `cube.tcl`.

`states_to_be_plotted`: The index of the states to be plotted by VMD. If there are many states considered for plotting they should be separated by comma.

`job_init_step`, `nsteps_this_job`, `njob`: Please leave these variables as they are. The program will automatically recognize and fill them.

Now that we have set some of the variables we need to run the Python code as `python -c`. First we need to import the function `step2_many_body` from `libra_py.workflows.nbra`. Then we have to set up some variables as a dictionary in `params`. The variables in `params` are as follows:

`es_software`: This defines the electronic structure calculations software name. Please, insert the software name in lower case such as `cp2k`.

`es_software_input_template`: The input template for the `es_software`. We recommend the use of the prepared input template as is mentioned above.

`es_software_exe`: The executable `es_software`. This is defined above as `$cp2k_exe`. 

`nprocs`: Number of processors to be used.

`project_name`: The project name. The cube file names is based on this variable.

`trajectory_xyz_filename`: The name of the trajectory `xyz` file. 

`logfile_directory`: The directory for CP2K output log files. 

`number_of_states`: The number of excited states to be considered. This value should not be larger than `NSTATE` in `TDDFPT` section.

`tolerance`: The tolerance factor will consider the excited states which the square of their configuration interaction coefficient is larger than this value. Set to zero here

`isUKS`: If this flag is set to **`1`**, then the program will consider the unrestricted spin calculations. If it is set to other values it will consider only the spin restricted case. Please, make sure for the `UKS .TRUE.` in the `cp2k_input_template.inp` if you have set this variable to **`1`**.

`min_band`: The minimum state index of the KS states. This will take the above value as `$min_band`.

`max_band`: The maximum state index of the KS states. This will take the above value as `$max_band`.

`ks_orbital_homo_index`: The index of the HOMO energy level in the KS basis. This will take the above value as `$ks_orbital_homo_index`.

`completion_level`: How much of the calculations to compute. 0 - compute KS overlaps and TDDFT calculations 1 - also compute Hvib in the MB basis. Here, we set this to zero, and create the Hvib in the MB basis in a post-process fashion. To compute only the overlaps and energies in the KS basis and not compute the TDDFT calculations, one has to turn off the option to compute TDDFT calculations in the cp2k input template

`do_phase_corrections`: The flag to perform phase correction. If this value is set to **`1`**, the program will apply the [phase correction algorithm](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b02826) to the overlap matrices.

`perform_state_reordering`: This flag will do the state tracking (reordering algorithm) if its value is set to **`1`**. There are different schemes of state reordering which can be defined by the variable `do_state_reordering`.

`state_reordering_alpha`: 

`istep`: The initial time step of the job to be considered. This value is set as `$job_init_step`. The program will start from the `istep` geometry in the trajectory `xyz` file.

`nsteps_this_job`: The number of steps for the job. This value is set as `$nsteps_this_job`.

`njob`: The job number which is defined as `$njob`.

`res_dir`: The `res` directory for storing the energies, overlap matrices and NACs. It is defined as `$res`.

`dt`: The time step used in the MD in **_atomic units_**. 

`do_cube_visualization`: If this flag is set to **`1`** the program will call VMD to plot the cube files. If you set this variable to **`1`**, pleae make sure to load VMD through `module load` or add the executable path to the `PATH` variable like `export PATH=/full/path/to/vmd/folder:$PATH`.

`path_to_tcl_file`: The full path to `tcl` file as the input of VMD. This variable is defined with `$path_to_tcl_file` as above.

`states_to_be_plotted`: The states to be plotted which is defined as `$states_to_be_plotted`.

`MO_images_directory`: The path to store the images of the molecular orbitals plotted by VMD. It is defined as `$MO_images_directory` as above.

Finally, we run the calculations using the `step2_many_body.run_step2_many_body( params )` for the parameters we have set into `params` dictionary variable.

## 3. Run all the jobs through `run.py`

`run.py` file is a file which has the burden to split the trajectory into multiple smaller trajectories and then submit the `slurm` or `pbs` bash files. It starts by making the `wd` folder which is the working directory. Then it will extract the `xyz` coordinates required for each step of each job. The code will also prepare the input files for each ste for each job based on the `cp2k_input_template.inp` file. This is done in each `job` folder in the `wd`.

The required inputs in the `run.py` file are as follows:

`trajectory_xyz_file`: The user must specify the full path of the trajectory `xyz` file or the name of the file if it is in the current folder. 

`es_software_input_template`: The input template for the electronic structure calculations. Here, the input template is `cp2k_input_template.inp` file.

`es_software`: The name of the electronic structure calculations software e.g. `cp2k`.

`istep`: The initial step from the _trajectory `xyz` file_. 

`fstep`: The final step from the _trajectory `xyz` file_. 

`njobs`: The total number of jobs that the user wants to submit. The maximum number of jobs must be less than half of the total steps i.e. `fstep-istep+1`.

`os.system("sbatch submit_"+str(njob)+".slm")`: The jobs are submitted through this line of code at the end of the `run.py`. Please, change it according to your HPC submission platform. For example if you use `pbs` files and you use `qsub`, after preparing the `submit_template.pbs` the same as `submit_template.slm` you can change this line to `os.system("qsub submit_"+str(njob)+".pbs")`.

**_Note_:** You can submit all your jobs by running only `python run.py` and the submission process of the jobs will be done on the local node. An alternative way for submitting the jobs is through submitting the `submit.slm` file which contains `python run.py`. This can be done if other nodes have the capability to perform the submission. Unless you have to use only `python run.py` on the local node or any other node that has the capability for submitting the jobs. Please also note that the current set up submits 800 jobs. To test the workflow it is advised to first set the following variables:
istep = 0
fstep = 1
njobs = 1
This way, the user can test the software for the first two steps to debug for potential errors in input, etc.  

**_NOTE_:** - Please note that the paths currently defined in these files may not be the correct paths for you. Please adjust all paths to your specific needs
