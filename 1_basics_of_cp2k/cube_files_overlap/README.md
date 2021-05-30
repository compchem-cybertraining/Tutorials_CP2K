# 1. Cutoff analysis for cube files overlaps

The last part in this manual is about the effect of CP2K parameters on the overlap matrices. CP2K uses the converged wavefunction to produce the cube files for each molecular 
orbital. It will append the numerical values of the wavefunction computed at each grid point in a cube file. The grid points is dependent on the cutoff value. 
The larger the cutoff the higher the number of grid points are. In fact, the cutoff value "defines" the number of grid points. 

Now, if we have a huge cutoff value, like 1000 Ry, and we want to append all the values of the plane waves in all the grid points, then each cube file size will rise up to
hundreds of MBs or even some GBs. This is also dependent on the system cell size as well and the larger the cell size the larger the cube files will be. There are different 
ways to avoid high computational cost for outputting the cube files by CP2K. Specially, if we are going to use a large number of states, for example 30 states in valence 
band and 20 in conduction band. Each cube file is 600 MB, for example, and therefore all the cube files will be up to 30 GB for only one time step and for two consecutive 
time step this would be 60 GB that in some cases may cause overflow of the disk.

One way to overcome this problem and avoiding the overflow of the disk space is to use the keyword `STRIDE` in the input file. The `STRIDE` keyword will use and evaluate 
the wave functions only at some grid points. To clarify what it means, it is better to give an example. The `STRIDE 1 1 1` will use all the grid points in all the three X, Y, 
and Z axis. The `STRIDE 2 2 2` will use half of the grid points on each of the X, Y, and Z axis and `STRIDE 3 3 3` will use one third of the grid points on each of the axes and 
so on. The choice of the `STRIDE` needs to be studied for each system separately. For example if one chooses a huge cutoff value, say 1500 Ry, the cube files will be massively 
large and the computational times will increase. It is worth noting that Libra uses multiprocessing for reading the cube files and it is fast enough and may not take a very long 
time but note that one should also be very careful about the memory of the computer. If it overflows the code will crash. Therefore, the use of `STRIDE` is necessary. Now, if we 
want to check which `STRIDE` we should use we can start with a few number of Kohn-Sham states, say only 5 in valence band and 5 in conduction band, and start from `STRIDE 1 1 1`, 
and save the overlap matrices obtained from this parameter. Then we can move to `STRIDE 2 2 2` and compare the computed overlaps with the ones obtained from the finer mesh points 
(`STRIDE 1 1 1`). This loop is continued until we get to a large difference between the overlap values that can change the matrix elements significantly. The latter is the optimum 
`STRIDE` to use. Note that if one uses `STRIDE n` it is equivalent to `STRIDE n n n` in which `n` is an integer value.

Here, we have adopted the `step2` files from [this repository](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP). What we need to do is to check the results for different 
cutoff values or time step. We have provided Python file, `cutoff_analysis.py`, which will change the `cp2k_input_template.inp` cutoff and 
runs the calculations by `python run.py`. This is done through `for new_cutoff in range(300,1401,100):` line. It will consider the cutoff values of 300,400,500,..,1400 Ry. 
Then, it will change the directory name to `res-cutoffvalue` and you can see the overlaps and time-overlaps for the above purpose. 
There are two equilibrated trajectories available. One is with 1 fs time step, `MD_BA2PbI4_dt1fs-pos-1.xyz` and the other with 0.1 fs time step, `MD_BA2PbI4_dt0.1fs-pos-1.xyz`.
To see the effect of time step on time overlaps, you can comment the current `trajectory_xyz_file` with 1 fs time step and uncomment the other one which is with 0.1 fs time 
step. Then you can run the calculations for different cutoff values as was mentioned above.

**Note**: In `run.py`, we have can use `os.system("sbatch submit_"+str(njob)+".slm")` which submits the job on CCR. But if one wants to run this file without submitting it can be changed 
to `os.system("sh submit_"+str(njob)+".slm")`.

This can also be done for different `STRIDE` values as well. You can either move the computed overlaps to a new directory and change the `STRIDE` in `cp2k_input_template.inp`, 
or you can write a new Python code that works based on changing the `STRIDE` keyword using the above Python code as an example. You can also set this up in two
for loops, one for `CUTOFF` and one for `STRIDE` as well.

One can see the changes of these values or any other parameter that is in CP2K input (such as `&VDW_POTENTIAL` for the effect of dispersion correction just as an example).
You can check the results of the `step3` and `step4` to see the changes on nonadiabatic couplings or nonadiabatic dynamics as well.
