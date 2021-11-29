import os
import sys
from libra_py import CP2K_methods

run_slurm = False
submit_template = 'submit_template.slm'
run_python_file = 'run_template.py'
istep = 1
fstep = 5
njobs = 1
os.system('rm -rf res job* all_logfiles all_pdosfiles')

print('Distributing jobs...')
CP2K_methods.distribute_cp2k_libint_jobs(submit_template, run_python_file, istep, fstep, njobs, run_slurm)

