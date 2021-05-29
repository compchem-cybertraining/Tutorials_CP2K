## This is inplace of the notebook for now -- to be converted to a notebook later

from liblibra_core import *

import multiprocessing as mp
import numpy as np
import math
import sys
import os


from libra_py import CP2K_methods
from libra_py import Gaussian_methods
from libra_py.workflows.nbra import step2_many_body

########################### creating directories ###############################
try:
    os.system("rm -r wd")
    os.mkdir('wd')
    os.system("rm -r res")
    os.mkdir('res')
except:
    sys.exit(0)

######################### creating and submitting jobs #########################
# For MD with time step of 1 fs
trajectory_xyz_file        = "MD_BA2PbI4_dt1fs-pos-1.xyz"

# For MD with time step of 0.1 fs
#trajectory_xyz_file        = "MD_BA2PbI4_dt0.1fs-pos-1.xyz"

es_software_input_template = "cp2k_input_template.inp"
es_software = "cp2k"
#es_software_input_template = "gaussian_input_template.gjf"
#es_software = "gaussian"
istep = 0
fstep = 1#3000 
njobs = 1#600
for njob in range( njobs ):

    job_init_step, job_final_step = step2_many_body.curr_and_final_step_job( istep, fstep, njobs, njob )
    nsteps_this_job = job_final_step - job_init_step + 1

    print( '\n1- We are before creating job ', njob, ' in directory:', os.getcwd() )
    # Here we create the jobs through auxiliary_functions.cp2k_distribute
    if es_software.lower() == "cp2k":
        CP2K_methods.cp2k_distribute( job_init_step, job_final_step, nsteps_this_job, trajectory_xyz_file, es_software_input_template, njob )

    elif es_software.lower() == "gaussian":
        Gaussian_methods.gaussian_distribute( job_init_step, job_final_step, nsteps_this_job, trajectory_xyz_file, es_software_input_template, njob )

    print('2- Finished distributing for job ',njob,'\nNow we are in directory', os.getcwd())
    os.chdir("wd/job"+str(njob)+"/")

    print('3- Now we have changed directory to :', os.getcwd(),'To submit job', njob)
    os.mkdir("cubefiles")
    os.mkdir("logfiles")
    
    print("4- The initial step for the job ", njob, " is: ", job_init_step, 'with final step: ', job_final_step)
    print("5- nsteps_this_job is: ", nsteps_this_job)
    St_ks_job = []
    S_ks_job  = []
    debug     = []

    ############################
    os.system("cp ../../submit_template.slm submit_"+str(njob)+".slm")
    # Now, open the submit_template file in this job folder
    # Add the values to the params for this job
    f = open( "submit_"+str(njob)+".slm" )
    submit_template_file = f.readlines()
    submit_template_file_size = len( submit_template_file )
    f.close()

    f = open( "submit_"+str(njob)+".slm" , 'w' ); f.close()
    for i in range( submit_template_file_size ):      

        f = open( "submit_"+str(njob)+".slm" , 'a' )

        submit_template_file_line = submit_template_file[i].split()
        if not submit_template_file_line:
            continue

        elif submit_template_file_line[0] == "job_init_step=":
            f.write( "declare -i job_init_step=%i" % (job_init_step) )
            f.write("\n")

        elif submit_template_file_line[0] == "nsteps_this_job=":
            f.write( "declare -i nsteps_this_job=%i" % (nsteps_this_job) )
            f.write("\n")

        elif submit_template_file_line[0] == "njob=":
            f.write( "declare -i njob=%i" % (njob) )
            f.write("\n")

        else:
            f.write( submit_template_file[i] )            
 
    #os.system("sbatch submit_"+str(njob)+".slm")
    os.system("sh submit_"+str(njob)+".slm")

    print("6- Now submitting the job in folder:  ", os.getcwd())

    # Change directory to the main directory
    os.chdir("../../")
    print("7- Finished submitting job, now we are back in directory:", os.getcwd(),'\n\n\n')

print("\nDone!\n")


