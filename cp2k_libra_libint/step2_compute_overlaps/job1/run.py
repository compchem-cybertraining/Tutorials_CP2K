import os
import sys
from libra_py.workflows.nbra import step2


path = '/home/97425008/compchem/Tutorials_CP2K/cp2k_libint_libra/step2_compue_overlaps'
params = {}
params['nprocs'] = 25
params['istep'] = 1
params['fstep'] = 5
params['init_state'] = 122-9
params['final_state'] = 122+9
params['isxTB'] = False
params['isUKS'] = True
params['is_spherical'] =  True
params['remove_molden'] = True
params['res_dir'] = path + '/res'
params['all_pdosfiles'] = path + '/all_pdosfiles'
params['all_logfiles'] = path + '/all_logfiles'
params['cp2k_exe'] = '/home/97425008/cp2k-v7/cp2k/exe/local/cp2k.popt'
params['cp2k_ot_input_template'] = ''
params['cp2k_diag_input_template'] = path + '/pbe.inp'
params['trajectory_xyz_filename'] = path + '/lic-pbe-pos-1.xyz'

step2.run_cp2k_libint_step2(params)

