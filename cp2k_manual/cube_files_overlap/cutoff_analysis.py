import os
import sys
import numpy as np

for new_cutoff in range(300,1401,100):

    # Open up the cp2k_input_template.inp file
    
    f1 = open('cp2k_input_template.inp','r')
    lines = f1.readlines()
    f1.close()
    
    # Finding the CUTOFF keyword line number
    cutoff_line_num = []
    for i in range(len(lines)):
        if 'cutoff' in lines[i].lower().split():
            cutoff_line_num.append(i)
    
    # Open up a new file which will replace the original one
    # This file contains the new grid cutoff
    
    f2 = open('cp2k_input_template_1.inp','w')
    
    for i in range(len(lines)):
        if i != cutoff_line_num[1]:
            f2.write(lines[i])
        else:
            f2.write('      CUTOFF %d\n'%new_cutoff)
    
    f2.close()
    os.system('rm cp2k_input_template.inp')
    os.system('mv cp2k_input_template_1.inp cp2k_input_template.inp')
    os.system('python run.py')
    print('Done with cutoff %d'%new_cutoff)
    os.system('mv res res-%dRy'%new_cutoff)


