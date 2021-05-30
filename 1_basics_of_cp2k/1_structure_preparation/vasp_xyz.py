import os
import sys
import numpy as np


def vasp_to_xyz(filename):
    
    f = open(filename,'r')
    lines = f.readlines()
    f.close()

    xyz_file_name = filename.replace('.vasp','')+'.xyz'
    f = open(xyz_file_name,'w')
    types = lines[5].split()
    n_types = [int(lines[6].split()[i]) for i in range(len(lines[6].split()))]
    print(n_types,types)
    coord = []
    for i in range(8,len(lines)):
        tmp_line = lines[i].split()
        x = float(tmp_line[0])
        y = float(tmp_line[1])
        z = float(tmp_line[2])
        coord.append([x,y,z])
    print('A',lines[2])
    print('B',lines[3])
    print('C',lines[4])
    
    f.write(str(np.sum(np.array(n_types)))+'\n\n')
    
    counter = 0
    for i in range(len(types)):
        for k in range(n_types[i]):
            f.write(types[i]+' '+str(coord[counter][0])+' '+str(coord[counter][1])+' '+str(coord[counter][2])+'\n')
            counter += 1
            
    f.close()

            


