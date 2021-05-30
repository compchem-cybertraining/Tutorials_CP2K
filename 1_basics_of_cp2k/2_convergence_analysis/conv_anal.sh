#!/bin/bash

# Remember you need to have the CP2K loaded. This is done either by `module load` or by exporting the path.
# This is explained in Readme.md file

for cutoff in {100..1500..100}; do
    for ngrids in {4,6,8,10,12,14}; do
       # Substituting the CUTOFF value 
       sed -i "/CUTOFF /c\     CUTOFF $cutoff" energy.inp
       # Substituting the NGRIDS value 
       sed -i "/NGRIDS /c\     NGRIDS $ngrids" energy.inp
       echo "----------$cutoff-------"
       echo "Running for cutoff $cutoff Ry and NGRIDS $ngrids..."
       mpirun -np 25 cp2k.psmp -i energy.inp -o OUT-conv-anal-$cutoff-$ngrids.log
       echo "Done!"
       # This command will 
       awk 'c&&!--c;/SCF WAVEFU/{c=4}' OUT-conv-anal-$cutoff-$ngrids.log | awk '{print $1 "      " $4 "      " $6}'
       # You can also use this command as well
       #cat OUT-conv-anal-$cutoff.log | grep -A 4 'SCF WAVEFUNC'
    
    done
done


