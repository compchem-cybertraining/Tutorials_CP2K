# Geometry optimization

We first need to obtain an optimized geometry of the structure. When we make a new structure, like heterostructures, we may also need to obtain an optimum cell as well. The cell 
optimization has only one more variable than geometry optimization which is pressure and the optimizer will let the cell changes its size after each cell optimization step. In 
CP2K it is better to perform both cell and geometry optimization together which is done by the keyword `TYPE DIRECT_CELL_OPT` in the  `&CELL_OPT` section of the input. We have the cell optimization inputs in `cell_optimization` folder. We will talk more on this in there.

For other structures, like supercells, usually obtained from `cif` files, the structure may be initially optimized when geometry optimization is performed but usually this is not the case 
because the energy calculations in each software package may be different. Therefore, you need to perform geometry optimization first. If the changes in the structure are massive
then you need to move to cell optimization as well. You can first use an initial cutoff value for this purpose, which is obtained through the convergence analysis for the initial structure,
and then after obtaining the optimized geometry one can perform the convergence analysis again to obtain another cutoff value and if the initial cutoff is not the same as the one obtained 
form convergence analysis you can redo the geometry optimization by the new cutoff value. This will not take longer than the previous geometry optimization and will finish sooner. 

Note that the CP2K optimizer works based on the forces for each of the geometry and cell optimization. In fact, the movement of the atoms is based on their computed forces. The 
lower the convergence criteria of the SCF cycle (`EPS_SCF`) the more accurate the computed forces are. This will lead to a better optimized geometry and will not confuse the 
optimizer. So, it is recommended that for the geometry or cell optimization to use a good cutoff value and a smaller target accuracy of SCF cycle defined by `EPS_SCF`, like `10E-8.0`. After you obtained the optimized 
geometry, for the MD you can use lower `EPS_SCF` values like `10E-6.0`. But the more accurate the forces in MD the more accurate the time-overlaps are and therefore the more accurate the nonadiabatic couplings will be. This is up to the user on which `EPS_SCF` value to choose and is totally dependent on the studied system and its computational cost.


The target forces and displacement for the geometry or cell optimizations are defined in the `&MOTION` and `&GEO_OPT` section with `MAX_FORCE` and `MAX_DR` keywords 
respectively. The force value unit is in `Hartree/Bohr` and the displacement unit is in `Bohr`. Here, the maximum force is set to 0.0003 Hartree/Bohr which is almost equivalent to 15 meV/A and the maximum displacement is set to 0.002 Bohr. When running the geometry optimization, it will print out the coordinates 
and forces in each step in `*-pos-1.xyz` and `*-frc-1.xyz` files. The last coordinates in the `*-pos-1.xyz` will be the optimized geometry. Also, `*.restart` files are 
produced that contain the information of the geometry optimization of the lsat step
and if the run is suddenly interrupted, you can change the extension to `.inp` by `mv` command and then run it again. This will continue the geometry optimization from the 
point it was interrupted. In this case, for faster SCF calculations of the first step after interruption, you can add the `WFN_RESTART_FILE_NAME` with the produced `wfn` file to the input and set the
`SCF_GUESS` to `RESTART`. The controls over the production of such files can be 
done in the [`&PRINT`](https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/PRINT.html) section of `&MOTION` section. The optimized geometry of the initial structure with a cutoff 
value of 500 Ry is obtained and uploaded above (`optimized_BA2_PbI4.xyz`). The `BFGS` algorithm is used for the geometry optimization in the inputs.

It is worth noting that it also possible to perform transition state geometry optimization with CP2K but we do not consider it here.
