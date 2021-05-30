# Obtaining the `xyz` coordinates from `cif` files

In order to obtain the structure coordinates, we will use [VESTA](https://jp-minerals.org/vesta/en/) software. For more information on how to use VESTA software we highly 
recommend you to watch and use the useful videos on the [Nickel and Copper YouTube Channel](https://www.youtube.com/channel/UCmOHJtv6B2IFqzGpJakANeg).

When loading the `cif` file into VESTA, you can simply expand the cell and make supercells. In order to use the coordinates of these files in CP2K you will need to first export
the structures data to `.vasp` format with Cartesian coordinates and then use the `vasp_xyz.py` file to generate the `xyz` coordinates file. Here is how you need to use this
Python file:
```
from vasp_xyz import *
vasp_to_xyz('unit_cell.vasp')
vasp_to_xyz('supercell.vasp')
```
You will obtain the `unit_cell.xyz` and `supercell.xyz` files which can be used in CP2K. Also, the cell size will be printed out (the same as below) and you can copy and
paste them in `&CELL` section of the CP2K input.
```
[32, 96, 8, 16, 4] ['C', 'H', 'N', 'I', 'Pb']
A         8.4280004501         0.0000000000         0.0000000000
B         0.0000000000         8.9860000610         0.0000000000
C         0.0000000000         0.0000000000        26.2329998016
```

**Note:** You can also use the `cif` file directly in CP2K by setting the name of the `cif` file in `COORD_FILE_NAME` in `&TOPOLOGY` section and change the format to `cif` file by
changing the `COORD_FILE_FORMAT` to `CIF`. It is also recommended to add this part to your input in the `&TOPOLOGY` section in case you use the raw `cif` file:
```
&GENERATE
   REORDER .TRUE.
&END GENERATE
```

