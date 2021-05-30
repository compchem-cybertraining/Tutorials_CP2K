# Hybrid functionals with CP2K

CP2K can perform hybrid functional calculations with a relatively good speed. You need to install CP2K with `LBINT` and `LIBXC` packages to be able to perform hybrid functional 
calculations. To access a library of exchange-correlation functionals you can use [LIBXC](https://www.tddft.org/programs/libxc/functionals/previous/libxc-5.0.0/). For different versions of CP2K one
needs a specific version of `LIBXC` or `LIBINT`. So, it is recommended to use the `./install_tool_chain.sh` to compile CP2K or look for the correct version of these libraries and compile them manually and then use them in compilation of CP2K (for more information take a look over [this link](https://xconfigure.readthedocs.io/en/latest/cp2k/)). 
The key to do the hybrid functional calculation speed is the initial guess for the SCF cycle. 
In order to do so, we need to first obtain an initial converged PBE `wfn` file and then use the PBE `wfn` file for the initial SCF guess of the hybrid functional. Here, we have
provided three files: `pbe.inp`, `hse06.inp`, and `b3lyp.inp`. The `pbe.inp` file runs a pure functional calculations and after the complete convergence, produces a `wfn` file. We use this file as an initial 
guess for the `hse06.inp` and `b3lyp.inp` using the `SCF_GUESS RESTART` and `WFN_RESTART_FILE_NAME BA2_PbI4_PBE-RESTART.wfn`.
For HSE06 caclulations we have adopted an input from [here](https://www.cp2k.org/_media/events:2018_summer_school:cp2k-uk-stfc-june-2018-sanliang-ling.pdf). We can use the
PBE potentials for HSE06 but for B3LYP, we need to use the BLYP potentials. The BLYP potentials are available for different atoms in [`GTH_POTENTIALS`](https://github.com/mkrack/cp2k-data/blob/master/potentials/Goedecker/cp2k/GTH_POTENTIALS) file.

In hybrid functional calculations, we need to add this part for HSE06 for the `&XC` section (the `&VDW_POTENTIAL` part is added in the input):
```
&XC
  &XC_FUNCTIONAL
  &XWPBE
    SCALE_X -0.25
    SCALE_X0 1.0 
    OMEGA 0.11
  &END
  &PBE
    SCALE_X 0.0
    SCALE_C 1.0
  &END PBE
  &END XC_FUNCTIONAL
  &HF
    &SCREENING
      EPS_SCHWARZ 1.0E-6
      EPS_SCHWARZ_FORCES 1.0E-5
      SCREEN_ON_INITIAL_P FALSE
    &END SCREENING
    &INTERACTION_POTENTIAL
      CUTOFF_RADIUS 10
      POTENTIAL_TYPE SHORTRANGE
      OMEGA 0.11
      !T_C_G_DATA t_c_g.DAT
    &END INTERACTION_POTENTIAL
    ! Defines the maximum amount of memory [MiB] to be consumed by the full HFX module.
    !&MEMORY
    !  MAX_MEMORY  10000
    !  EPS_STORAGE_SCALING 0.1
    !&END MEMORY
    FRACTION 0.25
  &END HF
&END XC
```
For hybrid functionals, CP2K uses Auxiliary Density Matrix Method (ADMM). This needs to be added to the `&DFT` section:
```
&AUXILIARY_DENSITY_MATRIX_METHOD
  ! recommended, i.e. use a smaller basis for HFX
  ! each kind will need an AUX_FIT_BASIS_SET.
  METHOD BASIS_PROJECTION
  ! recommended, this method is stable and allows for MD. 
  ! can be expensive for large systems
  ADMM_PURIFICATION_METHOD MO_DIAG
&END
```
Note that the `MO_DIAG` in `ADMM_PURIFICATION_METHOD` is good and stable for MD but also note that this method works only for `OT` method. In order to do the TD-DFT calculations you need to set this to `NONE`.
In the `&DFT` section, this part should also be added. The files `BASIS_ADMM` or `BASIS_ADMM_MOLOPT` include the auxiliary basis set for each atom type.
There is no problem to add multiple `BASIS_SET_FILE_NAME` to the input and the CP2K will concatenate them but one should not add multiple keywords in CP2K like `POTENTIAL_FILE_NAME`.
```
BASIS_SET_FILE_NAME BASIS_MOLOPT
BASIS_SET_FILE_NAME BASIS_ADMM
BASIS_SET_FILE_NAME BASIS_ADMM_MOLOPT
POTENTIAL_FILE_NAME GTH_POTENTIALS
```
For each atom type, two basis set should be added for hybrid functional calculations. The first one is the usual basis set and the second one is the auxiliary basis set
and is defined with `AUX_FIT` after `BASIS_SET`. Here is an example for Pb atom:
```
&KIND Pb
  BASIS_SET DZVP-MOLOPT-SR-GTH 
  BASIS_SET AUX_FIT cFIT6
  POTENTIAL GTH-PBE-q4
&END KIND
```
For B3LYP, we need to use `LIBXC` [library functionals](https://www.tddft.org/programs/libxc/functionals/previous/libxc-5.0.0/) which is `XC_HYB_GGA_XC_B3LYP`. You can also use other B3LYP functionals such as `HYB_GGA_XC_CAM_B3LYP` which is CAM-B3LYP. Here is the `&XC` section for B3LYP (the `&VDW_POTENTIAL` part is added in the input):
```
&XC
  &XC_FUNCTIONAL
    &LIBXC
      FUNCTIONAL XC_HYB_GGA_XC_B3LYP
    &END LIBXC
  &END XC_FUNCTIONAL
  &HF
    &SCREENING
      EPS_SCHWARZ 1.0E-10
    &END
    !&MEMORY
    ! This is the maximum memory for each processor in MBi, I just comment it but
    ! you can obtain it through computing the memory you ask in the slurm, pbs, or bash file
    ! divided by the number of processors.
    !  MAX_MEMORY  10000
    !  EPS_STORAGE_SCALING 0.1
    !&END MEMORY
    FRACTION 0.20
  &END
&END XC

```
In the `&KIND` section we need to add `BLYP` potentials from `GTH_POTENTIALS`. Here is the example for Pb atom:
```
&KIND Pb
  BASIS_SET DZVP-MOLOPT-SR-GTH
  BASIS_SET AUX_FIT cFIT6
  POTENTIAL GTH-BLYP-q4
&END KIND
```
By setting the `&MO_CUBES` and printing out only the energies by `WRITE_CUBE .FALSE.`, you can see that the energy gaps have increased compared to PBE functional calculations.

In this table you can see the timings needed for running the hybrid functional calculations with and without using the PBE converged `wfn` file. Here we used 25 processors (CP2K v7.1). Although the timing is not so much different but this will show itself for larger systems or some other systems which PBE convergence is easier and the hybrid functional will be converged in a couple of steps. Overall, it can be used as an alternative for the convergence of hybrid functionals.

|  Functional | Elapsed time (s) | 
|---|---|
|PBE|  171.1    |
|HSE06 with PBE `wfn` file   |   387.5 (Initial step: 139.4, 16 steps with 15.5)   |
|B3LYP with PBE `wfn` file   |  359.3 (Initial step: 57, 16 steps with 18.9)       |
|HSE06 without PBE `wfn` file|   408.2 (Initial step: 138.1, 26 steps with 15.7)    |
|B3LYP without PBE `wfn` file|   532 (Initial step: 56.9, 28 steps with 19)   |

Note that for TD-DFT with hybrid functionals we need to use the CP2K v7 or higher. In lower versions, there is a known problem with convergence of the TD-DFT calculations
with ADMM which seems to be fixed in v7. The same as before it is needed to add this section for TD-DFT calculations:
```
&PROPERTIES
  &TDDFPT
     NSTATES     20            # number of excited states
     MAX_ITER    200           # maximum number of Davidson iterations
     CONVERGENCE [eV] 1.0e-5   # convergence on maximum energy change between iterations
     &MGRID
        NGRIDS 16
        CUTOFF 500 # separate cutoff for TDDFPT calc
     &END
     ! Only in case you have a tdwfn file from previous calculations
     !RESTART     .TRUE.
     !WFN_RESTART_FILE_NAME RESTART.tdwfn
  &END TDDFPT
&END PROPERTIES
```
We have also added this part to the input although it is not needed since from now on for TD-DFT calculations with hybrid functionals we use the CP2K v7.1.
```
&XC_GRID
  XC_DERIV SPLINE2_SMOOTH
&END XC_GRID
```
The timings for the TD-DFT calculations for the B3LYP are shown in the table below. Here we used 25 processors (CP2K v7.1).
| #Excited states | Time for each TD-DFT cycle (s) |Total elapsed time (s)| Maximum excitation energy (eV) |
|---|---|---|---|
|20   | 363  |  6183 | 0.82 |
| 40  |  720 |  8639 | 0.98  |

The results of the TD-DFT calculations with PBE functional was shown in the [tddft] section. The band gaps for different functionals are shown in the following table. Note that with hybrid functionals the states energy gaps are higher than the pure functional. The maximum excitation energy with 20 states is 0.82 eV while for PBE functional it was 0.48 eV. So, dependent on the experiments, we can see that we can adopt even lower number of states for this functional.

|Functional | Band gap (eV) |
|---|---|
|PBE| 2.416|
|HSE06|  3.193  |
|B3LYP|  3.522  |

Finally, we have provided the B3LYP excitation analysis below. For the B3LYP hybrid functional, the TD-DFT calculations shows more combination of states and the excitonics effects are higher than PBE functional for the 2D perovskite of (BA)2PbI4. Also, less dark states can be observed compared to PBE functional. Another important point is that using B3LYP, the initial excitation energy has reduced from 3.522 eV to 2.778 eV. We can conclude that although the band gap has not a good agreement with experiments, in TDDFT for the first excited state, it will get closer to the experiments.
```
         State    Excitation        Transition dipole (a.u.)        Oscillator
         number   energy (eV)       x           y           z     strength (a.u.)
         ------------------------------------------------------------------------
 TDDFPT|      1       2.77853  -6.3878E-01  2.2010E-02  2.6291E-03   2.78096E-02
 TDDFPT|      2       2.77859  -4.5049E+00  2.2011E-01 -2.0072E-04   1.38478E+00
 TDDFPT|      3       2.88722   1.4421E-01  3.7284E+00 -3.2952E-03   9.84789E-01
 TDDFPT|      4       2.88734   1.3254E-01  3.1516E+00  3.6755E-03   7.03842E-01
 TDDFPT|      5       3.26202   1.2761E-04 -1.0053E-04  5.8915E-01   2.77390E-02
 TDDFPT|      6       3.26743  -1.3091E-02  1.3220E-02  3.6638E-03   2.87844E-05
 TDDFPT|      7       3.36752   7.3394E-01 -2.8251E-03  1.4920E-03   4.44425E-02
 TDDFPT|      8       3.36764   3.8735E-01 -1.4807E-03 -2.2663E-03   1.23797E-02
 TDDFPT|      9       3.38688  -1.0092E-02  1.1057E-04 -1.1777E-03   8.56662E-06
 TDDFPT|     10       3.38773  -1.0241E-02  2.6882E-05  7.5653E-04   8.75167E-06
 TDDFPT|     11       3.45930   2.1823E-02  5.1342E-02 -2.4491E-03   2.64271E-04
 TDDFPT|     12       3.45935  -1.0575E-02 -2.5333E-02 -4.0798E-03   6.52786E-05
 TDDFPT|     13       3.49901  -5.4491E-04 -7.9557E-02  2.2108E-02   5.84503E-04
 TDDFPT|     14       3.50016  -9.8226E-04  3.7558E-02  4.7076E-02   3.11090E-04
 TDDFPT|     15       3.50425  -4.9576E-01  1.7375E-02 -2.3333E-03   2.11267E-02
 TDDFPT|     16       3.50435   1.4074E-01 -4.1877E-03 -9.7436E-03   1.71020E-03
 TDDFPT|     17       3.53754  -4.0929E-01  3.8167E-02 -1.3196E-02   1.46595E-02
 TDDFPT|     18       3.53758   3.4443E-01 -3.3623E-02 -1.5747E-02   1.04012E-02
 TDDFPT|     19       3.59276  -3.1553E-01 -6.0260E-02 -1.2514E-02   9.09686E-03
 TDDFPT|     20       3.59288  -6.1659E-02 -1.7895E-02  6.0494E-02   6.84968E-04
 
 -------------------------------------------------------------------------------
 -                            Excitation analysis                              -
 -------------------------------------------------------------------------------
        State             Occupied              Virtual             Excitation
        number             orbital              orbital             amplitude
 -------------------------------------------------------------------------------
             1   2.77853 eV
                               196                  197              -0.697928
                               195                  198              -0.682179
                               195                  197               0.165071
             2   2.77859 eV
                               196                  198               0.703815
                               195                  197               0.682360
                               195                  198               0.114876
                               196                  197               0.080460
             3   2.88722 eV
                               196                  199               0.593180
                               195                  199              -0.563239
                               195                  200               0.448126
                               196                  200              -0.327656
             4   2.88734 eV
                               196                  200               0.630738
                               195                  200               0.518995
                               195                  199               0.424936
                               196                  199               0.359793
             5   3.26202 eV
                               194                  198              -0.699627
                               193                  197              -0.693322
                               192                  199              -0.050140
             6   3.26743 eV
                               194                  197              -0.701916
                               193                  198              -0.681339
                               192                  200              -0.062935
                               194                  198              -0.062620
                               191                  199              -0.060351
                               193                  197               0.052876
             7   3.36752 eV
                               193                  199               0.674225
                               194                  200               0.547033
                               194                  199              -0.428411
                               193                  200              -0.213522
             8   3.36764 eV
                               193                  200               0.653452
                               194                  199               0.557588
                               194                  200               0.456895
                               193                  199               0.190553
             9   3.38688 eV
                               196                  197               0.535026
                               195                  197               0.511871
                               196                  198              -0.485427
                               195                  198              -0.444621
                               191                  197              -0.057595
                               191                  198               0.051086
            10   3.38773 eV
                               195                  198               0.551273
                               196                  198              -0.499880
                               195                  197               0.474951
                               196                  197              -0.449680
                               191                  198               0.058523
                               191                  197               0.053854
            11   3.45930 eV
                               189                  197               0.645364
                               190                  198              -0.551130
                               189                  198              -0.299869
                               190                  197               0.285248
                               188                  198              -0.154828
                               183                  197              -0.131507
                               191                  198              -0.090199
                               188                  197               0.088598
                               192                  197              -0.087026
                               184                  198              -0.062239
                               183                  198               0.056571
                               191                  197               0.054277
            12   3.45935 eV
                               189                  198              -0.596198
                               190                  197               0.591213
                               190                  198               0.370301
                               189                  197              -0.222220
                               188                  197               0.156127
                               183                  198               0.128848
                               192                  198               0.089990
                               188                  198               0.089901
                               191                  197               0.088250
                               184                  197               0.064108
                               183                  197               0.053007
            13   3.49901 eV
                               196                  199               0.690140
                               195                  200              -0.562784
                               195                  199               0.389710
                               196                  200              -0.192716
                               191                  199              -0.071694
            14   3.50016 eV
                               196                  200              -0.661778
                               195                  199               0.573682
                               195                  200               0.438083
                               196                  199              -0.151507
                               191                  200               0.073505
            15   3.50425 eV
                               192                  197               0.611516
                               191                  198               0.578358
                               185                  198               0.261862
                               187                  197               0.239243
                               191                  197              -0.224724
                               192                  198              -0.113872
                               180                  197               0.104106
                               185                  197              -0.097532
                               190                  199              -0.092130
                               189                  200               0.087989
                               190                  198              -0.085328
                               189                  197               0.085308
                               182                  197               0.074427
                               181                  198               0.071020
                               188                  198               0.055425
                               189                  199              -0.051069
            16   3.50435 eV
                               192                  198              -0.621110
                               191                  197              -0.573200
                               185                  197              -0.261141
                               187                  198              -0.241626
                               191                  198              -0.180602
                               192                  197              -0.155998
                               180                  198              -0.103366
                               190                  200               0.098196
                               190                  197               0.088335
                               189                  199              -0.087177
                               185                  198              -0.086518
                               189                  198              -0.080555
                               182                  198              -0.075093
                               181                  197              -0.072811
                               187                  197              -0.058079
                               188                  197              -0.057431
            17   3.53754 eV
                               187                  197              -0.331335
                               185                  198              -0.320217
                               185                  197               0.309113
                               187                  198               0.242637
                               189                  199               0.234600
                               192                  197               0.227315
                               181                  198              -0.226898
                               190                  199               0.224525
                               181                  197               0.216329
                               180                  197              -0.215923
                               191                  198               0.214823
                               191                  197              -0.209240
                               182                  197              -0.204859
                               189                  200              -0.186122
                               192                  198              -0.165615
                               180                  198               0.160969
                               182                  198               0.150028
                               190                  200              -0.120818
                               188                  198              -0.076619
                               188                  197               0.071959
                               178                  198              -0.071424
                               178                  197               0.068473
                               188                  199               0.064068
                               170                  197              -0.053221
            18   3.53758 eV
                               187                  198               0.341406
                               185                  197               0.311387
                               185                  198               0.292689
                               187                  197               0.256478
                               190                  200              -0.254291
                               192                  198              -0.239955
                               181                  197               0.232006
                               181                  198               0.220805
                               182                  198               0.208160
                               180                  198               0.205338
                               191                  197              -0.205012
                               189                  200               0.201701
                               191                  198              -0.190813
                               192                  197              -0.181170
                               189                  199               0.167277
                               182                  197               0.156418
                               180                  197               0.151527
                               190                  199              -0.146759
                               188                  197               0.085269
                               188                  198               0.082213
                               178                  197               0.072937
                               178                  198               0.069099
                               188                  200              -0.064521
                               170                  198               0.057331
            19   3.59276 eV
                               189                  200               0.580489
                               190                  199              -0.544695
                               185                  198              -0.270409
                               190                  200              -0.257681
                               187                  197              -0.252870
                               188                  199              -0.152812
                               183                  200              -0.115595
                               181                  198              -0.113863
                               192                  200              -0.112791
                               191                  199              -0.107560
                               180                  197              -0.107266
                               182                  197              -0.085507
                               190                  198              -0.076137
                               187                  198              -0.069312
                               182                  200              -0.068839
                               189                  197               0.067710
                               183                  197               0.063282
                               188                  200              -0.062864
                               181                  199              -0.059075
                               184                  199              -0.057249
                               187                  200              -0.053742
            20   3.59288 eV
                               189                  199              -0.610926
                               190                  200               0.542339
                               185                  197               0.273392
                               187                  198               0.248315
                               190                  199              -0.175353
                               188                  200               0.153678
                               183                  199               0.117094
                               192                  199               0.115018
                               181                  197               0.113449
                               191                  200               0.113036
                               180                  198               0.106934
                               182                  198               0.084799
                               190                  197               0.074227
                               189                  198              -0.069277
                               182                  199               0.068032
                               183                  198              -0.065816
                               188                  199              -0.059236
                               181                  200               0.059022
                               187                  197              -0.058905
                               184                  200               0.056979
                               187                  199               0.055724
 -------------------------------------------------------------------------------
```

