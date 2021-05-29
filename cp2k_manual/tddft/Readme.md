# TD-DFT calculations in CP2K

In order to perform TD-DFT calculations in CP2K, you need to the following the `&FORCE_EVAL` section of the input:
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
The `NGRIDS` and cutoff are chosen the same as in the SCF calculations although you can run another convergence analysis for that too. One can also use the `RESTART` for TD-DFT calculations. To this end, you will need to add the `.tdwfn` file from previous calculations in front of `WFN_RESTART_FILE_NAME`. In CP2K v6.1 one also needs to add
this part to the `&XC` section. For higher versions, this isn't required.
```
&XC_GRID
  XC_DERIV SPLINE2_SMOOTH
&END XC_GRID
```
The results of the TD-DFT calculations for a different number of excited states (`NSTATES`) is shown in the following table.

| NGRIDS  | #Excited states  | #Processors  | Each TD-DFT cycle (s) |Total TD-DFT time (s)   | Maximum excitation energy (eV)
|---|---|---|---|---|---|
|16   |20   |25   | 147  |  1625 | 0.48   |
|16   |40   |25   |  290 |  2360 | 0.59  |
|16   |60   |25   |  445 | 5815  |  0.67 |

To plot the excitation analysis results or the TD-DFT spectrum please refer to [this link](https://github.com/AkimovLab/Project_Libra_CP2K) or [this link](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP). Here are the results of the excitation analysis for 20 states. You can see that the excitonic effects are quite abundant in the 2D perovskite of (BA)2PbI4 (see the excitation analysis for the state 3 and 8 with high oscillator strength).
```
          State     Excitation      Transition dipole (a.u.)        Oscillator
          number    energy (eV)       x         y         z       strength (a.u.)
         ------------------------------------------------------------------------
 TDDFPT|       1        2.42059     -0.0001    0.0000   -0.0000        0.00000
 TDDFPT|       2        2.42085     -0.0383   -0.0000   -0.0000        0.00009
 TDDFPT|       3        2.48567      3.9752   -0.0000   -0.0000        0.96233
 TDDFPT|       4        2.48591     -0.0003   -0.0000    0.0000        0.00000
 TDDFPT|       5        2.54556      0.0000   -0.1891    0.0002        0.00223
 TDDFPT|       6        2.54567     -0.0000   -0.0008   -0.0389        0.00009
 TDDFPT|       7        2.60655      0.0000   -0.0009   -0.0249        0.00004
 TDDFPT|       8        2.60666     -0.0000   -4.1670    0.0000        1.10887
 TDDFPT|       9        2.65817     -0.0000    0.0383   -0.0001        0.00010
 TDDFPT|      10        2.65843      0.0000   -0.0001   -0.0541        0.00019
 TDDFPT|      11        2.74497      0.0000   -0.0000    0.5740        0.02216
 TDDFPT|      12        2.75169     -0.0001    0.4633    0.0000        0.01447
 TDDFPT|      13        2.78312     -0.0000    0.0000   -0.0000        0.00000
 TDDFPT|      14        2.78324     -0.0040   -0.0000   -0.0000        0.00000
 TDDFPT|      15        2.85127      0.0000    0.0000    0.0000        0.00000
 TDDFPT|      16        2.85161      0.0000   -0.0000   -0.0000        0.00000
 TDDFPT|      17        2.85422     -0.5811    0.0001    0.0000        0.02362
 TDDFPT|      18        2.85440      0.0005    0.0000   -0.0000        0.00000
 TDDFPT|      19        2.89646      0.0929   -0.0000    0.0000        0.00061
 TDDFPT|      20        2.90257     -0.0000   -0.0000   -0.0000        0.00000
 
 Excitation analysis

   State        Occupied         Virtual        Excitation
   number        orbital         orbital        amplitude
 ---------------------------------------------------------
        1
                  196             197            -0.750469
                  195             198             0.660885
        2
                  196             198            -0.712322
                  195             197             0.701848
        3
                  195             197             0.710587
                  196             198             0.700071
        4
                  195             198             0.748843
                  196             197             0.659032
        5
                  196             199            -0.741485
                  195             200             0.670943
        6
                  196             200             0.724880
                  195             199            -0.688856
        7
                  195             199             0.722921
                  196             200             0.686773
        8
                  195             200            -0.739530
                  196             199            -0.668913
        9
                  194             197             0.744550
                  193             198            -0.667487
       10
                  194             198             0.717351
                  193             197            -0.696684
       11
                  193             197             0.700575
                  194             198             0.680315
                  190             199            -0.109787
                  189             200            -0.105539
                  176             198            -0.069008
                  175             197            -0.056353
       12
                  193             198             0.723561
                  194             197             0.646142
                  190             200            -0.127803
                  189             199            -0.121773
                  176             197            -0.075483
                  175             198            -0.061465
       13
                  194             199             0.740768
                  193             200            -0.671699
       14
                  194             200             0.729074
                  193             199            -0.684391
       15
                  192             197             0.762871
                  191             198             0.646378
       16
                  192             198            -0.717402
                  191             197            -0.696625
       17
                  193             199            -0.672983
                  194             200            -0.630668
                  190             197             0.293750
                  189             198             0.207061
                  185             199             0.054274
       18
                  193             200             0.686613
                  194             199             0.619826
                  190             198            -0.271967
                  189             197            -0.223628
                  185             200            -0.053717
       19
                  190             197             0.771161
                  189             198            -0.627909
                  193             199             0.073464
                  194             200             0.069153
       20
                  191             197            -0.705714
                  192             198             0.684736
                  188             199             0.088463
                  187             200            -0.086682
                  180             200             0.066596
                  179             199             0.058514
```


