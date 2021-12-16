# Compiling CP2K with Intel Parallel Studio


Here are the instruction for installing the CP2K software using Intel Parallel Studio compilers and are 
used for compiling the CP2K on UB-CCR cluster. 
First, we need to figure out what is the architecture of the CPU type that we want 
to run CP2K on. For example, the compilation might be successfully done on a node but 
it does not run on another node through submission of a CP2K job. Usually, the error is
`SIGILL Illegal instruction` when trying to run the CP2K compiled on a different node. This
can happen for many softwares, including VASP but it is not our focus at the moment.
In order to find the CPU architecture of a node, you can submit a job that runs the command `lscpu` or run a Python
script that prints the CPU type:

```
import platform
print(platform.processor())
```

Here we use `lscpu` to figure out what are the CPU architectures. As an example, the _Valhalla_ compute node prints the 
following data:

```
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                12
On-line CPU(s) list:   0-11
Thread(s) per core:    1
Core(s) per socket:    6
Socket(s):             2
NUMA node(s):          2
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 63
Model name:            Intel(R) Xeon(R) CPU E5-2620 v3 @ 2.40GHz
Stepping:              2
CPU MHz:               1378.857
CPU max MHz:           3200.0000
CPU min MHz:           1200.0000
BogoMIPS:              4794.43
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              256K
L3 cache:              15360K
NUMA node0 CPU(s):     0,2,4,6,8,10
NUMA node1 CPU(s):     1,3,5,7,9,11
Flags:                 fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm epb invpcid_single ssbd ibrs ibpb stibp tpr_shadow vnmi flexpriority ept vpid fsgsbase tsc_adjust bmi1 avx2 smep bmi2 erms invpcid cqm xsaveopt cqm_llc cqm_occup_llc dtherm ida arat pln pts md_clear spec_ctrl intel_stibp flush_l1d
```
