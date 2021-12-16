# Compiling CP2K with Intel Parallel Studio


Here are the instruction for installing the CP2K software using Intel Parallel Studio compilers and are 
used for compiling the CP2K on UB-CCR cluster. Before going into main instructions, it should be noted that we use the instruction
from [**XConfigure**](https://xconfigure.readthedocs.io/en/latest/cp2k/) for compilation using Intel Parallel Studio 20.2. XConfigure is very useful
for configuration and generating the `make` files and we will apply our own changes to the files downloaded from XConfigure so that it 
can be run our target nodes.
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

while for _general-compute_ node we get the following:

```
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                32
On-line CPU(s) list:   0-31
Thread(s) per core:    1
Core(s) per socket:    16
Socket(s):             2
NUMA node(s):          2
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 85
Model name:            Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz
Stepping:              4
CPU MHz:               1045.257
CPU max MHz:           3700.0000
CPU min MHz:           1000.0000
BogoMIPS:              4200.00
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              1024K
L3 cache:              22528K
NUMA node0 CPU(s):     0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30
NUMA node1 CPU(s):     1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31
Flags:                 fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc art arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm 3dnowprefetch epb cat_l3 cdp_l3 invpcid_single intel_ppin intel_pt ssbd mba ibrs ibpb stibp tpr_shadow vnmi flexpriority ept vpid fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm cqm mpx rdt_a avx512f avx512dq rdseed adx smap clflushopt clwb avx512cd avx512bw avx512vl xsaveopt xsavec xgetbv1 cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local dtherm ida arat pln pts pku ospke md_clear spec_ctrl intel_stibp flush_l1d
```

The `Flags` shows us what type of architecture the CPU supports. For example, the _Valhalla_ does not support the `avx512` while the _general-compute_ node support it. 
In fact, the _general-compute_ supports more flags than _Valhalla_. So, if we compile the CP2K on _general-compute_ node using the flags for `avx512` and we submit a 
job on _Valhall_ we will get `SIGILL Illegal instructions` error and the job will be terminated since it does not support `avx512`.

Therefore, we have to find the proper flags when we want to compile a software. Since we want to run the compiled CP2K on _Valhalla_ as well we need to install
every dependent library including `FFTW3`, `Libint`, `Libxc`, and any other like `ELPA` using the flags that it supports it.

Now, we go to the main instructions. You can load Intel libraries including mpi and mkl using the following commands:
```
source /util/academic/intel/20.2/compilers_and_libraries_2020.2.254/linux/bin/compilervars.sh intel64
source /util/academic/intel/20.2/compilers_and_libraries_2020.2.254/linux/mpi/intel64/bin/mpivars.sh
source /util/academic/intel/20.2/compilers_and_libraries_2020.2.254/linux/mkl/bin/mklvars.sh intel64
```
Other alternatives are through `module load` which dependent on the cluster type, it might be different. On UB-CCR this can be done using
```
module load intel/20.2
```
As was mentioned, the rest of the procedure will use the XConfigure instructions. Now, you can make a directory named `cp2k-intel` and do the rest of the procedure.

## 1. Compile dependent libraries

### 1.1 Compile `ELPA`

`ELPA` can eficiently increase the speed of calculations. For compiling it you need to run the following commands:
```
cd cp2k-intel
wget --content-disposition --no-check-certificate https://www.cp2k.org/static/downloads/elpa-2020.11.001.tar.gz
tar xvf elpa-2020.11.001.tar.gz
cd elpa-2020.11.001
wget --content-disposition --no-check-certificate https://github.com/hfp/xconfigure/raw/master/configure-get.sh
chmod +x configure-get.sh
./configure-get.sh elpa
```
Before running the configuration files downloaded from XConfigure, we need to make some changes to them. This will be used for all of the libraries that
we want to compile. In the file `configure-elpa-skx-omp.sh`, change the `TARGET="-xCORE-AVX512 -qopt-zmm-usage=high"` to `TARGET="-xavx2 -qopt-zmm-usage=high"`. 
Remember this procedure for other configurations of other libraries as well. Now, do the following:
```
./configure-elpa-skx-omp.sh
# Compile with 12 processors
make -j 12
make install
make clean
```
It wil install the library in the previous folder in `elpa`.

### 1.2 Compile `Libint`

To compile `Libint` we do exactly as above and in the configureation file we change the `-xCORE-AVX512` to `-xavx2`. The same is done for compilation of `Libxc`.
```
cd cp2k-intel
curl -s https://api.github.com/repos/cp2k/libint-cp2k/releases/latest \
| grep "browser_download_url" | grep "lmax-6" \
| sed "s/..*: \"\(..*[^\"]\)\".*/url \1/" \
| curl -LOK-
tar xvf libint-v2.6.0-cp2k-lmax-6.tgz
cd libint-v2.6.0-cp2k-lmax-6
wget --content-disposition --no-check-certificate https://github.com/hfp/xconfigure/raw/master/configure-get.sh
chmod +x configure-get.sh
./configure-get.sh libint
# Change the configure-libint-skx.sh file by replacing -xCORE-AVX512 to -xavx2.
./configure-libint-skx.sh
make -j
make install
make distclean
```
The compilation of `Libint` might take up to an hour so please wait until it gets completely done.
### 1.3 Compile `Libxc`

The instructions for compiling the `Libxc` is the same as above but note that for CP2K-v8.2, we need `Libxc` with versions higher than 5.1.
```
wget --content-disposition --no-check-certificate https://www.tddft.org/programs/libxc/down.php?file=5.1.7/libxc-5.1.7.tar.gz
tar xvf libxc-5.1.7.tar.gz
cd libxc-5.1.7

wget --content-disposition --no-check-certificate https://github.com/hfp/xconfigure/raw/master/configure-get.sh
chmod +x configure-get.sh
./configure-get.sh libxc
# Change the configure-libxc-skx.sh file by replacing -xCORE-AVX512 to -xavx2.
./configure-libxc-skx.sh
make -j
make install
make distclean
```
### 1.4 Download `Libxsmm`

For this step, you just need to download the `Libxsmm`. The compilation of this package will be done when trying to compile CP2K.
```
wget --content-disposition --no-check-certificate https://github.com/hfp/libxsmm/archive/1.16.1.tar.gz
tar xvf 1.16.1.tar.gz
```

## 2. Compile CP2K

Finally, we want to compile CP2K using the libraries that we compiled with `-xavx2` flag. Download CP2K:
```
git clone --recursive -b support/v8.2 https://github.com/cp2k/cp2k.git cp2k-8.2
cd cp2k-8.2
```
Then, download the arch files from XConfigure:
```
cd arch
wget https://github.com/hfp/cp2k/raw/master/arch/Linux-x86-64-intelx.arch
wget https://github.com/hfp/cp2k/raw/master/arch/Linux-x86-64-intelx.popt
wget https://github.com/hfp/cp2k/raw/master/arch/Linux-x86-64-intelx.psmp
wget https://github.com/hfp/cp2k/raw/master/arch/Linux-x86-64-intelx.sopt
wget https://github.com/hfp/cp2k/raw/master/arch/Linux-x86-64-intelx.ssmp
cd ..
```
Now, there is one more thing left to do. The same as above change the file `arch/Linux-x86-64-intelx.arch` by finding and adding the `xavx2` to `TARGET`:
```
    else ifeq (2,$(AVX))
      TARGET := -march=core-avx2

# replace the TARGET

    else ifeq (2,$(AVX))
      TARGET := -xavx2

```
Save the arch file and start compiling:
```
rm -rf exe obj lib
make -j 12 ARCH=Linux-x86-64-intelx VERSION=psmp  AVX=2 \
LIBINTROOT=/path/to/cp2k-intel/libint/intel-skx \
LIBXCROOT=/path/to/cp2k-intel/libxc/intel-skx \
ELPAROOT=/path/to/cp2k-intel/elpa/intel-skx-omp
```
The code will be compiled with `-xavx2` flag and you after compilation you can test and run it on the target node which in here was _Valhalla_.

