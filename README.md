# <span style="font-variant:small-caps;">Calliope</span>
<span style="font-variant:small-caps;">Calliope</span> is a three-dimensional pseudospectral code for solving magnetohydrodynamics (MHD). The code has two main features. The first is the adoption of the [P3DFFT](https://p3dfft.readthedocs.io) library to perform three-dimensional fast Fourier transforms (FFTs) with two-dimensional (pencil) decomposition of the computational domain. The second feature is that <span style="font-variant:small-caps;">Calliope</span> can solve turbulence driven by magnetorotational instability (MRI) in shearing coordinates using the remapping method (Rogallo 1981; Umurhan & Regev 2004). 

The physical models that <span style="font-variant:small-caps;">Calliope</span> can solve are (i) isothermal compressible MHD, (ii) incompressible MHD, and (iii) reduced MHD. The code structure is divided into two layers: common modules and model-dependent modules. One can easily modify and extend the above models by copying the model-dependent modules.


## Requirements

In addition to standard Fortran compiler and MPI libraries, <span style="font-variant:small-caps;">Calliope</span> requires the following:

- [FFTW3](http://www.fftw.org) (It may be OK with other FFT libraries, but we have not checked if it works.)
- [P3DFFT](https://p3dfft.readthedocs.io)
- [netCDF](https://docs.unidata.ucar.edu/netcdf-c/)
- [netCDF-Fortran](https://www.unidata.ucar.edu/software/netcdf/docs-fortran/)



## Getting the Code

```sh
git clone https://github.com/ykawazura/calliope.git
```



## Build Instructions

#### <u>Step1</u> : Install FFTW3
Follow the instruction linked above.
Let's say you installed it to `/usr/local/fftw3`.

#### <u>Step2</u> : Install P3DFFT
Download P3DFFT from the above link and install following the instruction. Specify FFTW3 for the FFT library. One can choose any build option (e.g., with or without OpenMP), except for stride-1 option which must be turned **<u>on</u>**.
Let's say you installed P3DFFT to `/usr/local/p3dfft`.

#### <u>Step3</u> : Install netCDF and netCDF-Fortran
Follow the instruction linked above.
Let's say you installed them to `/usr/local/netcdf`.

#### <u>Step4</u> : Create `arch/XXX.in`
You need to create a file having prefix `.in` in `arch/`. Let's say the name of this file is `XXX.in`. You then write make options in this file. For example, 

```makefile
FC = /usr/local/bin/mpif90

P3DFFT_HOME=/usr/local/p3dfft
P3DFFT_INC=-I$(P3DFFT_HOME)/include/
P3DFFT_LIB=-L$(P3DFFT_HOME)/lib/ -lp3dfft                  
                                                          
FFTW3_HOME=/usr/local/fftw3
FFTW3_LIB=-L$(FFTW3_HOME)/lib/ -lfftw3 -lfftw3f            
FFTW3_INC=-I$(FFTW3_HOME)/include/                         
                                                           
NETCDF_HOME=/usr/local/netcdf                              
NETCDF_INC=-I$(NETCDF_HOME)/include/                       
NETCDF_LIB=-L$(NETCDF_HOME)/lib -lnetcdf -lnetcdff         

F90FLAGS = -DGNU -O3 -Wunused                              
MPI_LIB = -lm -lmpi -lmpifort                              
LIB_OPENMP = -fopenmp # used only when USE_OPENMP=yes in ../Makefile.in
```

Templates for a number of machines can be found in `arch/`.

#### <u>Step5</u> : Edit `Makefile.in`
One needs to edit two lines in `Makefile.in`. 
First, you need to set `arch` to be the name of the file created above (`arch = XXX` for the above example).
Second, you need to choose `MODEL` from one of `RMHD`, `MHD_INCOMP`, or `MHD_COMP_ISOTH` depending on the model that you want to solve.

#### <u>Step6</u> : Do `make` to compile.
If the compilation is successful, you will find the execution file `calliope`.

#### <u>Step7</u> : Edit the input file `calliope.in`

#### <u>Step8</u> : Enjoy <span style="font-variant:small-caps;">Calliope</span>!



## Questions and/or Feedback

If you have questions and/or feedback about <span style="font-variant:small-caps;">Calliope</span>, please contact [Yohei Kawazura](https://sites.google.com/view/yoheikawazura/) ([kawazura@tohoku.ac.jp](mailto:kawazura@tohoku.ac.jp)).