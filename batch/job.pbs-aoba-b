#!/bin/sh
#PBS -q lx
#PBS -T openmpi
#PBS -l elapstim_req=0:10:00
#PBS -b 2

export LD_LIBRARY_PATH=/uhome/a01152/.local/fftw3_amd/lib:/uhome/a01152/.local/netcdf_amd/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=/uhome/a01152/.local/fftw3_amd/lib:/uhome/a01152/.local/netcdf_amd/lib:$LIBRARY_PATH
cd $PBS_O_WORKDIR
mpirun $NQSV_MPIOPTS -x LD_LIBRARY_PATH -x LIBRARY_PATH -np 160 ./calliope < calliope.in > calliope.out.std
