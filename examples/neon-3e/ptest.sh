#!/bin/bash
set -e
./pbuild.sh
prj=`basename "$PWD"`
module load plgrid/tools/cmake plgrid/libs/fftw/3.3.9 plgrid/libs/mkl/2021.3.0 plgrid/tools/intel/2021.3.0
echo "running [mpiexec -np 8 ./${prj}-im $@ -g -r]"
mpiexec -np 8 ./qsf-${prj}-im $@ -r
echo "running [mpiexec -np 8 ./qsf-${prj}-re $@ -g -r]"
mpiexec -np 8 ./qsf-${prj}-re $@ -r
