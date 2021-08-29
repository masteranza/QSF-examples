#!/bin/bash
set -e
./prometeusz-build.sh
module load plgrid/tools/cmake plgrid/libs/fftw/3.3.9 plgrid/libs/mkl/2021.3.0 plgrid/tools/intel/2021.3.0
echo "running [mpiexec -np 8 ./qsf-nitrogen-3e-im $@ -c]"
mpiexec -np 8 ./qsf-nitrogen-3e-im $@ -c
echo "running [mpiexec -np 8 ./qsf-nitrogen-3e-re $@ -c]"
mpiexec -np 8 ./qsf-nitrogen-3e-re $@ -c
