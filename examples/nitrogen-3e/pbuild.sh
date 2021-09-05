#!/bin/bash
set -ex
module load plgrid/tools/cmake plgrid/libs/fftw/3.3.9 plgrid/libs/mkl/2021.3.0 plgrid/tools/intel/2021.3.0
prj=`basename "$PWD"`
cd ../..;
rm -rf build; 
CC=gcc CXX=g++ cmake . -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -B build; 
cmake --build build --target ${prj}; 
cmake --build build --target ${prj}-im; 
cmake --build build --target ${prj}-re; 
cd -;