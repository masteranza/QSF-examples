#!/bin/bash
module load plgrid/tools/cmake plgrid/libs/fftw/3.3.9 plgrid/libs/mkl/2021.3.0 plgrid/tools/intel/2021.3.0
projectdir=`pwd`
projectname=`basename "$projectdir"`
cd ../..;
rm -rf build; 
CC=gcc CXX=g++ cmake . -DPROJECT=${projectname} -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -B build; 
cmake --build build --target install; 
cd -;