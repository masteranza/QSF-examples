module load plgrid/tools/impi/2018.5 plgrid/tools/intel/2021.3.0 plgrid/tools/cmake/3.19.3 plgrid/libs/mkl/2021.3.0

rm -rf build; 
LDFLAGS="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl" CC=mpicc CXX=mpicxx cmake -DFFTW_INCLUDE_DIRS=$MKLROOT/include/fftw -DFFTW_DOUBLE_MPI_LIB=$MKLPATH/libfftw3x_cdft_gnu_lp64.a -DFFTW_DOUBLE_LIB=$MKLPATH/libfftw3xc_gnu.a . -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -B build; 

cmake --build build --target install; 