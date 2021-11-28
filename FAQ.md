# Known problems

## MPI Window creation problems
```text
An error occurred in MPI_Win_create
reported by process [1909719041,0] on communicator
MPI_ERR_WIN: invalid window
MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort, and potentially your MPI job)
```

Try using a newer OPEN-MPI version: https://github.com/open-mpi/ompi
Install it from archive:

```bash
tar xf openmpi-<version>.tar.bz2
cd openmpi-<version>
./configure --prefix=<path> |& tee config.out
make -j 8 |& tee make.out
make install |& tee install.out
```
and create custom `configure-...` script pointing to the new `mpicc` and `mpicxx`