#!/bin/bash -l
## Nazwa zlecenia
#SBATCH -J qsf-nitrogen-3e
## Liczba alokowanych węzłów
#SBATCH -N 8
## Liczba zadań per węzeł (domyślnie jest to liczba alokowanych rdzeni na węźle)
#SBATCH --ntasks-per-node=16
## Ilość pamięci przypadającej na jeden rdzeń obliczeniowy (domyślnie 5GB na rdzeń)
#SBATCH --mem-per-cpu=4Gb
## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=72:00:00
## Nazwa grantu do rozliczenia zużycia zasobów
#SBATCH -A plgjonizacja5
## Specyfikacja partycji
#SBATCH -p plgrid
## Plik ze standardowym wyjściem
#SBATCH --output="output.out"
## Plik ze standardowym wyjściem błędów
#SBATCH --error="error.err"

srun /bin/hostname

module load plgrid/tools/impi/2018.5 plgrid/tools/intel/2021.3.0 plgrid/tools/cmake/3.19.3 plgrid/libs/mkl/2021.3.0

export LDFLAGS="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl"

cd $SLURM_SUBMIT_DIR

mpiexec ./qsf-nitrogen-3e-re -b
