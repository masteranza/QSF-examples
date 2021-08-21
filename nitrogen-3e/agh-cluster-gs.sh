#!/bin/bash -l
## Nazwa zlecenia
#SBATCH -J qsf-nitrogen-3e-gs
## Liczba alokowanych węzłów
#SBATCH -N 8
## Liczba zadań per węzeł (domyślnie jest to liczba alokowanych rdzeni na węźle)
#SBATCH --ntasks-per-node=4
## Ilość pamięci przypadającej na jeden rdzeń obliczeniowy (domyślnie 5GB na rdzeń)
#SBATCH --mem-per-cpu=18Gb
## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=48:00:00
## Nazwa grantu do rozliczenia zużycia zasobów
#SBATCH -A plgjonizacja5
## Specyfikacja partycji
#SBATCH -p plgrid
## Plik ze standardowym wyjściem
#SBATCH --output="output.out"
## Plik ze standardowym wyjściem błędów
#SBATCH --error="error.err"

srun /bin/hostname

module load plgrid/tools/cmake plgrid/libs/fftw/3.3.9 plgrid/libs/mkl/2021.3.0 plgrid/tools/intel/2021.3.0

cd $SLURM_SUBMIT_DIR

mpiexec ./qsf-nitrogen-3e-im -b
