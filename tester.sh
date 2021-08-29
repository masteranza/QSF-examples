#!/bin/bash
## declare an array variable
declare -a nodes=("512" "1024" "2048" "256")
declare -a borders=("4" "8" "16")
declare -a deltas=("0.05" "0.1" "0.2" "0.4")
declare -a fields=("0.12", "0.15", "0.18", "0.22", "0.27", "0.32", "0.40", "0.50")
set -e

cd nitrogen-3e;

for ns in "${nodes[@]}"
do
   echo "running ./qsf-nitrogen-3e-im -n ${ns}"
   mpirun --host localhost:4 ./qsf-nitrogen-3e-im -n ${ns}
   for bs in "${borders[@]}"
   do
      for ds in "${deltas[@]}"
      do
         mpirun --host localhost:4 ./qsf-nitrogen-3e-re -n ${ns} -b ${bs} -t ${ds} -f 0.4
         # for fs in "${fields[@]}"
         # do
         #    echo "running ./qsf-nitrogen-3e-re -n ${ns} -b ${bs} -d ${ds} -f ${fs}"
         #    mpirun --host localhost:4 ./qsf-nitrogen-3e-re -n ${ns} -b ${bs} -d ${ds} -f ${fs}
         # done
      done
   done
done
