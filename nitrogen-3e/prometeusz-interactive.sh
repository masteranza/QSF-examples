#!/bin/bash
# Starts an interactive session on Prometeusz AGH cluster
srun -N 1 --cpus-per-task=8 --mem-per-cpu=1Gb -p plgrid-testing -t 0:30:0 --pty /bin/bash