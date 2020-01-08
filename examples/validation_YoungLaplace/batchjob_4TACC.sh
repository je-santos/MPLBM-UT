#!/bin/bash

#SBATCH -J conv_1e0         # job name
#SBATCH -o conv_1e0.o       # name of the outfile file
#SBATCH -N 4                # number of nodes requested
#SBATCH -n 192               # total number of mpi tasks requested
#SBATCH -p normal           # queue (partition) -- normal, development, etc.
#SBATCH -t 48:00:00         # run time (hh:mm:ss) - 1.5 hours

# Slurm email notifications are now working on Lonestar 5
#SBATCH --mail-user=jesantos@utexas.edu
#SBATCH --mail-type=end     # email me when the job finishes
#SBATCH -A pge-fracture     # allocation
# run the executable
ibrun ../../src/2-phase_LBM/ShanChen input_tubes.xml
