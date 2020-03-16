#!/bin/bash

#SBATCH -J unsteadyflodding	 # job name
#SBATCH -o sp_un.o           # name of the outfile file
#SBATCH -N 4                 # number of nodes requested
#SBATCH -n 192               # total number of mpi tasks requested
#SBATCH -p normal            # queue (partition) -- normal, development, etc.
#SBATCH -t 24:00:00          # run time (hh:mm:ss) - 1.5 hours

# Slurm email notifications are now working on Lonestar 5
#SBATCH --mail-user=jesantos@utexas.edu
#SBATCH --mail-type=end     # email me when the job finishes
#SBATCH -A pge-fracture     # allocation
# run the executable
ibrun ../../src/2-phase_LBM/ShanChen input_spherepack.xml
