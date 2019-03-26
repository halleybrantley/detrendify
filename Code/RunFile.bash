#!/bin/tcsh

#BSUB -W 1200
#BSUB -n 2
module load gurobi
module load R/gcc_4.8.5_R-3.5.1
R CMD BATCH --vanilla filename

