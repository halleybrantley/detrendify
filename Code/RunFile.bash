#!/bin/tcsh

#BSUB -W 1200
#BSUB -n 2
source ~/R-3.4.2.csh
R CMD BATCH --vanilla filename
#BSUB -o out.%J
#BSUB -e err.%J
