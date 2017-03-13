#! /bin/sh
#$ -S /bin/sh
#$ -j y
#$ -N nomjob2
#$ -cwd
R CMD BATCH --no-save resolution_poly_2.R
