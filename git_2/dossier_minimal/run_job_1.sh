#! /bin/sh
#$ -S /bin/sh
#$ -j y
#$ -N nomjob1
#$ -cwd
R CMD BATCH --no-save resolution_poly_1.R
