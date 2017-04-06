#! /bin/sh
#$ -S /bin/sh
#$ -j y
#$ -N datanalysis
#$ -cwd
R CMD BATCH verif_replicats.R
