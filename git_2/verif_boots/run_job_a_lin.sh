#! /bin/sh
#$ -S /bin/sh
#$ -j y
#$ -N alphaalinoll
#$ -cwd
#$ -pe orte 16
R CMD BATCH res_as_boots.R
