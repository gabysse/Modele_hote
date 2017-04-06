#! /bin/sh
#$ -S /bin/sh
#$ -j y
#$ -N as_boots
#$ -cwd 
R CMD BATCH lancement_boots.R
