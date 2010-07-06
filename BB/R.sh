#!/bin/sh
#$ -cwd
#$ -pe orte 4
#$ -N R-BB-test
#$ -V

. /usr/modules/init/bash
module load openmpi
module load R

R < /mnt/home/cblair/cluster_tests/R_tests/BB/supplement_A-par.R --no-save -i data/003SCF.dat

