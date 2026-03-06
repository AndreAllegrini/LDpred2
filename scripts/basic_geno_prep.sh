#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=02:00:00
#$ -l mem=8G
#$ -l tmpfs=8G
#$ -pe smp 10
#$ -N basic_geno_prep
#$ -m be

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

module -f unload compilers mpi gcc-libs
module load r/recommended
export R_LIBS=~/Rlibs

wdPath="/path/to/wd"
cd ${wdPath}

echo "Running job $JOB_ID in $PWD"

R --file=${wdPath}/prep_hapmap3plus_impute.r --vanilla