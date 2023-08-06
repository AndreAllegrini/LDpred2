#!/bin/bash -l
#SBATCH --output=/cluster/projects/p471/people/andrea/LDpred2/out/rds.N200k.out
#SBATCH --error=/cluster/projects/p471/people/andrea/LDpred2/out/rds.N200k.err
#SBATCH --account=p471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8GB

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

## Set up job enviroment:
source /cluster/bin/jobsetup
module purge
module add R/3.5.0

wdPath=/cluster/projects/p471/people/andrea/LDpred2/scripts

cd $wdPath

echo I am job $SLURM_JOBID

R --file=prep_hapmap3plus_impute.r --vanilla
