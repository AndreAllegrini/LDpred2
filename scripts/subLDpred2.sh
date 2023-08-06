#!/bin/bash -l
#SBATCH --output=/cluster/projects/p471/people/andrea/LDpred2/out/scoring.LDpred2.array_%A_%a.out
#SBATCH --error=/cluster/projects/p471/people/andrea/LDpred2/out/scoring.LDpred2.array_%A_%a.err
#SBATCH --account=p471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8GB
#SBATCH --array=1-2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

## Set up job enviroment:
source /cluster/bin/jobsetup
module purge
module add R/3.5.0

set -o errexit

WORKdir=/cluster/projects/p471/people/andrea/LDpred2/scripts

cd $WORKdir

echo I am job $SLURM_JOBID

sumstat=$(awk -F',' "NR==$SLURM_ARRAY_TASK_ID {print \$1}" ${WORKdir}/GWAS.csv)                                       #take first column of excel input file
type=$(awk -F',' "NR==$SLURM_ARRAY_TASK_ID {print \$2}" ${WORKdir}/GWAS.csv)                                         #take second column of excel input file

Rscript --vanilla ldpred2_auto_inf_qc.r -s $sumstat -t $type 
