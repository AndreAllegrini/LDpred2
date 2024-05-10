#!/bin/bash -l
#SBATCH --output=scoring.LDpred2.array_%A_%a.out
#SBATCH --error=scoring.LDpred2.array_%A_%a.err
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

echo I am job $SLURM_JOBID

sumstat=$(awk -F',' "NR==$SLURM_ARRAY_TASK_ID {print \$1}" ${WORKdir}/GWAS.csv) #take first column of input file
type=$(awk -F',' "NR==$SLURM_ARRAY_TASK_ID {print \$2}" ${WORKdir}/GWAS.csv) #take second column of input file

Rscript --vanilla ldpred2_auto_inf_qc.r -s $sumstat -t $type 
