#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=06:00:00
#$ -l mem=8G
#$ -l tmpfs=8G
#$ -pe smp 8
#$ -t 1-16
#$ -N subLDpred2
#$ -m be

# Limit threads for reproducibility
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export R_LIBS=~/Rlibs

# Load R
module -f unload compilers mpi gcc-libs
module load r/recommended

# Stop on error
set -o errexit

# Working directory
WORKDIR="/path/to/wd"
cd "$WORKDIR"

echo "Running job $JOB_ID in $PWD for SGE task $SGE_TASK_ID"

# Get sumstat and type from CSV
SUMSTAT=$(awk -F',' "NR==$SGE_TASK_ID {print \$1}" "${WORKDIR}/sumstats_list.csv")
TYPE=$(awk -F',' "NR==$SGE_TASK_ID {print \$2}" "${WORKDIR}/sumstats_list.csv")

# Paths
LDREF="/path/to/ldref"          # LD reference
GENOFILE="/path/to/genofile.rds" #  bigsnpr object
OUTFILE="/path/to/output"        # out dir
SUMSTATDIR="/path/to/sumstats"   # sumstats dir

# Run LDpred2
Rscript --vanilla ldpred2_auto_inf_qc_lift2.r \
  -s "$SUMSTAT" \
  -t "$TYPE" \
  -m "$LDREF" \
  -l hg19 \
  -g "$GENOFILE" \
  -o "$OUTFILE" \
  -c 8 \
  -d "$SUMSTATDIR"