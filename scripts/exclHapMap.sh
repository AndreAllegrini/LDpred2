#!/bin/bash
#SBATCH --output=/cluster/projects/p471/people/andrea/LDpred2/out/imp.hapmap3plus.out
#SBATCH --error=/cluster/projects/p471/people/andrea/LDpred2/out/imp.hapmap3plus.err
#SBATCH --account=p471
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB

module purge
module load plink/1.90b6.2
set -o errexit

IN="/cluster/projects/p471/data/genetic_data/MoBaPsychGen_v1/"
OUT="/cluster/projects/p471/people/andrea/LDpred2/geno_data/"

plink --bfile ${IN}MoBaPsychGen_v1-ec-eur-batch-basic-qc --extract /cluster/projects/p471/people/andrea/LDpred2/misc/hapmap3plus/hapmap3plus.snplist --make-bed --out ${OUT}MoBaPsychGen_v1-ec-eur-batch-basic-qc-exc-hapmap3plus
