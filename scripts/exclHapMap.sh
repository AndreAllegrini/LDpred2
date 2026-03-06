#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=00:30:00
#$ -l mem=3G
#$ -pe smp 1
#$ -N exclHapMap
#$ -m be

module load plink/1.90b3.40
set -o errexit

IN="/data/genotypes/"          # dir with .bed/.bim/.fam
OUT="/data/output/"
GENO="geno_data"              # prefix of .bed/.bim/.fam 
LIST="/data/maps/map_hm3_plus.list" # SNP list


plink --bfile ${IN}${GENO} --extract ${LIST} --make-bed --out ${OUT}${GENO}-hapmap-only