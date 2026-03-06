# LDpred2 Workflow

Workflow to generate `LDpred2` (-auto and infinitesimal) scores.

Before using check out the tutorial by the developer: <https://privefl.github.io/bigsnpr/articles/LDpred2.html>

------------------------------------------------------------------------

## 1. Files Description

| Script | Description |
|-----------------------------|-------------------------------------------|
| `scripts/exclHapMap.sh` | Subset genotype dataset to HapMap3+ SNPs (requires `map_ref/map_hm3_plus.rds`). |
| `scripts/prep_hapmap3plus_impute.r` | Prepare dataset for LDpred2: convert to `.rds` and impute missing SNPs. |
| `scripts/basic_geno_prep.sh` | Job submission script for `prep_hapmap3plus_impute.r`. |
| `scripts/ldpred2_auto_inf_qc_lift2.r` | Main R script to generate PRS. Use `Rscript --vanilla ldpred2_auto_inf_qc_lift2.r -h` for help. |
| `scripts/subLDpred2.sh` | Job submission script for `ldpred2_auto_inf_qc_lift2.r` (array job). |

------------------------------------------------------------------------

## 2. Input Files

1.  **list of GWAS summary statistics .CSV** (`GWAS.csv`)

-   First column: list of sumstat filenames\

-   Second column: trait type (`TRUE` = binary, `FALSE` = continuous)\

-   No header required

Example for two traits:

```         
neuroticism.gz,FALSE
MDD.gz,TRUE
```

3.  **Genotype file**: `.rds` file prepared with `prep_hapmap3plus_impute.r`\
4.  **LD reference folder**: HapMap3+ LD reference by chromosome

Expected structure:

```         
hapmap3plus/
├── LDref/
│   ├── LD_with_blocks_chr1.rds
│   ├── LD_with_blocks_chr2.rds
│   └── ...
└── map_hm3_plus.rds
```

------------------------------------------------------------------------

## 3. Main R Script Options (`ldpred2_auto_inf_qc_lift2.r`)

| Flag | Description | Default |
|-------------------|------------------------------|-----------------------|
| `-s, --sumstats` | Name of GWAS summary statistics file | \- |
| `-t, --type` | Trait type (`TRUE` = binary, `FALSE` = continuous) | TRUE |
| `-g, --geno` | Path to `.rds` genotype file | \- |
| `-o, --out` | Output directory | \- |
| `-d, --sdir` | Directory containing sumstats | \- |
| `-m, --misc` | LD reference folder (HapMap3+) | \- |
| `--maf` | MAF threshold for QC | 0.01 |
| `--info` | INFO threshold for QC | 0.6 |
| `-l, --lift` | Genome build (`hg19` default; `hg18` or `hg38` if lift-over needed) | FALSE |
| `-c, --cores` | Number of cores | 8 |
| `-h, --help` | Show help message | \- |

------------------------------------------------------------------------

## 4. Output Files

| Model         | Weights File           | Polygenic Score File |
|---------------|------------------------|----------------------|
| Infinitesimal | `_beta_inf.txt`        | `_pred_inf.txt`      |
| Auto          | `_final_beta_auto.txt` | `_pred_auto.txt`     |

-   **Log file** (`*.log`) documenting QC steps and summary statistics details.\
-   All output files are saved in `-o`.

------------------------------------------------------------------------

## 5. Job Array Example (SGE)

Update the array line to match the number of PRS to generate:

``` bash
#$ -t 1-2
```

Example script to extract sumstat and type for each array task:

``` bash
#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=06:00:00
#$ -l mem=8G
#$ -pe smp 8
#$ -t 1-2
#$ -N subLDpred2
#$ -m be

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export R_LIBS=~/Rlibs

module -f unload compilers mpi gcc-libs
module load r/recommended

set -o errexit

WORKDIR="/path/to/code"
cd "$WORKDIR"

echo "Running job $JOB_ID in $PWD for task $SGE_TASK_ID"

SUMSTAT=$(awk -F',' "NR==$SGE_TASK_ID {print \$1}" "${WORKDIR}/GWAS.csv")
TYPE=$(awk -F',' "NR==$SGE_TASK_ID {print \$2}" "${WORKDIR}/GWAS.csv")

GENOFILE="/path/to/genofile.rds"
OUTFILE="/path/to/output"
SUMSTATDIR="/path/to/sumstats"
LDREF="/path/to/ldref"

Rscript --vanilla ldpred2_auto_inf_qc_lift2.r \
  -s "$SUMSTAT" \
  -t "$TYPE" \
  -g "$GENOFILE" \
  -o "$OUTFILE" \
  -d "$SUMSTATDIR" \
  -m "$LDREF" \
  -c 8
```

------------------------------------------------------------------------

-   UKB HapMap3+ LD reference: <https://figshare.com/articles/dataset/LD_reference_for_HapMap3_/21305061>
