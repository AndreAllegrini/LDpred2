# LDpred2

Workflow to generate `LDpred2` (-auto and infinitesimal) scores.

Please before using check out the tutorial by the developer: <https://privefl.github.io/bigsnpr/articles/LDpred2.html>

------------------------------------------------------------------------

## Files description:

-   `scripts/exclHapMap.sh` : subset genetic dataset to HapMap3+ SNPs (map_ref/map_hm3_plus.rds)

-   `scripts/prep_hapmap3plus_impute.r` : prepare dataset in a digestable format for LDpred2 (convert to .rds and impute SNPs)

-   `scripts/basic_geno_prep.sh` : job submission script for prep_hapmap3plus_impute.r

-   `scripts/subLDpred2.sh` : job submission script for ldpred2_auto_inf_qc.R

An array job script assuming you have a `.csv` file named 'GWAS.csv' in your working directory.

The first column should include a list of summary statistics names for which you wish to generate PRS, the second column should include values '`TRUE`' for case/control traits and '`FALSE`' for continuous traits. No header.

Example GWAS.csv for two summary statistics:

> neuroticism.gz,FALSE
>
> MDD.gz,TRUE

The job script will have to be modified manually for the number of PRS you wish to generate `#$ -t 1-[number of PRS]` e.g. `#$ -t 1-2`

Example usage:

> #!/bin/bash -l
>
> #SBATCH --time=6:00:00
>
> [...]
>
> #SBATCH --array=1-2
>
> sumstat=$(awk -F',' "NR==$SLURM_ARRAY_TASK_ID {print \\\$1}" \${WORKdir}/GWAS.csv) #take first column of input file
>
> type=$(awk -F',' "NR==$SLURM_ARRAY_TASK_ID {print \\\$2}" \${WORKdir}/GWAS.csv) #take second column of input file
>
> Rscript --vanilla ldpred2_auto_inf_qc.r -s \$sumstat -t \$type

-   `scripts/ldpred2_auto_inf_qc.r` : this is the actual LDpred2 script including the recommended QC + some optional light QC. Please look at flags and defaults below:

    -   Usage: `Rscript --vanilla ldpred2_auto_inf_qc.r -h`

```         
Options: -s CHARACTER, --sumstats=CHARACTER Name of GWAS summary statistics.

            Note sumstats should have *at least* the following header:

            case/control traits: CHR BP A2 A1 NCAS NCON BETA SE 

            continuous traits: CHR BP A2 A1 N BETA SE


-g CHARACTER, --geno=CHARACTER
    (Path to) Genetic dataset in .rds format.

          [default = /cluster/projects/p471/people/andrea/LDpred2/geno_data/genoHapMap3plus_N200k.rds]

-t LOGICAL, --type=LOGICAL
    Whether GWAS trait is case/control.

          [default = TRUE]

-o CHARACTER, --out=CHARACTER
    (Path to) Output directory.

          [default = /cluster/projects/p471/people/andrea/LDpred2/out/]

-d CHARACTER, --Sdir=CHARACTER
    Sumstats directory.

          [default = /cluster/p/p471/cluster/people/andrea/LDpred2/sumstats/]

-m CHARACTER, --misc=CHARACTER
    (Path to) Directory including LD reference info file(.rds) with ld scores, and LD matrices by chromosome.

          [default= /cluster/projects/p471/people/andrea/LDpred2/misc/hapmap3plus/]

          Note the directory should have  the following structure:

          hapmap3plus
          ├── LDref
          │   ├── LD_with_blocks_chr1.rds
          │   ├── LD_with_blocks_chr2.rds
          │   ├── LD_with_blocks_chr3.rds
          etc...
          └── map_hm3_plus.rds

-c NUMBER, --cores=NUMBER
    Number of cores. 

          [default = 32]

-h, --help
    Show this help message and exit
```

------------------------------------------------------------------------

## Output:

The script generates a folder within the specified output directory containing the following files:

*Weights*

-   infinitesimal model: \*\_beta_inf.txt

-   auto model: \*\_final_beta_auto.txt

*Polygenic scores*

-   infinitesimal model: \*\_pred_inf.txt

-   auto model: \*\_pred_auto.txt

*Log File*

-   A log file (\*.log) documenting QC steps and details about the summary statistics.
