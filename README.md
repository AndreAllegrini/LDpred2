# LDpred2

Workflow to generate LDpred2 (-auto and infinitesimal) scores. 

Please before using check out the tutorial by the developer: https://privefl.github.io/bigsnpr/articles/LDpred2.html

Here you can find a session from the London Genetic Network Methods where Florian Privé discusses the method and options in more detail: https://osf.io/4cqwf/

For questions/bugs please raise an issue or contact: a.allegrini [at] ucl.ac.uk 


Files description: 

* scripts/exclHapMap.sh : subset genetic dataset to HapMap3+ SNPs (map_ref/map_hm3_plus.rds) 

* scripts/prep_hapmap3plus_impute.r : prepare dataset in a digestable format for LDpred2 (convert to .rds and impute SNPs)

* scripts/basic_geno_prep.sh : job submission script for  prep_hapmap3plus_impute.r

* scripts/subLDpred2.sh : job submission script for ldpred2_auto_inf_qc.R

This is an Array job script assuming you have a .csv file named 'GWAS.csv' in your working directory. 

The first column should include a list of summary statistics names for which you wish to generate PRS, the second column should include values 'TRUE' for case/control traits and 'FALSE' for continuous traits. 

The job script will have to be modified manually for the number of PRS you wish to generate  #$ -t 1-total number of PRS  e.g. #$ -t 1-10

* scripts/ldpred2_auto_inf_qc.R : this is the actual LDpred2 workflow including the recommended QC + some optional light QC. Please look at flags and defaults below: 

Options:


	-s CHARACTER, --sumstats=CHARACTER
		Name of GWAS summary statistics.

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
		
		
		
NOTE: I use UKB as default LD reference panel, avialable at: https://figshare.com/articles/dataset/LD_reference_for_HapMap3_/21305061