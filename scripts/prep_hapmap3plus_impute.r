library(bigstatsr)
library(bigsnpr)

geno_path = '/cluster/projects/p471/people/andrea/LDpred2/geno_data/'
out_path = '/cluster/projects/p471/people/andrea/LDpred2/geno_data/'


cat('reading .bed...')
(rds <- snp_readBed2(bedfile = paste0(geno_path,'MoBaPsychGen_v1-ec-eur-batch-basic-qc-exc-hapmap3plus.bed'),
                     backingfile = paste0(out_path,"genoHapMap3plus_N200k"), 
                     ncores = 32))
obj.bigsnp <- snp_attach(rds)

str(obj.bigsnp, max.level = 2)

#impute missings
G <- obj.bigsnp$genotypes

counts <- big_counts(G)
counts[, 1:8]

cat('Imputing SNPs...')
G2 <- snp_fastImputeSimple(G, method = "mean2", ncores = NCORES)

big_counts(G2, ind.col = 1:8)
big_counts(G, ind.col = 1:8)

obj.bigsnp$genotypes <- G2
snp_save(obj.bigsnp)