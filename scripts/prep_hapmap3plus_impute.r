library(bigstatsr)
library(bigsnpr)

geno_path <- "/geno/path/"
out_path <- "/out/path/"
bed_file <- "geno.bed"
NCORES <- 32

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

cat('Reading .bed file...\n')
rds <- snp_readBed2(
  bedfile = paste0(geno_path, bed_file),
  backingfile = paste0(out_path, bed_file), 
  ncores = NCORES
)
obj.bigsnp <- snp_attach(rds)
str(obj.bigsnp, max.level = 2)

G <- obj.bigsnp$genotypes
counts <- big_counts(G)
counts[, 1:8]

cat('Imputing SNPs...\n')
G2 <- snp_fastImputeSimple(G, method = "mean2", ncores = NCORES)

big_counts(G2, ind.col = 1:8)

obj.bigsnp$genotypes <- G2
snp_save(obj.bigsnp)

cat("Imputed genotypes saved.\n")