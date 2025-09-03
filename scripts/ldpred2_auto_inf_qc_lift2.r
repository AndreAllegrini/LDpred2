#! /usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(bigsnpr)
library(bigreadr)
library(data.table)

# this scripts uses info from the following resources:
#recommended GWAS QC
#see  https://github.com/privefl/paper-misspec/tree/main/code/prepare-sumstats
#see https://privefl.github.io/bigsnpr/articles/LDpred2.html
# ew LD ref - hap-map + 
#see also https://github.com/privefl/paper-infer/blob/main/code/example-with-provided-LD.R


option_list = list(
  make_option(c("-s", "--sumstats"), type="character", default=NULL, 
              help=
                "Name of GWAS summary statistics.\n
                Note sumstats should have *at least* the following header:\n
                case/control traits: CHR BP A2 A1 NCAS NCON BETA SE \n
                continuous traits: CHR BP A2 A1 N BETA SE\n", metavar="character"),
  
  make_option(c("-g", "--geno"), type="character", 
              default="/cluster/projects/p471/people/andrea/LDpred2/geno_data/genoHapMap3plus_N200k.rds", 
              help="(Path to) Genetic dataset in .rds format.\n 
              [default = %default]", metavar="character"),
  
  make_option(c("-t", "--type"), type="logical", 
              default = TRUE,
              help="Whether GWAS trait is case/control.\n 
              [default = %default]", metavar="logical"),
  
  make_option(c("-o", "--out"), type="character", 
              default="/cluster/projects/p471/people/andrea/LDpred2/out/", 
              help="(Path to) Output directory.\n 
              [default = %default]", metavar="character"),
  
  make_option(c("-d", "--sdir"), type="character", 
              default="/cluster/p/p471/cluster/people/andrea/LDpred2/sumstats/", 
              help="Sumstats directory.\n 
              [default = %default]", metavar="character"),
  
  make_option(c("-m", "--misc"), type="character", 
              default = "/cluster/projects/p471/people/andrea/LDpred2/misc/hapmap3plus/", 
              help = "(Path to) Directory including LD reference info file(.rds) with ld scores, and LD matrices by chromosome.\n
              [default= %default]\n
              Note the directory should have  the following structure:\n 
              hapmap3plus
              ├── LDref
              │   ├── LD_with_blocks_chr1.rds
              │   ├── LD_with_blocks_chr2.rds
              │   ├── LD_with_blocks_chr3.rds
              etc...
              └── map_hm3_plus.rds", metavar="character"),
  
  make_option(c("--maf"), type = "double", default = 0.01,
              help = "MAF threshold for QC. [default = %default]",
              metavar = "numeric"),
  
  make_option(c("--info"), type = "double", default = 0.6,
              help = "INFO threshold for QC. [default = %default]",
              metavar = "numeric"),
  
  make_option(c("-l", "--lift"), type = "character",
              default = FALSE,
              help = "Assuming test data genome build hg19. (checks for consistency)
              If not do lift. Values: 'hg18' OR 'hg38'. [default = %default]",
              metavar = "character"),
  
  make_option(c("-c", "--cores"), type="integer", 
              default = max(1, parallel::detectCores() - 1), 
              help="Number of cores. \n 
              [default = %default]", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print help
print_help(opt_parser)

#print options
print(opt)

out_path = opt$out
misc_path = opt$misc
NCORES = opt$cores
geno = opt$geno

#checks

# Sumstats not found
sumstats_path <- file.path(opt$sdir, opt$sumstats)
if (!file.exists(sumstats_path)) {
  stop("Sumstats file not found at: ", sumstats_path, call. = FALSE)
}

lift_val <- opt$lift
if (lift_val == "FALSE") lift_val <- FALSE
if (!isFALSE(lift_val)) {
  lift_val <- tolower(lift_val)# only if not FALSE
}


if (!isFALSE(lift_val) && !lift_val %in% c("hg18","hg38")) {
  stop("--lift must be one of: 'hg18', 'hg38', FALSE. You gave: ", opt$lift, call. = FALSE)
}

#type
opt$type <- isTRUE(opt$type)


start <- Sys.time() #start time

file_log <- paste0(opt$sumstats,".log") 

file.create(file_log) #create log file

cat(" -------------------------------------------------","\n",
    " Generate LDpred2-auto and infinitesimal scores. ","\n",
    " bugs and questions: a.allegrini@ucl.ac.uk ","\n",
    "-------------------------------------------------","\n",
    " ","\n", sep = " ", file=file_log, append=TRUE)

#print  options
sink(file = file_log, append = TRUE)
#print_help(opt_parser)

cat("#####","SELECTED OPTIONS:","#####\n"," ", "\n", file=file_log, append=TRUE)
opt
cat("##########\n"," ", "\n", file=file_log, append=TRUE)
sink()


cat(paste0("Analyses started at ", start),"\n"," ","\n", file=file_log,append=TRUE)

cat("Reading: ", opt$sumstats, " sumstats.", "\n"," ","\n", sep='',file=file_log,append=TRUE) 

# load sumstats and convert to LDpred header format + calculate effective sample size 

sumstats <- bigreadr::fread2(input = sumstats_path)

#sumstats <- bigreadr::fread2(input = paste0(sumstatDir,opt$s))

cat("Loaded sumstats have: ", dim(sumstats)[1], " rows and ", dim(sumstats)[2], " columns.", "\n"," ","\n", sep="",file=file_log,append=TRUE) 


cat("Sumstats header is: ", paste(names(sumstats), collapse = " ", sep = " "), "\n"," ","\n", sep="",file=file_log,append=TRUE) 

# basic QC 

#to add discard sample N whe .6 < .9 quant 
#aslo add per variant effective sample size if more than .5 discarded
#https://github.com/privefl/bigsnpr/issues/281


cat("Starting sumstats QC...","\n"," ","\n", file=file_log,append=TRUE)  


# filter out duplicates
dups <- vctrs::vec_duplicate_detect(sumstats[, c("CHR", "BP")])

sumstats <- sumstats[!dups, ]

#dup_count <- nrow(sumstats[duplicated(sumstats[, c("chr","pos")]), ])


cat("There were: ", sum(dups)," duplicated physical positions in GWAS data.\n" ," ","\n", sep='',file=file_log,append=TRUE)


# filter on MAF if present 
if ("MAF" %in% colnames(sumstats)) {
  
  # Ensure MAF ≤ 0.5
  sumstats$MAF <- ifelse(sumstats$MAF <= 0.5, sumstats$MAF, 1 - sumstats$MAF)
  
  n_maf_disc <- sum(sumstats$MAF < opt$maf, na.rm = TRUE)
  
  cat("N = ", n_maf_disc,
      " variants discarded because MAF < ", opt$maf, "\n",
      sep = '', file = file_log, append = TRUE)
  
  sumstats <- subset(sumstats, MAF >= opt$maf | is.na(MAF))
  
} else {
  cat("No MAF column provided.\n", file = file_log, append = TRUE)
}

# filter on INFO if present 
if ("INFO" %in% colnames(sumstats)) {
  
  n_info_disc <- sum(sumstats$INFO < opt$info, na.rm = TRUE)
  
  cat("N = ", n_info_disc,
      " variants discarded because INFO < ", opt$info, "\n",
      sep = '', file = file_log, append = TRUE)
  
  sumstats <- subset(sumstats, INFO >= opt$info | is.na(INFO))
  
} else {
  cat("No INFO column provided.\n"," ","\n", file = file_log, append = TRUE)
}


# effective sample size
if(opt$type == T){ 
  
  #check colummn names are correct
  if (!all(c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE") %in% names(sumstats))) {
    stop(paste0("Sumstats header is not correct. Minium header required for case/control traits is:\n        ",
                paste(c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE"), collapse=", ")))}
  
  #select columns and rename
  sumstats <- sumstats[,c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE")]
  names(sumstats) <- c("chr", "pos", "a0", "a1", "Ncas","Ncon","beta","beta_se")
  
  Ncas <- sumstats$Ncas 
  Ncon <- sumstats$Ncon
  
  sumstats$n_eff <- 4 / (1 / Ncas + 1 / Ncon)
  
}else{ #if not
  
  #check colummn names are correct
  if (!all(c("CHR","BP","A2","A1","N","BETA","SE") %in% names(sumstats))) {
    stop(paste0("Sumstats header is not correct. Minium header required for continuous traits is:\n        ",
                paste(c("CHR","BP","A2","A1","N","BETA","SE"), collapse=", ")))}
  
  #select columns and rename
  sumstats <- sumstats[,c("CHR","BP","A2","A1","N","BETA","SE")]
  names(sumstats) <- c("chr", "pos", "a0", "a1", "n_eff","beta","beta_se")
  
}


#sometimes this is an issue:
if(!is.integer(sumstats$chr)){
  sumstats$chr <- as.integer(sumstats$chr)
}


map_ldref <- readRDS(paste0(misc_path,"map_hm3_plus.rds")) #read reference map


# update build if needed (to hg18 or hg38) using positions for hapmap3+ set provided in ref file

# check overlap with all builds

overlap_hg18 <- sum(paste0(sumstats$chr, ":", sumstats$pos) %in% paste0(map_ldref$chr, ":", map_ldref$pos_hg18))
overlap_hg19 <- sum(paste0(sumstats$chr, ":", sumstats$pos) %in% paste0(map_ldref$chr, ":", map_ldref$pos))
overlap_hg38 <- sum(paste0(sumstats$chr, ":", sumstats$pos) %in% paste0(map_ldref$chr, ":", map_ldref$pos_hg38))

overlaps <- c(hg18 = overlap_hg18, hg19 = overlap_hg19, hg38 = overlap_hg38)
max_build <- names(which.max(overlaps))

cat("Sumstats overlap most with build:", max_build, "\n"," ","\n", 
    "Overlaps: hg18 =", overlap_hg18, ", hg19 =", overlap_hg19, ", hg38 =", overlap_hg38, "\n"," ","\n",
    file = file_log, append = TRUE)

# target set build
if (isFALSE(lift_val)) {
  target_build <- "hg19"  #no lifting
} else {
  target_build <- lift_val  # user specified build
  map_ldref$pos_hg19 <- map_ldref$pos # create hg19 pos column for consistency - easy fix otherwise merge breaks
}

# Harmonize sumstats to target_build
if (max_build != target_build) {
  # Merge sumstats with reference to lift positions
  
  cat("Sumstats build (",max_build ,") different from target build. Lifting to:", target_build, "\n"," ","\n", 
      file = file_log, append = TRUE)
  
  if (target_build != "hg19") {
    
    sumstats <- merge(
      sumstats,
      map_ldref[, c("chr", "pos", "pos_hg18","pos_hg19", "pos_hg38")],
      by.x = c("chr", "pos"),
      by.y = c("chr", paste0("pos_", max_build)),
      all = FALSE,
      sort = FALSE
    )
    
    sumstats$pos <- sumstats[[paste0("pos_", target_build)]]
    
  }
  
  # Update map_ldref positions to target_build
  map_ldref$pos <- if (target_build == "hg19") map_ldref$pos else map_ldref[[paste0("pos_", target_build)]]
}

info_snp <- tibble::as_tibble(snp_match(sumstats, map_ldref, return_flip_and_rev = TRUE))


cat("There were: ", sum( vctrs::vec_duplicate_detect(sumstats[, c("chr","pos")]))," duplicated physical positions in GWAS data.\n" ," ","\n", sep='',file=file_log,append=TRUE)

# match sumstats with reference data
cat(" N = ", dim(info_snp)[1], "SNPs have been matched with reference data (i.e. HapMap3 + )\n",
    "N =  ",sum(info_snp$`_FLIP_`), "SNPs were flipped\n",
    "N =  ",sum(info_snp$`_REV_`), "were reversed.\n",
    "\n",file=file_log,append=TRUE)

#drop NAs and make sure order is the same with SD file below
info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp))

#chi-squared GWAS
chi2 <- with(info_snp, (beta / beta_se)^2)
cat("Mean Chi^2 = ",mean(chi2,na.rm=T),".\n",sep='',file=file_log,append=TRUE)  


#lambda gc 
lgc <- median(chi2,na.rm=T)/qchisq(0.5,1)
cat("Lambda GC = ",lgc,".\n",sep='',file=file_log,append=TRUE)  

#SD REFERENCE

sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))

if(opt$type == T){  #if GWAS trait is binary 
  
  cat("treating GWAS trait as binary.\n",sep='',file=file_log,append=TRUE)  
  
  sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2)) #sumstats sd
  
}else{ #if it is continuous 
  
  cat("treating GWAS trait as continuous.\n",sep='',file=file_log,append=TRUE)  
  
  #consistent with e.g.: https://github.com/privefl/paper-misspec/blob/main/code/prepare-sumstats-bbj/height.R
  #although assumes sd(y) = 1 (beta std 1)
  #if not below sd(y) is reestimated 
  
  sd_ss = with(info_snp, 1 / sqrt(n_eff * beta_se^2 + beta^2))
  
  sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)
  
}

is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

cat(" ","\n", file=file_log,append=TRUE)

#make (temp) out dir
dir.create(paste0(out_path,opt$sumstats,'_',format(Sys.time(), '%d%B%Y')))

tmp <- paste0(out_path,opt$sumstats,'_',format(Sys.time(), '%d%B%Y'),'/')

#plot SD SS vs SD REF
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave(paste0(tmp,"sd_", opt$sumstats,".png"), width = 10, height = 7)


df_beta <- info_snp[!is_bad, ] #remove bad SNPs

cat(sum(is_bad, na.rm = T)," SNPs are bad.\n"," ","\n", sep='',file=file_log,append=TRUE)

cat("After QC there are: ", dim(df_beta)[1], " SNPs.\n"," ","\n",sep='',file=file_log,append=TRUE)

#bit to be fixed/double checked:

# stop if more than 50% variants have discordant SD
# see also https://github.com/privefl/bigsnpr/issues/281
if(sum(is_bad, na.rm = T) > (length(is_bad)*0.5)) {
  
  cat("WARNING: More than half the variants had a discordant SD. Imputing Neff.\n
      Double check your input sumstats: reference population, per-variant sample size, effective N.\n\n",file=file_log,append=TRUE) # see: https://github.com/privefl/bigsnpr/issues/281.
  if(opt$type == T){
    #see eq 4 and 5 here  https://www.sciencedirect.com/science/article/pii/S2666247722000525?via%3Dihub#sec3.2
    
    info_snp$n_eff_imp <- (4 / sd_ldref^2 - info_snp$beta^2) / info_snp$beta_se^2
    
    #correcting difference in per allele effect size
    info_snp$n_eff_imp <- with(info_snp, ifelse(n_eff_imp >= .67 * quantile(n_eff_imp,.90), n_eff_imp, NA))  #doi: 10.1093/bioinformatics/btw613
    
    info_snp$n_eff <- median(info_snp$n_eff_imp, na.rm=T)
    
    sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
    
  }else{
    
    cat("WARNING: Be extra sceptic of median imputation for continuous traits.\n
      Double check source code.\n\n",file=file_log,append=TRUE) 
    
    info_snp$n_eff_imp <- (1 / (sd_ldref^2 - info_snp$beta^2)) / info_snp$beta_se^2
    
    info_snp$n_eff_imp <- with(info_snp, ifelse(n_eff_imp >= .67 * quantile(n_eff_imp,.90), n_eff_imp, NA))  
    
    info_snp$n_eff <- median(info_snp$n_eff_imp, na.rm=T)
    
    sd_y = with(info_snp, sqrt(quantile(0.5 * (n_eff * beta_se^2 + beta^2), 0.01))) # https://github.com/privefl/bigsnpr/issues/349
    
    sd_ss = with(info_snp, sd_y / sqrt(n_eff * beta_se^2 + beta^2))
    
    sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)
  }
  
  cat("N = ",sum(is.na(info_snp$n_eff_imp)), " SNPs have been filtered to correct for low sample size.\n"," ","\n",sep='',file=file_log,append=TRUE)
  
  cat("new iputed median N = ", median(info_snp$n_eff_imp, na.rm=T) ,"\n", file=file_log,append=TRUE) 
  
  is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
  
  cat(sum(is_bad, na.rm = T)," SNPs are bad.\n"," ","\n", sep='',file=file_log,append=TRUE)
  
  
  #plot SD SS vs SD REF
  qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
    theme_bigstatsr() +
    coord_equal() +
    scale_color_viridis_d(direction = -1) +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from allele frequencies of the LD reference",
         y = "Standard deviations derived from the summary statistics",
         color = "Removed?")
  
  ggsave(paste0(tmp,"sd_", opt$sumstats,".medianN.png"), width = 10, height = 7)
  
  df_beta <- info_snp[!is_bad, ] #remove bad SNPs
  
  cat("After QC there are: ", dim(df_beta)[1], " SNPs.\n"," ","\n",sep='',file=file_log,append=TRUE)
  
}

#load test data
obj.bigsnp <- snp_attach(paste0(geno))

# shortcut for geno test set
G <- obj.bigsnp$genotypes

map_test <- dplyr::transmute(obj.bigsnp$map,
                             chr = as.integer(chromosome), pos = physical.pos,
                             a0 = allele2, a1 = allele1)

#remove duplicates (if any) in test data 
dups <- vctrs::vec_duplicate_detect(map_test[, c("chr","pos")])
if (any(dups)) {
  cat("There are: ", sum( vctrs::vec_duplicate_detect(map_test[, c("chr","pos")])), " duplicated physical positions in test data.\n" ," ","\n", sep='',file=file_log,append=TRUE)
  cat("removing duplicated...\n"," ","\n",file=file_log,append=TRUE) 
  map_test <-  map_test[!dups, ] }

# match with variants in test data and prepare files for matrix multiplication later
cat("Matching with test data...\n"," ","\n",file=file_log,append=TRUE)

map_pgs <- df_beta[1:4]; map_pgs$beta <- 1

map_pgs2 <- snp_match(map_pgs, map_test, return_flip_and_rev = T)

# match sumstats with reference data
cat(" N = ", dim(map_pgs2)[1], "SNPs have been matched with test data\n",
    "N =  ",sum(map_pgs2$`_FLIP_`), "SNPs were flipped\n",
    "N =  ",sum(map_pgs2$`_REV_`), "were reversed.\n",
    "\n",file=file_log,append=TRUE)

#in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_test[, c("chr", "pos")])
in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_pgs2[, c("chr", "pos")])

#cat("in_test vector = TRUE is: ", sum(in_test), " long.\n"," ","\n",sep='',file=file_log,append=TRUE)

df_beta <- df_beta[in_test, ]

cat("There are: ", dim(df_beta)[1], " SNPs in common with test data.\n"," ","\n",sep='',file=file_log,append=TRUE)

# LDSC
cat("Running LDSC...\n"," ","\n",file=file_log,append=TRUE)

(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))

h2_est <- ldsc[["h2"]]

cat(paste(names(ldsc),ldsc, sep ="=", collapse="; "),"\n"," ","\n", file=file_log, append=TRUE)

#ldsc ratio 
ldsc_int <- ldsc[["int"]]
mchi2 <- mean(chi2,na.rm=T)

ldscRATIO <- (ldsc_int - 1) /  (mchi2 - 1)

cat("Ratio = ",ldscRATIO,".\n"," ","\n",file=file_log,append=TRUE)  

cat("Running sparse matrix...\n"," ","\n",file=file_log,append=TRUE)

#sparse matrix
for (chr in 1:22) {
  
  cat(chr, ".. ", sep = "")
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
  
  corr_chr <- readRDS(paste0(misc_path,"LDref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, paste0(tmp,"corr_chr"), compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

file.size(corr$sbk) / 1024^3  # file size in GB



cat("Running LDpred2-inf...\n"," ","\n",file=file_log,append=TRUE)

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)

beta_inf2 <- cbind(df_beta[,1:4], beta_inf) #bind with chr:pos:a0:a1

#save weights
write.table(beta_inf2, paste0(tmp,opt$sumstats,"_beta_inf.txt"), col.names=T,row.names=F,quote=F)

pred_inf <- big_prodVec(G, beta_inf * map_pgs2$beta,
                        ind.col = map_pgs2[["_NUM_ID_"]],
                        ncores = NCORES)


final_pred_inf <- cbind(obj.bigsnp$fam,pred_inf)		#bind with test map				

#save final scores
write.table(final_pred_inf, paste0(tmp,opt$sumstats,"_pred_inf.txt"), col.names=T, row.names=F, quote = F)


cat("Running LDpred2-auto...\n"," ","\n", file=file_log,append=TRUE)

# LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               ncores = NCORES)

#save weights
saveRDS(multi_auto, file = paste0(tmp,opt$sumstats,"_multi_beta_auto.rds"))

#check that chains converged
auto <- multi_auto[[1]]  # first chain
plot_grid(
  qplot(y = auto$path_p_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)

ggsave(paste0(tmp,"auto_chains_", opt$sumstats,".png"), width = 12, height = 10)


# Filter outlier betas and average remaining ones - new recommended way
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- (range > (0.95 * quantile(range, 0.95))))
final_beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

#save final betas
final_beta_auto2 <- cbind(df_beta[,1:4], final_beta_auto) #bind with chr:pos:a0:a1
write.table(final_beta_auto2, paste0(tmp,opt$sumstats,"_final_beta_auto.txt"), col.names=T,row.names=F,quote=F)

#calculate score (matrix multiplication)
final_pred_auto <- big_prodVec(G, final_beta_auto * map_pgs2$beta, 
                               ind.col = map_pgs2[["_NUM_ID_"]],
                               ncores = NCORES)

final_pred_auto <- cbind(obj.bigsnp$fam,final_pred_auto)		#bind with test fam				

write.table(final_pred_auto, paste0(tmp,opt$sumstats,"_pred_auto.txt"), col.names=T, row.names=F, quote = F)

#remove sbk
file.remove(paste0(tmp, "corr_chr.sbk"))

end <- Sys.time()

cat(paste0("Analyses ended at ", end)," ","\n", file=file_log,append=TRUE)

cat(paste0("(Analyses took: ", round(as.numeric(difftime(end, start , units="mins")),digits=1) ," minutes)"),"\n"," ","\n", file=file_log,append=TRUE)

cat("###END###\n",file=file_log,append=TRUE)

#move logs and figures in out folder 
system(paste('mv ', file_log ,tmp, sep = " "))

