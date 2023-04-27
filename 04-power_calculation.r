#!/usr/bin/env Rscript

# Load relevant libraries------
list_of_packages <- c(
  "foreach",
  "doParallel",
  "kableExtra",
  "data.table",
  "devtools",
  "tidyr",
  "remotes",
  "progressr",
  "TwoSampleMR",
  "ieugwasr",
  "MRPRESSO",
  "moloc",
  "MVMR"
)

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

if(length(new_packages) > 0){
  install.packages(new_packages, dep=TRUE)
}

# Load packages
for(package.i in list_of_packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}
rm(list = c("list_of_packages", "new_packages", "package.i"))

load("#.RData")

# <<< Define functions >>>
# Proportion of variance explained by SNP
pve_snp <- function(beta, se, af, n) {
  num <- 2*(beta^2)*af*(1-af)
  den <- 2*(beta^2)*af*(1-af) + (se^2)*2*n*af*(1-af)
  pve <- num/den
  return(pve)
}
# 2S UVMR Conditional F-statistic------
fstat_fun <- function(res_df, trait_data) {
  snp_dat <- res_df[res_df$mr_keep == TRUE, ]
  snp_dat <- as.data.table(snp_dat)
  snp_dat[ , `:=` (exposure = trait_data[1], outcome = trait_data[2])]
  cols_to_keep <- c("exposure", "outcome", "SNP", "beta.exposure", "se.exposure", "eaf.exposure", "samplesize.exposure", "pval.exposure")
  snp_dat <- snp_dat[, ..cols_to_keep]
  snp_dat$t_stat <- snp_dat$beta.exposure/snp_dat$se.exposure
  snp_dat$f_stat <- snp_dat$t_stat^2
  snp_dat[, weak_snp := ifelse(f_stat < 10, "Yes", "No")]
  snp_dat[, maf := ifelse(eaf.exposure > 0.50, 1 - eaf.exposure, eaf.exposure)]
  snp_dat[, pve := pve_snp(beta.exposure, se.exposure, maf, samplesize.exposure)]
  snp_dat[, pve := pve*100 ]
  return(snp_dat)
}
# 2S UVMR harmonised SNP stats------
snp_dat_fun <- function(res_df, trait_data) {
  snp_dat <- res_df[res_df$mr_keep == TRUE, ]
  snp_dat <- as.data.table(snp_dat)
  snp_dat[ , `:=` (exposure = trait_data[1], outcome = trait_data[2])]
  cols_to_keep <- c("exposure", "outcome", "SNP", 
                    "beta.exposure", "se.exposure", "eaf.exposure", "samplesize.exposure", "pval.exposure", 
                    "beta.outcome", "se.outcome", "eaf.outcome", "samplesize.outcome", "pval.outcome")
  snp_dat <- snp_dat[, ..cols_to_keep]
  return(snp_dat)
}


# 2S UVMR Conditional F-statistic------
myCluster <- makeCluster(4, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

uvmr_fstat_list <- foreach(x = 1:48) %dopar% {
  fstat_fun(
    res_df = mr_list[[x]][[7]],
    trait_data = c(mr_list[[x]][[2]]$exposure[1], mr_list[[x]][[2]]$outcome[1])
  )
}

stopCluster(myCluster)

uvmr_fstat_all <- rbindlist(uvmr_fstat_list)
fwrite(uvmr_fstat_all, "#.txt", quote = F, row.names = F, sep = "\t", na = NA)


# 2S UVMR SNP description------
myCluster <- makeCluster(4, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

uvmr_desc_list <- foreach(x = 1:48) %dopar% {
  snp_dat_fun(
    res_df = mr_list[[x]][[7]],
    trait_data = c(mr_list[[x]][[2]]$exposure[1], mr_list[[x]][[2]]$outcome[1])
  )
}

stopCluster(myCluster)

uvmr_desc_all <- rbindlist(uvmr_desc_list)
fwrite(uvmr_desc_all, "#.txt", quote = F, row.names = F, sep = "\t", na = NA)


# For MVMR, proportion of variance explained in a pair-wise manner
# So basically look at neutrophil SNPs and see how much PVE they explain in bas, eos, lym, mon
# Same for other WBC traits
summary(uvmr_fstat_all[exposure == "neu" & outcome == "overall", pve])
sum(uvmr_fstat_all[exposure == "eos" & outcome == "overall", pve])

# function for pulling out pve
uvmr_fstat_all <- fread("#.txt")

wbc_trait_vec <- c("bas", "eos", "lym", "mon", "neu", "wbc_all")
m1 <- matrix(nrow = 6, ncol = 6, dimnames = list(wbc_trait_vec, wbc_trait_vec))

for(x in 1:length(wbc_trait_vec)) {
  wbc_out_trait <- wbc_trait_vec[x]
  wbc_out_dat <- fread(paste0("#", wbc_out_trait, ".txt"))

  for(y in 1:length(wbc_trait_vec)) {
    wbc_exp_trait <- wbc_trait_vec[y]
    print(paste0("Studying instruments for ", wbc_exp_trait, " inside ", wbc_out_trait, "."))
    snps_wbc_trait <- uvmr_fstat_all[exposure == wbc_trait_vec[y] & outcome == "overall", ]

    sum_pve_trait1 <- sum(snps_wbc_trait$pve)

    if(wbc_out_trait == wbc_exp_trait) {
      sum_pve_trait2 <- sum_pve_trait1
      m1[y,x] <- sum_pve_trait2
    } else {
      wbc_mvmr_dat <- wbc_out_dat[wbc_out_dat$ID %in% snps_wbc_trait$SNP, ]
      wbc_mvmr_dat <- wbc_mvmr_dat[!duplicated(wbc_mvmr_dat$ID)]
      wbc_mvmr_dat[, maf := ifelse(eaf > 0.50, 1 - eaf, eaf)]
      wbc_mvmr_dat[, pve := pve_snp(beta, se, maf, n_samples)]
      wbc_mvmr_dat[, pve := pve*100]
      sum_pve_trait2 <- sum(wbc_mvmr_dat$pve)
      m1[y,x] <- sum_pve_trait2
    }
  }
}

snps_wbc_trait[!duplicated(snps_wbc_trait$SNP)]
wbc_mvmr_dat[!duplicated(wbc_mvmr_dat$ID)]

m1 <- as.data.table(m1, keep.rownames=TRUE)

fwrite(m1, "#.txt", quote = F, row.names = F, sep = "\t", na = NA)


