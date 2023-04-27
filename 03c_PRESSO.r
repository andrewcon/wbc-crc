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
  "moloc"
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

# <<< MR-PRESSO >>>
# Load RData file with clumped summary stats for MVMR analysis
load("#.RData")

# Basically I only need the UVMR results which I will have to compile into a data frame and then input into the MR PRESSO
# First, UVMR. Make function just for this.
# Second, MVMR compared to wbc_all. Third, MVMR compared between immune cells. Make function using old MVMR.Rdata for this.

# Function and analysis for UVMR
uvmr_presso_fun <- function(res_df, trait_data) {
  snp_dat <- res_df[res_df$mr_keep == TRUE, ]
  cols_to_keep <- c("SNP", "beta.exposure", "se.exposure", "pval.exposure", "beta.outcome", "se.outcome", "pval.outcome")
  snp_dat <- snp_dat[, cols_to_keep]
  presso_res <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = snp_dat, NbDistribution = 3000,  SignifThreshold = 0.05)
  presso_df <- as.data.table(presso_res[[1]])
  presso_df[, Exposure := NULL]
  presso_global_pval <- presso_res[[2]]$`Global Test`$Pvalue
  outlier_nr <- length(presso_res[[2]]$`Distortion Test`$`Outliers Indices`)
  outlier_snp <- snp_dat[presso_res[[2]]$`Distortion Test`$`Outliers Indices`,]$SNP
  outlier_snp <- paste(outlier_snp, collapse = ";")
  presso_distortion_pval <- presso_res[[2]]$`Distortion Test`$Pvalue
  presso_df[ , `:=` (exposure = trait_data[1], outcome = trait_data[2], presso_global_pval = presso_global_pval, 
    outlier_nr = outlier_nr, outlier_snp = outlier_snp, presso_distortion_pval = presso_distortion_pval)]
  return(presso_df)
}

myCluster <- makeCluster(12, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

presso_list <- foreach(x = 1:48) %dopar% {
  uvmr_presso_fun(
    res_df = mr_list[[x]][[7]],
    trait_data = c(mr_list[[x]][[2]]$exposure[1], mr_list[[x]][[2]]$outcome[1])
  )
}

stopCluster(myCluster)

save(presso_list, file = "#.RData")
