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
  "rlang",
  "ieugwasr",
  "MRPRESSO",
  "moloc",
  "ggplot2",
  "ggforestplot",
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

library(TwoSampleMR)
library("ieugwasr")

# Functions------
gecco.mvmr <- function(exp_path, out_path, out_name, reverse = FALSE) {
  # Format exposure data
  set.seed(1647540)
  if(!reverse) {
    print("Format exposure data")
    exposure_dat <- mv_extract_exposures_local(
      exp_path,
      sep = "\t",
      phenotype_col = "phenotype",
      snp_col = "ID",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "reference_allele",
      other_allele_col = "other_allele",
      pval_col = "p-value",
      ncase_col = "ncase",
      samplesize_col = "n_samples",
      min_pval = 1e-200,
      log_pval = FALSE,
      pval_threshold = 5e-08,
      clump_r2 = 0.001,
      clump_kb = 10000,
      harmonise_strictness = 2
    )

    # Format outcome data
    print("Format outcome data")
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = out_path,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "StdErr",
      effect_allele_col = "Allele1",
      other_allele_col = "Allele2",
      eaf_col = "Freq1",
      pval_col = "P.value",
      samplesize_col = "Neff"
    )
  } else {
    print("Format exposure data")
    exposure_dat <- mv_extract_exposures_local(
      exp_path,
      sep = "\t",
      phenotype_col = "phenotype",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "StdErr",
      eaf_col = "Freq1",
      effect_allele_col = "Allele1",
      other_allele_col = "Allele2",
      pval_col = "P.value",
      ncase_col = "ncase",
      samplesize_col = "Neff",
      min_pval = 1e-200,
      log_pval = FALSE,
      pval_threshold = 5e-08,
      clump_r2 = 0.001,
      clump_kb = 10000,
      harmonise_strictness = 2
    )

    # Format outcome data
    print("Format outcome data")
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = out_path,
      sep = "\t",
      snp_col = "ID",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "reference_allele",
      other_allele_col = "other_allele",
      eaf_col = "eaf",
      pval_col = "p-value",
      samplesize_col = "n_samples"
    )
  }

  # Harmonise the exposure and outcome data
  print("Harmonising data")
  mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

  # Perform the MVMR
  print("Performing MV MR analysis")
  res <- mv_multiple(mvdat)
  res <- res[[1]]
  res <- res[, c(2,4:8)]
  res$outcome <- out_name
  print("Exporting results...")
  list_mvmr <- list(res, mvdat)
  return(list_mvmr)
}

# Assign path to WBC and CRC data------
wbc_paths <- c(list.files("#", full.names=T), list.files("#", full.names=T)) 
names(wbc_paths) <- rep(c("bas", "eos", "lym", "mon", "neu", "wbc_all"), 2)
crc_paths <- c(list.files("#", full.names = T), list.files("#", full.names=T))
names(crc_paths) <- rep(c("colon", "distal", "female", "left", 
                      "male", "overall", "proximal", "rectal"), 2)

# Create list which will contain all MR data------
mr_list <- list(
  list(wbc_paths[1:5], crc_paths[1], names(crc_paths[1]), FALSE),
  list(wbc_paths[1:5], crc_paths[2], names(crc_paths[2]), FALSE),
  list(wbc_paths[1:5], crc_paths[3], names(crc_paths[3]), FALSE),
  list(wbc_paths[1:5], crc_paths[4], names(crc_paths[4]), FALSE),
  list(wbc_paths[1:5], crc_paths[5], names(crc_paths[5]), FALSE),
  list(wbc_paths[1:5], crc_paths[6], names(crc_paths[6]), FALSE),
  list(wbc_paths[1:5], crc_paths[7], names(crc_paths[7]), FALSE),
  list(wbc_paths[1:5], crc_paths[8], names(crc_paths[8]), FALSE),
  list(wbc_paths[c(1,6)], crc_paths[1], names(crc_paths[1]), FALSE),
  list(wbc_paths[c(2,6)], crc_paths[1], names(crc_paths[1]), FALSE),
  list(wbc_paths[c(3,6)], crc_paths[1], names(crc_paths[1]), FALSE),
  list(wbc_paths[c(4,6)], crc_paths[1], names(crc_paths[1]), FALSE),
  list(wbc_paths[c(5,6)], crc_paths[1], names(crc_paths[1]), FALSE),
  list(wbc_paths[c(1,6)], crc_paths[2], names(crc_paths[2]), FALSE),
  list(wbc_paths[c(2,6)], crc_paths[2], names(crc_paths[2]), FALSE),
  list(wbc_paths[c(3,6)], crc_paths[2], names(crc_paths[2]), FALSE),
  list(wbc_paths[c(4,6)], crc_paths[2], names(crc_paths[2]), FALSE),
  list(wbc_paths[c(5,6)], crc_paths[2], names(crc_paths[2]), FALSE),
  list(wbc_paths[c(1,6)], crc_paths[3], names(crc_paths[3]), FALSE),
  list(wbc_paths[c(2,6)], crc_paths[3], names(crc_paths[3]), FALSE),
  list(wbc_paths[c(3,6)], crc_paths[3], names(crc_paths[3]), FALSE),
  list(wbc_paths[c(4,6)], crc_paths[3], names(crc_paths[3]), FALSE),
  list(wbc_paths[c(5,6)], crc_paths[3], names(crc_paths[3]), FALSE),
  list(wbc_paths[c(1,6)], crc_paths[4], names(crc_paths[4]), FALSE),
  list(wbc_paths[c(2,6)], crc_paths[4], names(crc_paths[4]), FALSE),
  list(wbc_paths[c(3,6)], crc_paths[4], names(crc_paths[4]), FALSE),
  list(wbc_paths[c(4,6)], crc_paths[4], names(crc_paths[4]), FALSE),
  list(wbc_paths[c(5,6)], crc_paths[4], names(crc_paths[4]), FALSE),
  list(wbc_paths[c(1,6)], crc_paths[5], names(crc_paths[5]), FALSE),
  list(wbc_paths[c(2,6)], crc_paths[5], names(crc_paths[5]), FALSE),
  list(wbc_paths[c(3,6)], crc_paths[5], names(crc_paths[5]), FALSE),
  list(wbc_paths[c(4,6)], crc_paths[5], names(crc_paths[5]), FALSE),
  list(wbc_paths[c(5,6)], crc_paths[5], names(crc_paths[5]), FALSE),
  list(wbc_paths[c(1,6)], crc_paths[6], names(crc_paths[6]), FALSE),
  list(wbc_paths[c(2,6)], crc_paths[6], names(crc_paths[6]), FALSE),
  list(wbc_paths[c(3,6)], crc_paths[6], names(crc_paths[6]), FALSE),
  list(wbc_paths[c(4,6)], crc_paths[6], names(crc_paths[6]), FALSE),
  list(wbc_paths[c(5,6)], crc_paths[6], names(crc_paths[6]), FALSE),
  list(wbc_paths[c(1,6)], crc_paths[7], names(crc_paths[7]), FALSE),
  list(wbc_paths[c(2,6)], crc_paths[7], names(crc_paths[7]), FALSE),
  list(wbc_paths[c(3,6)], crc_paths[7], names(crc_paths[7]), FALSE),
  list(wbc_paths[c(4,6)], crc_paths[7], names(crc_paths[7]), FALSE),
  list(wbc_paths[c(5,6)], crc_paths[7], names(crc_paths[7]), FALSE),
  list(wbc_paths[c(1,6)], crc_paths[8], names(crc_paths[8]), FALSE),
  list(wbc_paths[c(2,6)], crc_paths[8], names(crc_paths[8]), FALSE),
  list(wbc_paths[c(3,6)], crc_paths[8], names(crc_paths[8]), FALSE),
  list(wbc_paths[c(4,6)], crc_paths[8], names(crc_paths[8]), FALSE),
  list(wbc_paths[c(5,6)], crc_paths[8], names(crc_paths[8]), FALSE)
)

# Run the MR analysis------
myCluster <- makeCluster(4, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

mvmr_between_five <- foreach(x = 1:8) %dopar% {
  gecco.mvmr(
    exp_path = mr_list[[x]][[1]],
    out_path = mr_list[[x]][[2]],
    out_name = mr_list[[x]][[3]],
    reverse = mr_list[[x]][[4]]
  )
}

save(mvmr_between_five, file = "#.RData")
gc()

mvmr_against_wbc_all <- foreach(x = 9:48) %dopar% {
  gecco.mvmr(
    exp_path = mr_list[[x]][[1]],
    out_path = mr_list[[x]][[2]],
    out_name = mr_list[[x]][[3]],
    reverse = mr_list[[x]][[4]]
  )
}

stopCluster(myCluster)

save(mvmr_against_wbc_all, file = "#.RData")

# 2S MVMR instrument strength------
# Load MVMR data
load("#.RData")

# Create pheontypic matrix for each combination of wbc type to overall WBC count
# Order is "bas" "neu" "eos" "lym" "mon"
wbc_matrix <- matrix(c(1, 0.1, 0.1, 0.2, 0.1,
                       0.1, 1, 0.1, 0.2, 0.2,
                       0.1, 0.1, 1, 0.2, 0.2,
                       0.2, 0.2, 0.2, 1, 0.3,
                       0.1, 0.2, 0.2, 0.3, 1), 
                       nrow = 5, ncol = 5, 
                       dimnames = list(c("bas", "neu", "eos", "lym", "mon"),
                                       c("bas", "neu", "eos", "lym", "mon")))


fstat_mvmr_between <- function(x) {
  r_input <- format_mvmr(
  BXGs = mvmr_between_five[[x]][[2]][[1]],
  BYG = mvmr_between_five[[x]][[2]][[4]],
  seBXGs = mvmr_between_five[[x]][[2]][[3]],
  seBYG = mvmr_between_five[[x]][[2]][[6]],
  RSID = row.names(mvmr_between_five[[x]][[2]][[1]]))

  mvmr_ple <- pleiotropy_mvmr(r_input)
  mvmr_strength <- strength_mvmr(r_input)

  exposures <- mvmr_between_five[[x]][[1]][[1]]

  mvmr_weak <- qhet_mvmr(r_input, wbc_matrix, CI=FALSE)
  bas_weak_beta <- mvmr_weak[[1]][[1]]
  weak <- "Yes"

  mvmr_snp_dat <- data.table(
    exposure = exposures, 
    exposure_Fstat = c(c(mvmr_strength)[[1]], c(mvmr_strength)[[2]], c(mvmr_strength)[[3]], c(mvmr_strength)[[4]], c(mvmr_strength)[[5]]),
    outcome = mvmr_between_five[[x]][[1]][[2]][[1]],
    weak = c("Yes", "No", "No", "No", "No"),
    q_pval = mvmr_ple[[2]],
    weak_beta = c(bas_weak_beta, NA, NA, NA, NA)
  )
  mvmr_snp_dat <- mvmr_snp_dat[order(exposure)]
  return(mvmr_snp_dat)
}

myCluster <- makeCluster(4, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

mvmr_fstat_list <- foreach(x = 1:8) %dopar% {
  fstat_mvmr_between(
    x = x
  )
}

stopCluster(myCluster)

mvmr_fstat_all <- rbindlist(mvmr_fstat_list)

# Format table and output it as document
mvmr_fstat_all_fin <- mvmr_fstat_all
mvmr_fstat_all_fin[, or_weak := ifelse(is.na(weak_beta), NA, exp(weak_beta)) ]
mvmr_fstat_all_fin <- mvmr_fstat_all_fin[outcome != "left", ]
names(mvmr_fstat_all_fin) <- c("Exposure", "Exposure_Fstat", "Outcome", "Weak", "Het P-value", "Beta weak", "OR weak")
# mvmr_fstat_all_fin[, Weak := rep(c("Yes", "No", "No", "No", "No"), 7)]

mvmr_fstat_all_fin <- mvmr_fstat_all_fin[, Exposure := fcase(
  Exposure == "bas", "Basophil",
  Exposure == "eos", "Eosinophil",
  Exposure == "lym", "Lymphocyte",
  Exposure == "mon", "Monocyte",
  Exposure == "neu", "Neutrophil"
)][, Outcome := fcase(
  Outcome == "overall", "Overall",
  Outcome == "colon", "Colon",
  Outcome == "distal", "Distal",
  Outcome == "female", "Female",
  Outcome == "male", "Male",
  Outcome == "proximal", "Proximal",
  Outcome == "rectal", "Rectal"
)]
mvmr_fstat_all_fin$`Het P-value` <- format.pval(mvmr_fstat_all_fin$`Het P-value`)
mvmr_fstat_all_fin <- mvmr_fstat_all_fin %>% mutate_at(vars(`Beta weak`, `OR weak`, `Exposure_Fstat`), list(~ signif(., 2)))

# Combine with MVMR results to compare estimates
mvmr_between_five <- lapply(mvmr_between_five, `[[`, 1)
mvmr_between_five <- rbindlist(mvmr_between_five)
#mvmr_between_five$nsnp <- NULL
mvmr_between_five$method <- "IVW MVMR"
mvmr_between_five <- mvmr_between_five[, c("exposure", "outcome", "method", "b", "se", "pval")]
mvmr_between_five <- mvmr_between_five[outcome != "left"]
mvmr_between_five <- mvmr_between_five[, exposure := fcase(
  exposure == "bas", "Basophil",
  exposure == "eos", "Eosinophil",
  exposure == "lym", "Lymphocyte",
  exposure == "mon", "Monocyte",
  exposure == "neu", "Neutrophil"
)][, outcome := fcase(
  outcome == "overall", "Overall",
  outcome == "colon", "Colon",
  outcome == "distal", "Distal",
  outcome == "female", "Female",
  outcome == "male", "Male",
  outcome == "proximal", "Proximal",
  outcome == "rectal", "Rectal"
)]
mvmr_between_five$OR_IVW <- exp(mvmr_between_five$b)
mvmr_between_five <- mvmr_between_five[, c("exposure", "outcome", "OR_IVW")]
mvmr_fstat_all_fin <- merge(mvmr_fstat_all_fin, mvmr_between_five, by.x=c("Exposure", "Outcome"), by.y=c("exposure", "outcome"))
mvmr_fstat_all_fin <- mvmr_fstat_all_fin %>% mutate_at(vars(`Beta weak`, `OR weak`, `OR_IVW`, `Exposure_Fstat`), list(~ signif(., 2)))


# Save as txt table
fwrite(mvmr_fstat_all_fin, "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Save as docx table
set_flextable_defaults(
  font.size = 11, font.family = "Arial",
  padding = 4)
mvmr_fstat_ft <- flextable(mvmr_fstat_all_fin) %>% 
    theme_zebra() %>% 
    autofit() %>% set_table_properties(layout = "autofit")
save_as_docx(mvmr_fstat_ft, path = "#.docx")

# Figures for the MV MR analysis------
# Use Nightingale forestplot package
# First do it when comparing
# Load vs_wbc_all data
source("00_ggforestplot.r")
analyses_list_all <- fread("#.txt")
analyses_list_all <- analyses_list_all[model == "Model 1" & exposure != "wbc",]
analyses_list_all <- analyses_list_all[type %in% c("univariable", "multivariable_between"), ]
analyses_list_all$method <- ifelse(analyses_list_all$type == "univariable", "Observational", "Observational MV") 
analyses_list_all <- analyses_list_all[, c("exposure", "outcome", "method", "b", "se", "pval")]

load("#.RData")
uvmr_dat <- lapply(mr_list, `[[`, 1)
uvmr_dat <- rbindlist(uvmr_dat)
uvmr_dat <- uvmr_dat[, c("exposure", "outcome", "method", "b", "se", "pval")]
uvmr_dat <- uvmr_dat[! exposure %in% c("overall", "colon", "distal", "female", "left", "male", "proximal", "rectal"), ]
uvmr_dat <- uvmr_dat[method == "Inverse variance weighted", ]
uvmr_dat <- uvmr_dat[exposure != "wbc_all", ]
uvmr_dat$method <- "IVW UVMR"

load("#.RData")
mvmr_between_five <- lapply(mvmr_between_five, `[[`, 1)

curate_mvmr <- function(df) {
  df <- df[order(df$exposure),]
  df$effect_type <- c("direct")
  df <- df[, c(1,2,7,3,4,5,6)]
  return(df)
}
mvmr_between_five <- lapply(mvmr_between_five, curate_mvmr)

mvmr_between_five <- rbindlist(mvmr_between_five)
mvmr_between_five$method <- "MVMR between WBC subtypes"
mvmr_between_five <- mvmr_between_five[, c("exposure", "outcome", "effect_type", "method", "b", "se", "pval")]
mvmr_between_five$method <- "IVW MVMR (direct)"
mvmr_between_five$effect_type <- NULL
mvmr_between_five <- rbind(uvmr_dat, mvmr_between_five)

# Now combine all the data
mvmr_all_res <- rbind(mvmr_between_five, analyses_list_all)
mvmr_all_res <- mvmr_all_res[outcome != "left"]
mvmr_all_res <- mvmr_all_res[, exposure := fcase(
  exposure == "bas", "Basophil",
  exposure == "eos", "Eosinophil",
  exposure == "lym", "Lymphocyte",
  exposure == "mon", "Monocyte",
  exposure == "neu", "Neutrophil"
)][, outcome := fcase(
  outcome == "overall", "Overall",
  outcome == "colon", "Colon",
  outcome == "distal", "Distal",
  outcome == "female", "Female",
  outcome == "male", "Male",
  outcome == "proximal", "Proximal",
  outcome == "rectal", "Rectal"
)]
mvmr_all_res
mvmr_all_res <- mvmr_all_res[method != "Observational MV",]
mvmr_all_res$outcome <- factor(mvmr_all_res$outcome, 
                             levels=c('Overall','Male','Female','Colon','Proximal','Distal','Rectal'))
mvmr_all_res$method <- factor(mvmr_all_res$method, 
                               levels=rev(c('Observational', 'IVW UVMR', 'IVW MVMR (direct)')))

mvmr_ivw_res <- mvmr_all_res[method=="IVW UVMR" | method=="IVW MVMR (direct)",]
mvmr_ivw_res <- generate_odds_ratios(mvmr_ivw_res)
fwrite(mvmr_ivw_res, "#.txt", quote=F, row.names=F, na=NA, sep="\t")


source("00_ggforestplot.r")

#Custom palette for keeping colours consistent

custom_palette <- c(`IVW MVMR (direct)` = "#00C0C0", `IVW UVMR` = "#FF4A4F", `Observational` = "#323D5A")

# overall CRC only, for paper
all_tf_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 1.8), xtickbreaks = seq(0.0, 1.8, 0.1)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "Multivariable method", shape = "Multivariable method") + facet_grid(. ~ outcome) + 
  theme(legend.position="top") + guides(color=guide_legend(nrow=2, byrow=FALSE, reverse = TRUE), shape=guide_legend(nrow=2, byrow=FALSE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette) 

ggsave(filename = "#.png", 
       plot = all_tf_mvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

# Overall, male, female
all_tf_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall", "Male", "Female") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 2.5), xtickbreaks = seq(0.0, 2.0, 0.4)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "Multivariable method", shape = "Multivariable method") + facet_grid(. ~ outcome) + 
  theme(legend.position="top") + guides(color=guide_legend(nrow=2, byrow=FALSE, reverse = TRUE), shape=guide_legend(nrow=2, byrow=FALSE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette) 

ggsave(filename = "#.png", 
       plot = all_tf_mvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

# Second do the overall, colon, proximal, distal, rectal
all_sp_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall", "Colon", "Proximal", "Distal", "Rectal") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 2.4), xtickbreaks = seq(0.0, 2.0, 0.4)
) + theme(axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "Multivariable method", shape = "Multivariable method") + facet_grid(. ~ outcome) + 
  theme(legend.position="top") + guides(color=guide_legend(nrow=2, byrow=TRUE, reverse = TRUE), shape=guide_legend(nrow=2, byrow=TRUE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette) 

ggsave(filename = "#.png", 
       plot = all_sp_mvmr,
       scale = 1, width = 270,
       device = "png",
       height = 210, units = c("mm"),
       dpi = 250, bg = "white")



####### Observational last
# All figures, but with observational as last analysis
#Custom palette for keeping colours consistent
# Now combine all the data
mvmr_all_res <- rbind(mvmr_between_five, analyses_list_all)
mvmr_all_res <- mvmr_all_res[outcome != "left"]
mvmr_all_res <- mvmr_all_res[, exposure := fcase(
  exposure == "bas", "Basophil",
  exposure == "eos", "Eosinophil",
  exposure == "lym", "Lymphocyte",
  exposure == "mon", "Monocyte",
  exposure == "neu", "Neutrophil"
)][, outcome := fcase(
  outcome == "overall", "Overall",
  outcome == "colon", "Colon",
  outcome == "distal", "Distal",
  outcome == "female", "Female",
  outcome == "male", "Male",
  outcome == "proximal", "Proximal",
  outcome == "rectal", "Rectal"
)]
mvmr_all_res
mvmr_all_res[method == "IVW MVMR (direct)", method := "IVW MVMR"]
mvmr_all_res <- mvmr_all_res[method != "Observational MV",]
mvmr_all_res$outcome <- factor(mvmr_all_res$outcome, 
                             levels=c('Overall','Male','Female','Colon','Proximal','Distal','Rectal'))
mvmr_all_res$method <- factor(mvmr_all_res$method, 
                               levels=rev(c('IVW UVMR', 'IVW MVMR', 'Observational')))

custom_palette <- c(`Observational` = "#323D5A", `IVW MVMR` = "#00C0C0", `IVW UVMR` = "#FF4A4F")

# overall CRC only, for paper
all_tf_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 1.8), xtickbreaks = seq(0.0, 1.8, 0.1)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "Method", shape = "Method") + facet_grid(. ~ outcome) + 
  theme(legend.position="top") + guides(color=guide_legend(nrow=2, byrow=FALSE, reverse = TRUE), shape=guide_legend(nrow=2, byrow=FALSE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette)

ggsave(filename = "#.png", 
       plot = all_tf_mvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

# Overall, male, female
all_tf_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall", "Male", "Female") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 2.5), xtickbreaks = seq(0.0, 2.0, 0.4)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "Method", shape = "Method") + facet_grid(. ~ outcome) + 
  theme(legend.position="top") + guides(color=guide_legend(nrow=2, byrow=FALSE, reverse = TRUE), shape=guide_legend(nrow=2, byrow=FALSE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette) 

ggsave(filename = "#.png", 
       plot = all_tf_mvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

# Second do the overall, colon, proximal, distal, rectal
all_sp_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall", "Colon", "Proximal", "Distal", "Rectal") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 2.4), xtickbreaks = seq(0.0, 2.0, 0.4)
) + theme(axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "Method", shape = "Method") + facet_grid(. ~ outcome) + 
  theme(legend.position="top") + guides(color=guide_legend(nrow=2, byrow=TRUE, reverse = TRUE), shape=guide_legend(nrow=2, byrow=TRUE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette) 

ggsave(filename = "#.png", 
       plot = all_sp_mvmr,
       scale = 1, width = 270,
       device = "png",
       height = 210, units = c("mm"),
       dpi = 250, bg = "white")


####### Observational last
# All figures, but with observational as last analysis
# These have the legend on the right
# overall CRC only, for paper
all_tf_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 1.6), xtickbreaks = seq(0.0, 1.6, 0.1)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.text.x=element_text(angle=45,hjust=1)) +
  labs(colour = "Method", shape = "Method") + facet_grid(. ~ outcome) + 
  theme(legend.position="right") + guides(color=guide_legend(nrow=3, byrow=FALSE, reverse = TRUE), shape=guide_legend(nrow=3, byrow=FALSE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette)

ggsave(filename = "#.png", 
       plot = all_tf_mvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

# Overall, male, female
all_tf_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall", "Male", "Female") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 2.5), xtickbreaks = seq(0.0, 2.0, 0.4)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.text.x=element_text(angle=45,hjust=1)) +
  labs(colour = "Method", shape = "Method") + facet_grid(. ~ outcome) + 
  theme(legend.position="right") + guides(color=guide_legend(nrow=3, byrow=FALSE, reverse = TRUE), shape=guide_legend(nrow=3, byrow=FALSE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette)

ggsave(filename = "#.png", 
       plot = all_tf_mvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

# Second do the overall, colon, proximal, distal, rectal
all_sp_mvmr <- forestplot_mod(
  df = mvmr_all_res[outcome %in% c("Overall", "Colon", "Proximal", "Distal", "Rectal") & method != "Observational MV",], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 2.4), xtickbreaks = seq(0.0, 2.0, 0.4)
) + theme(axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.text.x=element_text(angle=45,hjust=1)) +
  labs(colour = "Method", shape = "Method") + facet_grid(. ~ outcome) + 
  theme(legend.position="right") + guides(color=guide_legend(nrow=3, byrow=FALSE, reverse = TRUE), shape=guide_legend(nrow=3, byrow=FALSE, reverse = TRUE)) +
  scale_color_manual(values=custom_palette) + scale_fill_manual(values=custom_palette)

ggsave(filename = "#.png", 
       plot = all_sp_mvmr,
       scale = 1, width = 270,
       device = "png",
       height = 210, units = c("mm"),
       dpi = 250, bg = "white")



