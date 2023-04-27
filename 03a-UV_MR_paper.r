#!/usr/bin/env Rscript

# Load relevant libraries------
list_of_packages <- c(
  "foreach",
  "doParallel",
  "kableExtra",
  "data.table",
  "TwoSampleMR",
  "ieugwasr",
  "tidyr",
  "rlang"
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

# Functions------
gecco_mr <- function(exp_path, exp_name, out_path, out_name, reverse = FALSE) {
  # Load data and format it
  if(!reverse) {
    exp_dat <- fread(exp_path, stringsAsFactors = F)
    out_dat <- fread(out_path, stringsAsFactors = F)
    exp_dat <- format_data(
      exp_dat, type="exposure", snp_col = "ID", beta_col = "beta", se_col = "se",
      eaf_col = "eaf", effect_allele_col = "reference_allele", other_allele_col = "other_allele", 
      pval_col = "p-value", chr_col = "CHR", pos_col = "BP", 
      samplesize_col = "n_samples"
    )
    # Perform clumping
    exp_dat <- clump_data(exp_dat)
    # Format outcome data
    out_dat <- format_data(
      out_dat, type="outcome", snp_col = "SNP", beta_col = "Effect", se_col = "StdErr",
      eaf_col = "Freq1", effect_allele_col = "Allele1", other_allele_col = "Allele2", 
      pval_col = "P.value", chr_col = "Chromosome", pos_col = "Position", 
      samplesize_col = "Neff", snps = exp_dat$SNP
    )
  }
  # Harmonize data
  dat <- harmonise_data(
    exposure_dat = exp_dat, 
    outcome_dat = out_dat
  )
  # Perform MR
  res <- mr(dat)
  res <- generate_odds_ratios(res)
  res$exposure <- exp_name
  res$outcome <- out_name
  het <- mr_heterogeneity(dat)
  het$exposure <- exp_name
  het$outcome <- out_name
  ple <- mr_pleiotropy_test(dat)
  ple$exposure <- exp_name
  ple$outcome <- out_name
  res_single <- mr_singlesnp(dat)
  res_single$exposure <- exp_name
  res_single$outcome <- out_name
  res_loo <- mr_leaveoneout(dat)
  res_loo$exposure <- exp_name
  res_loo$outcome <- out_name
  dir_test <- directionality_test(dat)
  dir_test$exposure <- exp_name
  dir_test$outcome <- out_name
  res_all <- list(res, het, ple, res_single, res_loo, dir_test, dat)
  print(paste0("Finished ", exp_name, " on ", out_name, "."))
  return(res_all)
}

# Assign path to WBC and CRC data------
wbc_paths <- c(list.files("#", full.names=T), list.files("#", full.names=T)) 
names(wbc_paths) <- rep(c("bas", "eos", "lym", "mon", "neu", "wbc_all"), 2)
crc_paths <- c(list.files("#", full.names = T), list.files("#", full.names=T))
names(crc_paths) <- rep(c("colon", "distal", "female", "left", 
                      "male", "overall", "proximal", "rectal"), 2)

# Create list which will contain all MR data------
mr_list <- list()
exposure <- expand_grid(wbc_paths[1:6], crc_paths[1:8])
outcome  <- expand_grid(wbc_paths[7:12], crc_paths[9:16])
names(exposure) <- c("wbc", "crc")
names(outcome) <- c("wbc", "crc") 
mr_info <- as.data.table(rbind(exposure, outcome))
mr_info$wbc_name <- rep(names(wbc_paths), each=8)
mr_info$crc_name <- rep(names(crc_paths), 6)
mr_info <- mr_info[, c(1,3,2,4)]
mr_info$reverse <- c(rep(FALSE, 48), rep(TRUE, 48))

# Run the MR analysis------
myCluster <- makeCluster(6, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

mr_list <- foreach(x = 1:96) %dopar% {
  gecco.mr(
    exp_path = mr_info[[x,1]],
    exp_name = mr_info[[x,2]],
    out_path = mr_info[[x,3]],
    out_name = mr_info[[x,4]],
    reverse = mr_info[[x,5]]
  )
}

stopCluster(myCluster)

save(mr_list, file = "#.RData")

### Table and figure generation------
# Load relevant libraries
list_of_packages <- c(
  "data.table",
  "TwoSampleMR",
  "ieugwasr",
  "ggplot2",
  "ggpubr",
  "tidyverse",
  "ggforestplot",
  "rlang"
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
rm(list = c("new_packages", "package.i", "list_of_packages"))

load("#.RData")
load("#.RData")
analyses_list_all <- fread("#.txt")

source("00_ggforestplot.r")

### Figures for the UV MR analysis
# Load main MR data
main_mr_dat <- lapply(mr_list, `[[`, 1)
main_mr_dat <- rbindlist(main_mr_dat)
main_mr_dat <- main_mr_dat[, c("exposure", "outcome", "method", "b", "se", "pval")]
main_mr_dat <- main_mr_dat[! exposure %in% c("overall", "colon", "distal", "female", "left", "male", "proximal", "rectal"), ]
presso_res <- rbindlist(presso_list, fill = TRUE)
presso_res <- presso_res[`MR Analysis` != "Raw",]
presso_res$method <- "MR PRESSO"
presso_res <- presso_res[, c("exposure", "outcome", "method", "Causal Estimate", "Sd", "P-value")]
setnames(presso_res, c("Causal Estimate", "Sd", "P-value"), c("b", "se", "pval"))
analyses_list_all <- analyses_list_all[model == "Model 1",]
analyses_list_all <- analyses_list_all[type == "univariable", c("exposure", "outcome", "method", "b", "se", "pval")]

all_mr_dat <- rbind(main_mr_dat, presso_res, analyses_list_all)
all_mr_dat <- all_mr_dat[order(exposure, outcome, method)]
all_mr_dat <- all_mr_dat[outcome != "left"]
all_mr_dat <- all_mr_dat[method != "Simple mode"]
all_mr_dat <- all_mr_dat[, exposure := fcase(
  exposure == "bas", "Basophil",
  exposure == "eos", "Eosinophil",
  exposure == "lym", "Lymphocyte",
  exposure == "mon", "Monocyte",
  exposure == "neu", "Neutrophil",
  exposure == "wbc_all" | exposure == "wbc", "Overall WBC count"
)][, outcome := fcase(
  outcome == "overall", "Overall",
  outcome == "colon", "Colon",
  outcome == "distal", "Distal",
  outcome == "female", "Female",
  outcome == "male", "Male",
  outcome == "proximal", "Proximal",
  outcome == "rectal", "Rectal"
)]
all_mr_dat$outcome <- factor(all_mr_dat$outcome, 
                             levels=c('Overall','Male','Female','Colon','Proximal','Distal','Rectal'))
all_mr_dat$method <- factor(all_mr_dat$method, 
                             levels=rev(c('Observational', 'Inverse variance weighted', 'MR Egger', 'Weighted median', 'Weighted mode', 'MR PRESSO')))

# Use Nightingale forestplot package
# First do the overall CRC for the paper
all_fp_uvmr <- forestplot_mod(
  df = all_mr_dat[outcome %in% c("Overall") & method %in% c('Observational', 'Inverse variance weighted') & !exposure %in% c('Overall WBC count'),], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 1.5), xtickbreaks = seq(0.0, 2.4, 0.4)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "MR method", shape = "MR method") + facet_grid(. ~ outcome) + 
  theme(legend.position="top", legend.justification = c(1,0))
  #+ theme(legend.position="top") + guides(colour=guide_legend(reverse = TRUE, nrow=2, byrow=TRUE), shape=guide_legend(reverse = TRUE, nrow=2,byrow=TRUE))

ggsave(filename = "#.png", 
       plot = all_fp_uvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

# Second do the overall CRC, male and female
all_fp_uvmr <- forestplot_mod(
  df = all_mr_dat[outcome %in% c("Overall", "Male", "Female"),], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 1.5), xtickbreaks = seq(0.0, 2.4, 0.4)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "MR method", shape = "MR method") + facet_grid(. ~ outcome) + 
  theme(legend.position="top", legend.justification = c(1,0))
  #+ theme(legend.position="top") + guides(colour=guide_legend(reverse = TRUE, nrow=2, byrow=TRUE), shape=guide_legend(reverse = TRUE, nrow=2,byrow=TRUE))

ggsave(filename = "#.png", 
       plot = all_fp_uvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

# Second do the overall, colon, proximal, distal, rectal
all_mr_dat_site <- all_mr_dat[outcome %in% c("Overall", "Colon", "Proximal", "Distal", "Rectal") & method != "Observational",]
all_mr_dat_site$outcome <- factor(all_mr_dat_site$outcome, 
                             levels=c('Overall','Male','Female','Colon','Proximal','Distal','Rectal'))
all_mr_dat_site$method <- factor(all_mr_dat_site$method, 
                             levels=rev(c('Inverse variance weighted', 'MR Egger', 'Weighted median', 'Weighted mode', 'MR PRESSO')))

all_sp_uvmr <- forestplot_mod(
  df = all_mr_dat_site, name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for colorectal cancer (95% CI)
  per 1−SD increment in WBC count",
  xlim = c(0.6, 1.5), xtickbreaks = seq(0.0, 2.4, 0.4)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "MR method", shape = "MR method") + facet_grid(. ~ outcome) + 
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) + 
  theme(legend.position="top", legend.justification = c(1,0), legend.box="horizontal", legend.margin=margin())

ggsave(filename = "#.png", 
       plot = all_sp_uvmr,
       scale = 1, width = 210,
       device = "png",
       height = 250, units = c("mm"),
       dpi = 250, bg = "white")

### Table with het and ple MR results
## Concatenate with main MR results and generate a table
main_mr_dat <- rbindlist(lapply(mr_list, `[[`, 1))
main_mr_dat <- main_mr_dat[, c("exposure", "outcome", "method", "b", "se", "pval", "nsnp")]
main_mr_dat <- main_mr_dat[! exposure %in% c("overall", "colon", "distal", "female", "left", "male", "proximal", "rectal"), ]
presso_res <- rbindlist(presso_list, fill = TRUE)
presso_res <- presso_res[`MR Analysis` != "Raw",]
presso_res$method <- "MR PRESSO"
presso_res <- presso_res[, c("exposure", "outcome", "method", "Causal Estimate", "Sd", "P-value")]
setnames(presso_res, c("Causal Estimate", "Sd", "P-value"), c("b", "se", "pval"))
presso_res <- merge(presso_res, main_mr_dat[method == "Inverse variance weighted", c("exposure", "outcome", "nsnp")][1:48,], by=c("exposure", "outcome"))

all_mr_dat <- rbind(main_mr_dat, presso_res)
all_mr_dat <- all_mr_dat[order(exposure, outcome, method)]
all_mr_dat <- all_mr_dat[outcome != "left"]
all_mr_dat <- all_mr_dat[method != "Simple mode"]
all_mr_dat <- generate_odds_ratios(all_mr_dat)

all_mr_dat <- all_mr_dat[, exposure := fcase(
  exposure == "bas", "Basophil",
  exposure == "eos", "Eosinophil",
  exposure == "lym", "Lymphocyte",
  exposure == "mon", "Monocyte",
  exposure == "neu", "Neutrophil",
  exposure == "wbc_all" | exposure == "wbc", "Overall WBC count"
)][, outcome := fcase(
  outcome == "overall", "Overall",
  outcome == "colon", "Colon",
  outcome == "distal", "Distal",
  outcome == "female", "Female",
  outcome == "male", "Male",
  outcome == "proximal", "Proximal",
  outcome == "rectal", "Rectal"
)]

all_mr_dat <- all_mr_dat[, -c("b", "se", "lo_ci", "up_ci")]
all_mr_dat <- all_mr_dat[, c("exposure", "outcome", "method", "nsnp", "or", "or_lci95", "or_uci95", "pval")]
names(all_mr_dat) <- c("Exposure", "Outcome", "Method", "No SNPs", "OR", "LCI", "UCI", "P-value")
all_mr_dat <- all_mr_dat %>% mutate_at(vars(OR, LCI, UCI), list(~ round(., 2)))
all_mr_dat[, `P-value` := format.pval(`P-value`)]

set_flextable_defaults(
  font.size = 11, font.family = "Arial",
  padding = 4)
uvmr_ft <- flextable(all_mr_dat) %>% 
    theme_zebra() %>% 
    color(~ `P-value` < 0.05, ~ `P-value`, color = "red") %>%
    autofit()
save_as_docx(uvmr_ft, path = "#.docx")

## Make MR-PRESSO table
presso_res <- rbindlist(presso_list, fill = TRUE)
presso_res <- presso_res[`MR Analysis` != "Raw",]
presso_res$method <- "MR PRESSO"
presso_res <- merge(presso_res, main_mr_dat[method == "Inverse variance weighted", c("exposure", "outcome", "nsnp")][1:48,], by=c("exposure", "outcome"))
presso_res <- presso_res[outcome != "left",]
presso_res <- presso_res[, exposure := fcase(
  exposure == "bas", "Basophil",
  exposure == "eos", "Eosinophil",
  exposure == "lym", "Lymphocyte",
  exposure == "mon", "Monocyte",
  exposure == "neu", "Neutrophil",
  exposure == "wbc_all" | exposure == "wbc", "Overall WBC count"
)][, outcome := fcase(
  outcome == "overall", "Overall",
  outcome == "colon", "Colon",
  outcome == "distal", "Distal",
  outcome == "female", "Female",
  outcome == "male", "Male",
  outcome == "proximal", "Proximal",
  outcome == "rectal", "Rectal"
)]

# PRESSO summary table
presso_res_summary <- presso_res[, c("exposure", "outcome", "Causal Estimate", "Sd", "P-value", "presso_global_pval", "presso_distortion_pval", "nsnp", "outlier_nr")]
names(presso_res_summary) <- c("Exposure", "Outcome", "Beta", "SE", "P-value", "Global P-value", "Distortion P-value", "No SNPs", "No outliers")
presso_res_summary <- presso_res_summary %>% mutate_if(is.numeric, ~round(., 2))
presso_res_summary[, `Global P-value` := ifelse(`Global P-value` == "<0.000333333333333333", "<3.33e-4", as.character(round(as.numeric(`Global P-value`), 2)))]
presso_summary_ft <- flextable(presso_res_summary) %>% 
    theme_zebra()
save_as_docx(presso_summary_ft, path = "#.docx")

# PRESSO outlier table 
presso_res_outliers <- presso_res[, c("exposure", "outcome", "outlier_snp")]
presso_res_outliers <- as.data.table(separate_rows(presso_res_outliers, outlier_snp, sep = ";", convert = FALSE))
presso_res_outliers <- presso_res_outliers %>% 
  setcolorder(c("outlier_snp", "exposure", "outcome")) %>%
  arrange(outlier_snp, .by_group = FALSE) %>% 
  group_by(outlier_snp) %>% as.data.table()
presso_res_outliers[presso_res_outliers=="" | presso_res_outliers=="NA"] <- NA
presso_res_outliers <- na.omit(presso_res_outliers)
presso_res_outliers[, exposure_outcome := paste0(exposure, "_", outcome)]

newdf <- dcast(presso_res_outliers, outlier_snp ~ exposure,
               value.var = "exposure", fun.aggregate = length)
newdf[newdf==0] <- NA
setnames(newdf, "outlier_snp", "Outlier SNP")

newdf_ft <- flextable(newdf)

colourer <- scales::col_numeric(
  palette = c("transparent", "red"),
  na.color = "#FFFFFF",
  domain = c(0, 10))

newdf_ft <- newdf_ft %>% bg(
    bg = colourer,
    j = ~ . -`Outlier SNP`,
    part = "body") %>% 
    theme_vanilla() %>%
    autofit()
save_as_docx(newdf_ft, path = "#.docx")

## Combine het, ple and dir_test
mr_het <- rbindlist(lapply(mr_list, `[[`, 2))[1:96,]
mr_ple <- rbindlist(lapply(mr_list, `[[`, 3))[1:48,]
mr_dir_test <- rbindlist(lapply(mr_list, `[[`, 6))[1:48,]

mr_sens <- merge(mr_het, mr_ple, by=c("exposure", "outcome"))
mr_sens <- merge(mr_sens, mr_dir_test, by=c("exposure", "outcome"))
mr_sens <- mr_sens[, c("exposure", "outcome", "method", "Q_pval", "egger_intercept", "pval", "correct_causal_direction", "steiger_pval")]
mr_sens <- mr_sens[outcome != "left",]
names(mr_sens) <- c("Exposure", "Outcome", "Het method", "Het P", "Ple intercept", "Ple P", "Correct direction", "Steiger P")
mr_sens <- mr_sens[, Exposure := fcase(
  Exposure == "bas", "Basophil",
  Exposure == "eos", "Eosinophil",
  Exposure == "lym", "Lymphocyte",
  Exposure == "mon", "Monocyte",
  Exposure == "neu", "Neutrophil",
  Exposure == "wbc_all" | Exposure == "wbc", "Overall WBC count"
)][, Outcome := fcase(
  Outcome == "overall", "Overall",
  Outcome == "colon", "Colon",
  Outcome == "distal", "Distal",
  Outcome == "female", "Female",
  Outcome == "male", "Male",
  Outcome == "proximal", "Proximal",
  Outcome == "rectal", "Rectal"
)][, `Het method` := fcase(
  `Het method` == "MR Egger", "MR Egger",
  `Het method` == "Inverse variance weighted", "IVW"
)]

mr_sens <- mr_sens %>% mutate_at(vars(`Het P`, `Ple P`, `Ple intercept`, `Steiger P`), list(~ format.pval(.,)))

mr_sens_ft <- flextable(mr_sens) %>% 
    theme_zebra() %>% 
    color(~ `Het P` < 0.05, ~ `Het P`, color = "orange") %>%
    color(~ `Ple P` < 0.05, ~ `Ple P`, color = "red") %>% 
    autofit() %>% set_table_properties(layout = "autofit")
mr_sens_ft <- height(mr_sens_ft, height = .5)
save_as_docx(mr_sens_ft, path = "#.docx")

## Work with single and loo SNP data
mr_single <- rbindlist(lapply(mr_list, `[[`, 4))
mr_single <- mr_single[exposure %in% c("bas", "eos", "lym", "mon", "neu", "wbc_all"), ]
mr_single <- mr_single[outcome != "left",]
mr_single <- mr_single[, exposure := fcase(
  exposure == "bas", "Basophil",
  exposure == "eos", "Eosinophil",
  exposure == "lym", "Lymphocyte",
  exposure == "mon", "Monocyte",
  exposure == "neu", "Neutrophil",
  exposure == "wbc_all" | exposure == "wbc", "Overall WBC count"
)][, outcome := fcase(
  outcome == "overall", "Overall",
  outcome == "colon", "Colon",
  outcome == "distal", "Distal",
  outcome == "female", "Female",
  outcome == "male", "Male",
  outcome == "proximal", "Proximal",
  outcome == "rectal", "Rectal"
)]
mr_single[, exposure_outcome := paste0(exposure, "_", outcome)]

for(x in unique(mr_single$exposure_outcome)) {
  p2 <- mr_forest_plot(mr_single[exposure_outcome == x,])
  p2 <- p2[[1]] + labs_pubr() + theme(axis.text.y = element_text(size = 5)) 
  ggsave(filename = paste0("#", x, "_singlesnp_res.png"), 
        plot = p2,
        scale = 1, width = 210,
        device = "png",
        height = 280, units = c("mm"),
        dpi = 250, bg = "white")
}

