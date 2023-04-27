# This script is made to work on an interactive terminal

library("data.table")
library("reshape")
library("stringr")
library("tidyr")
library("dplyr")
library("lubridate")
library("RNOmni")
library("ggplot2")
library("ggpubr")
library("flextable")
library("officer")
library("broom")
library("psych")
library("janitor")
library("gtsummary")
library("ggcorrplot")
source("00-UKBB_functions.r")

set_flextable_defaults(
  font.size = 11, font.family = "Arial",
  padding = 4)

covars <- fread("#", stringsAsFactors = F, na.strings=c(""," ","NA"))
info_vec <- readLines("#")

### Create descriptive statistics on all WBC variables
ds_wbc <- covars %>% 
        select(bas_count, eos_count, lym_count, mon_count, neu_count, wbc_count) %>% 
        summarise_all(describe) %>% pivot_longer(cols = everything(),
                names_sep = "_",
                names_to  = c("variable", ".value")) %>% as.data.table %>% 
 mutate_if(is.numeric, round)
ds_wbc[, variable := c("Basophil count", "Eosinophil count", "Lymphocyte count", "Monocyte count", "Neutrophil count", "Overall WBC count")]
colnames(ds_wbc) <- c("Trait", "cv", "N", "Mean", "SD", "Median", "Trimmed", "MAD", "Min", "Max", "Range", "Skew", "Kurtosis", "SE")
ds_wbc <- ds_wbc[, -c("cv", "Trimmed", "SE", "N")]
ds_wbc_ft <- flextable(ds_wbc) %>% autofit()
save_as_docx(ds_wbc_ft, path = "#.docx")

### Create descriptive statistics on all relevant variables

# First convert all relevant variables
covars$incident_colorectal_cancer <- factor(covars$incident_colorectal_cancer, 
                                levels = c(1, 2), 
                                labels = c("Control", "Case")
)
covars$sex <- factor(covars$sex, 
                     levels = c(0, 1), 
                     labels = c("Female", "Male")
)
covars$smoking_status <- factor(covars$smoking_status, 
                                levels = c(0, 1, 2), 
                                labels = c("Never", "Previous", "Current")
)
covars$drinker_status <- factor(covars$drinker_status, 
                                levels = c(0, 1, 2), 
                                labels = c("Never", "Previous", "Current")
)

# Create table
ds_vars <- covars %>% 
  select(incident_colorectal_cancer, sex, age, BMI, smoking_status, drinker_status, townsendDI, bas_count, eos_count, lym_count, mon_count, neu_count, wbc_count) %>% 
  tbl_summary(     
    by = incident_colorectal_cancer,                                               # stratify entire table by outcome
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = all_categorical() ~ "categorical",                 # force all categorical levels to display
    label  = list(                                              # display labels for column names
      incident_colorectal_cancer ~ "Outcome",                           
      sex ~ "Sex",
      age ~ "Age (years)",
      BMI ~ "Body mass index",
      smoking_status ~ "Smoking status",
      drinker_status ~ "Alcohol drinker status",
      townsendDI ~ "Townsend deprivation index",
      bas_count ~ "Basophil count",
      eos_count ~ "Eosinophil count",
      lym_count ~ "Lymphocyte count",
      mon_count ~ "Monocyte count",
      neu_count ~ "Neutrophil count",
      wbc_count  ~ "Overall WBC count"),
    missing_text = "Missing"
  )
ds_vars_ft <- as_flex_table(ds_vars) %>% autofit()
save_as_docx(ds_vars_ft, path = "#.docx")

### Create correlation table and matrix for WBC count
corr <- covars %>% 
  select(bas_count, eos_count, lym_count, mon_count, neu_count) %>% 
  rename_at(1:5, ~ c("Basophil", "Eosinophil", "Lymphocyte", "Monocyte", "Neutrophil")) %>%
  cor %>% round(1)

p.mat <- covars %>% 
  select(bas_count, eos_count, lym_count, mon_count, neu_count) %>% 
  rename_at(1:5, ~ c("Basophil", "Eosinophil", "Lymphocyte", "Monocyte", "Neutrophil")) %>%
  cor_pmat

wbc_corrplot <- ggcorrplot(corr, hc.order = FALSE,
    type = "lower", show.diag=TRUE, p.mat = p.mat, lab = TRUE)
ggsave(filename="#.png", 
  plot = wbc_corrplot, device = "png", width = 210, height = 210, units = "mm")

### Variance explained by different variables
covars$sample_year <- factor(covars$sample_year, 
                             levels = c(2007, 2008, 2009, 2010), 
                             labels = c(2007, 2008, 2009, 2010)
)
covars$sample_month <- factor(covars$sample_month, 
                              levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
                              labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
)

wbc_traits <- c("bas_count", "eos_count", "lym_count", "mon_count", "neu_count", "wbc_count")
names(wbc_traits) <- c("Basophil", "Eosinophil", "Lymphocyte", "Monocyte", "Neutrophil", "Overall WBC")
wbc_anova_df <- list()

for(z in 1:length(wbc_traits)) {
    lm_formula <- c("assessment_center", "sample_device", "sample_year", "sample_month", "sample_day", "sample_day_minutes", 
        "sex", "age", paste0("PC", 1:10),
        "BMI", "townsendDI", "smoking_status", "drinker_status")
    lm_formula <- paste0(lm_formula, collapse = " + ")
    lm_formula <- as.formula(paste0(wbc_traits[z], " ~ ", lm_formula))
    covars_lm <- lm(lm_formula, data = covars)
    # Univariate analysis
    covars.univariate <- data.table()
    x <- c("assessment_center", "sample_device", "sample_year", "sample_month", "sample_day", "sample_day_minutes", 
        "sex", "age", paste0("PC", 1:10),
        "BMI", "townsendDI", "smoking_status", "drinker_status")
    for(v in x) {
        temp_formula <- as.formula(paste0(wbc_traits[z], " ~ ", v))
        temp_lm <- lm(temp_formula, data = covars)
        temp_aov <- aov(temp_lm)
        temp_summary <- summary(temp_aov)
        afss <- sum(temp_summary[[1]][[2]])
        temp.univariate <- data.table(
            "rn" = v, 
            "DF" = temp_summary[[1]][1,1],
            "Sum Sq" = temp_summary[[1]][1,2],
            "Mean Sq" = temp_summary[[1]][1,3],
            "F value" = temp_summary[[1]][1,4],
            "Pr(>F)" = temp_summary[[1]]$'Pr(>F)'[1],
            "PctExp" = temp_summary[[1]][1,2]/afss*100,
            "type" = "univariate"
        )
        covars.univariate <- rbind(covars.univariate, temp.univariate)
    }
    # Type I ANOVA
    covars.I.aov <- aov(covars_lm)
    covars.I.aov <- summary(covars.I.aov)[[1]]
    covars.I.aov <- setDT(covars.I.aov, keep.rownames = TRUE)[]
    afss <- covars.I.aov$"Sum Sq"
    covars.I.aov[, PctExp := afss/sum(afss)*100]
    covars.I.aov$type <- "anova_one"
    covars.I.aov[, "rn" := lapply(.SD, trimws), .SDcols = "rn"]
    # Type II ANOVA
    covars.II.aov <- car::Anova(covars_lm, type = 2)
    covars.II.aov <- setDT(covars.II.aov, keep.rownames = TRUE)[]
    afss <- covars.II.aov$"Sum Sq"
    covars.II.aov[, PctExp := afss/sum(afss)*100]
    covars.II.aov$type <- "anova_two"
    covars.II.aov$'Mean Sq' <- NA
    covars.II.aov <- covars.II.aov[, c("rn", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)", "PctExp", "type")]
    # Merge results
    covars.all.aov <- rbind(covars.univariate, covars.I.aov, covars.II.aov, fill = TRUE)
    covars.all.aov$PctExp <- signif(covars.all.aov$PctExp, digits = 3)
    covars.all.aov <- covars.all.aov[!startsWith(rn, "Residuals")]
    covars.all.aov <- covars.all.aov[dim(covars.all.aov)[1]:1,]
    covars.all.aov$type <- factor(covars.all.aov$type, 
                                levels = c("anova_two", "anova_one", "univariate"),
                                labels = c("ANOVA II", "ANOVA I (top-bottom order)", "Univariate"))
    covars.all.aov$trait <- names(wbc_traits)[z]
    covars.all.aov <- covars.all.aov[, c(ncol(covars.all.aov), 1:(ncol(covars.all.aov) - 1)), with = F]
    wbc_anova_df[[z]] <- covars.all.aov
}

wbc_anova_df <- rbindlist(wbc_anova_df)

# Plot variance explained
# Define scale
aov_scale_manual <- scale_fill_manual(
  "Analysis type",
  values = c("#4DAF4A", "#377EB8", "#E41A1C"),
  guide = guide_legend(reverse = TRUE)
)
aov_scale_x <- scale_x_discrete(
  breaks = c("assessment_center", "sample_device", "sample_year", "sample_month", "sample_day", "sample_day_minutes", 
        "sex", "age", paste0("PC", 1:10),
        "BMI", "townsendDI", "smoking_status", "drinker_status"),
  labels = c("Assessment center", "Sampling device", "Sample year", "Sample month", "Sample day", "Sample day\nminutes",
             "Sex", "Age", 
             "PC1", "PC2", "PC3", "PC4", "PC5",
             "PC6", "PC7", "PC8", "PC9", "PC10",
             "Body mass index", "Townsend\ndepravation index", "Smoking status", "Drinker status")
)
v_barplot <- ggbarplot(
  data = wbc_anova_df,
  x = "rn",
  y = "PctExp",
  combine = TRUE,
  facet.by = "trait",
  color = "black",
  fill = "type",
  palette = NULL,
  orientation = "horiz",
  size = 0.25,
  width = 0.8,
  title = NULL,
  xlab = "",
  ylab = "Variance explained %",
  label = TRUE,
  lab.col = "black",
  lab.size = 2,
  lab.pos = c("out"),
  lab.vjust = 0.6,
  lab.hjust = -0.2,
  lab.nb.digits = 2,
  sort.by.groups = TRUE,
  position = position_dodge(0.9)
) + aov_scale_manual + aov_scale_x +
  theme(axis.text.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.background = element_rect(colour = "black", size = 0.5)) + 
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) + 
  guides(colour = guide_legend(reverse=TRUE))

ggsave(filename="#.png", 
  plot = v_barplot, device = "png", width = 210, height = 297, units = "mm")

save(covars, file = "#.RData")
