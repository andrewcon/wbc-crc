# Load packages------
load_pkg <- rlang::quos(tidyverse, data.table, ggplot2, ggpubr, ggrepel, RColorBrewer, gghighlight, foreach, doParallel, RNOmni, flextable)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))
source("#")

# Load data------
load("#.RData")

set_flextable_defaults(
  font.size = 11, font.family = "Arial",
  padding = 4)

### Create RIN of residuals from a lm for each WBC count
# FIRST, LOG10 of each WBC count
# 1) Adjust for sex, age, age-squared, the first 10 principal components, and other cohort-specific covariates
# 2) Additionally adjusting for BMI, smoking status, drinker status and townsendDI

# Double check that residuals are ordered to the same rows in dataset
m <- lm(wbc_count ~ age + neu_count, data=covars)
cor((predict(m) + residuals(m)), covars$wbc_count) 
cor(residuals(m), covars$wbc_count) 

# Log10 of all WBC count
covars$bas_log <- log(covars$bas_count + 1)
covars$eos_log <- log(covars$eos_count + 1)
covars$lym_log <- log(covars$lym_count + 1)
covars$mon_log <- log(covars$mon_count + 1)
covars$neu_log <- log(covars$neu_count + 1)
covars$wbc_log <- log(covars$wbc_count + 1)
shapiro.dt(covars$bas_log, "BC_log")
shapiro.dt(covars$eos_log, "EC_log")
shapiro.dt(covars$lym_log, "LC_log")
shapiro.dt(covars$mon_log, "MC_log")
shapiro.dt(covars$neu_log, "NC_log")
shapiro.dt(covars$wbc_log, "WC_log")


## 1)
covars$bas_residuals_rin <- RankNorm( residuals( lm(bas_log ~ sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$eos_residuals_rin <- RankNorm( residuals( lm(eos_log ~ sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$lym_residuals_rin <- RankNorm( residuals( lm(lym_log ~ sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$mon_residuals_rin <- RankNorm( residuals( lm(mon_log ~ sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$neu_residuals_rin <- RankNorm( residuals( lm(neu_log ~ sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$wbc_residuals_rin <- RankNorm( residuals( lm(wbc_log ~ sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
shapiro.dt(covars$bas_residuals_rin, "BC_residuals_rin")
shapiro.dt(covars$eos_residuals_rin, "EC_residuals_rin")
shapiro.dt(covars$lym_residuals_rin, "LC_residuals_rin")
shapiro.dt(covars$mon_residuals_rin, "MC_residuals_rin")
shapiro.dt(covars$neu_residuals_rin, "NC_residuals_rin")
shapiro.dt(covars$wbc_residuals_rin, "WC_residuals_rin")
covars$bas_residuals_rin_mv <- RankNorm( residuals( lm(bas_log ~ wbc_log + sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$eos_residuals_rin_mv <- RankNorm( residuals( lm(eos_log ~ wbc_log + sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$lym_residuals_rin_mv <- RankNorm( residuals( lm(lym_log ~ wbc_log + sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$mon_residuals_rin_mv <- RankNorm( residuals( lm(mon_log ~ wbc_log + sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$neu_residuals_rin_mv <- RankNorm( residuals( lm(neu_log ~ wbc_log + sex + age + age^2 + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
# Do the same, but divided by sex
# Univariable
covars_female <- covars[sex == "Female", ]
covars_male <- covars[sex == "Male", ]
covars_female$bas_female_residuals_rin <- RankNorm(residuals(lm(bas_log ~ age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$eos_female_residuals_rin <- RankNorm(residuals(lm(eos_log ~ age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$lym_female_residuals_rin <- RankNorm(residuals(lm(lym_log ~ age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$mon_female_residuals_rin <- RankNorm(residuals(lm(mon_log ~ age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$neu_female_residuals_rin <- RankNorm(residuals(lm(neu_log ~ age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$wbc_female_residuals_rin <- RankNorm(residuals(lm(wbc_log ~ age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_male$bas_male_residuals_rin <- RankNorm(residuals(lm(bas_log ~ age + age^2 + assessment_center +  sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$eos_male_residuals_rin <- RankNorm(residuals(lm(eos_log ~ age + age^2 + assessment_center +  sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$lym_male_residuals_rin <- RankNorm(residuals(lm(lym_log ~ age + age^2 + assessment_center +  sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$mon_male_residuals_rin <- RankNorm(residuals(lm(mon_log ~ age + age^2 + assessment_center +  sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$neu_male_residuals_rin <- RankNorm(residuals(lm(neu_log ~ age + age^2 + assessment_center +  sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$wbc_male_residuals_rin <- RankNorm(residuals(lm(wbc_log ~ age + age^2 + assessment_center +  sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
# Multivariable
covars_female$bas_female_residuals_rin_mv <- RankNorm(residuals(lm(bas_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$eos_female_residuals_rin_mv <- RankNorm(residuals(lm(eos_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$lym_female_residuals_rin_mv <- RankNorm(residuals(lm(lym_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$mon_female_residuals_rin_mv <- RankNorm(residuals(lm(mon_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$neu_female_residuals_rin_mv <- RankNorm(residuals(lm(neu_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_male$bas_male_residuals_rin_mv <- RankNorm(residuals(lm(bas_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$eos_male_residuals_rin_mv <- RankNorm(residuals(lm(eos_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$lym_male_residuals_rin_mv <- RankNorm(residuals(lm(lym_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$mon_male_residuals_rin_mv <- RankNorm(residuals(lm(mon_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$neu_male_residuals_rin_mv <- RankNorm(residuals(lm(neu_log ~ wbc_log + age + age^2 + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)

## 2)
covars$bas_residuals_rin2 <- RankNorm( residuals( lm(bas_log ~ sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$eos_residuals_rin2 <- RankNorm( residuals( lm(eos_log ~ sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$lym_residuals_rin2 <- RankNorm( residuals( lm(lym_log ~ sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$mon_residuals_rin2 <- RankNorm( residuals( lm(mon_log ~ sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$neu_residuals_rin2 <- RankNorm( residuals( lm(neu_log ~ sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$wbc_residuals_rin2 <- RankNorm( residuals( lm(wbc_log ~ sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$bas_residuals_rin_mv2 <- RankNorm( residuals( lm(bas_log ~ wbc_log + sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$eos_residuals_rin_mv2 <- RankNorm( residuals( lm(eos_log ~ wbc_log + sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$lym_residuals_rin_mv2 <- RankNorm( residuals( lm(lym_log ~ wbc_log + sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$mon_residuals_rin_mv2 <- RankNorm( residuals( lm(mon_log ~ wbc_log + sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
covars$neu_residuals_rin_mv2 <- RankNorm( residuals( lm(neu_log ~ wbc_log + sex + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_device + sample_year + sample_month + sample_day + sample_day_minutes + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=covars ) ) , k = 0.375)
# Do the same, but divided by sex
# Univariable
covars_female$bas_female_residuals_rin2 <- RankNorm(residuals(lm(bas_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$eos_female_residuals_rin2 <- RankNorm(residuals(lm(eos_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$lym_female_residuals_rin2 <- RankNorm(residuals(lm(lym_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$mon_female_residuals_rin2 <- RankNorm(residuals(lm(mon_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$neu_female_residuals_rin2 <- RankNorm(residuals(lm(neu_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$wbc_female_residuals_rin2 <- RankNorm(residuals(lm(wbc_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_male$bas_male_residuals_rin2 <- RankNorm(residuals(lm(bas_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$eos_male_residuals_rin2 <- RankNorm(residuals(lm(eos_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$lym_male_residuals_rin2 <- RankNorm(residuals(lm(lym_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$mon_male_residuals_rin2 <- RankNorm(residuals(lm(mon_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$neu_male_residuals_rin2 <- RankNorm(residuals(lm(neu_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$wbc_male_residuals_rin2 <- RankNorm(residuals(lm(wbc_log ~ age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
# Multivariable
covars_female$bas_female_residuals_rin_mv2 <- RankNorm(residuals(lm(bas_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$eos_female_residuals_rin_mv2 <- RankNorm(residuals(lm(eos_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$lym_female_residuals_rin_mv2 <- RankNorm(residuals(lm(lym_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$mon_female_residuals_rin_mv2 <- RankNorm(residuals(lm(mon_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_female$neu_female_residuals_rin_mv2 <- RankNorm(residuals(lm(neu_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_female)), k = 0.375)
covars_male$bas_male_residuals_rin_mv2 <- RankNorm(residuals(lm(bas_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$eos_male_residuals_rin_mv2 <- RankNorm(residuals(lm(eos_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$lym_male_residuals_rin_mv2 <- RankNorm(residuals(lm(lym_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$mon_male_residuals_rin_mv2 <- RankNorm(residuals(lm(mon_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)
covars_male$neu_male_residuals_rin_mv2 <- RankNorm(residuals(lm(neu_log ~ wbc_log + age + age^2 + BMI + townsendDI + smoking_status + drinker_status + assessment_center + sample_year + sample_month + sample_day + sample_day_minutes + sample_device + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covars_male)), k = 0.375)


### Logistic regression
## WBC + covariates on CRC risk
analysis_plan <- data.table(
    exposure = c(rep(c("bas", "eos", "lym", "mon", "neu", "wbc"), 3), rep(c("bas", "eos", "lym", "mon", "neu"), 3)),
    outcome = c(rep("overall", 6), rep("male", 6), rep("female", 6), rep("overall", 5), rep("male", 5), rep("female", 5)),
    method = c(rep("univariable", 18), rep("multivariable (vs. WBC count)", 15))
)

gecco_obs <- function (x) {
    my_exposure <- analysis_plan[[x,1]]
    my_outcome <- analysis_plan[[x,2]]
    my_method <- analysis_plan[[x,3]]
    if(my_outcome == "overall") {
        if(my_method == "univariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_residuals_rin"))
            my_glm <- glm(my_formula, data = covars, family = "binomial")
        } else if (my_method == "multivariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_residuals_rin_mv"))
            my_glm <- glm(my_formula, data = covars, family = "binomial")
        }
    } else if (my_outcome == "male") {
        if(my_method == "univariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_male_residuals_rin"))
            my_glm <- glm(my_formula, data = covars_male, family = "binomial")
        } else if (my_method == "multivariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_male_residuals_rin_mv"))
            my_glm <- glm(my_formula, data = covars_male, family = "binomial")
        }
    } else if (my_outcome == "female") {
        if(my_method == "univariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_female_residuals_rin"))
            my_glm <- glm(my_formula, data = covars_female, family = "binomial")
        } else if (my_method == "multivariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_female_residuals_rin_mv"))
            my_glm <- glm(my_formula, data = covars_female, family = "binomial")
        }
    }
    my_glm_summmary <- as.data.table(summary(my_glm)[[12]], keep.rownames=TRUE)
    my_result_dt <- data.table(
        exposure = my_exposure,
        outcome = my_outcome,
        model = "Model 1",
        type = my_method,
        method = "Observational",
        b = my_glm_summmary[[2, 2]],
        se = my_glm_summmary[[2, 3]],
        pval = my_glm_summmary[[2, 5]],
        exp_var_name = my_glm_summmary[[2, 1]]
    )
    return(my_result_dt)
}

# 1) Univariable analysis 
# 2) Multivariable analysis with log of immune cell counts
# 3) Compile results in a data.table
# Run the Observational analysis------
myCluster <- makeCluster(6, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

analyses_list <- foreach(x = 1:33) %dopar% {
  gecco_obs(
    x = x
  )
}

analyses_list <- rbindlist(analyses_list)

stopCluster(myCluster)

analyses_list$or <- exp(analyses_list$b)
analyses_list$or_lci <- exp(analyses_list$b - 1.96*analyses_list$se)
analyses_list$or_uci <- exp(analyses_list$b + 1.96*analyses_list$se)

### Now do it for Model 2
gecco_obs2 <- function (x) {
    my_exposure <- analysis_plan[[x,1]]
    my_outcome <- analysis_plan[[x,2]]
    my_method <- analysis_plan[[x,3]]
    if(my_outcome == "overall") {
        if(my_method == "univariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_residuals_rin2"))
            my_glm <- glm(my_formula, data = covars, family = "binomial")
        } else if (my_method == "multivariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_residuals_rin_mv2"))
            my_glm <- glm(my_formula, data = covars, family = "binomial")
        }
    } else if (my_outcome == "male") {
        if(my_method == "univariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_male_residuals_rin2"))
            my_glm <- glm(my_formula, data = covars_male, family = "binomial")
        } else if (my_method == "multivariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_male_residuals_rin_mv2"))
            my_glm <- glm(my_formula, data = covars_male, family = "binomial")
        }
    } else if (my_outcome == "female") {
        if(my_method == "univariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_female_residuals_rin2"))
            my_glm <- glm(my_formula, data = covars_female, family = "binomial")
        } else if (my_method == "multivariable") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ ", my_exposure, "_female_residuals_rin_mv2"))
            my_glm <- glm(my_formula, data = covars_female, family = "binomial")
        }
    }
    my_glm_summmary <- as.data.table(summary(my_glm)[[12]], keep.rownames=TRUE)
    my_result_dt <- data.table(
        exposure = my_exposure,
        outcome = my_outcome,
        model = "Model 2",
        type = my_method,
        method = "Observational",
        b = my_glm_summmary[[2, 2]],
        se = my_glm_summmary[[2, 3]],
        pval = my_glm_summmary[[2, 5]],
        exp_var_name = my_glm_summmary[[2, 1]]
    )
    return(my_result_dt)
}

# 1) Univariable analysis 
# 2) Multivariable analysis with log of immune cell counts
# 3) Compile results in a data.table
# Run the Observational analysis------
myCluster <- makeCluster(6, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

analyses_list2 <- foreach(x = 1:33) %dopar% {
  gecco_obs2(
    x = x
  )
}

analyses_list2 <- rbindlist(analyses_list2)

stopCluster(myCluster)

analyses_list2$or <- exp(analyses_list2$b)
analyses_list2$or_lci <- exp(analyses_list2$b - 1.96*analyses_list2$se)
analyses_list2$or_uci <- exp(analyses_list2$b + 1.96*analyses_list2$se)

# Multivariable when adjusting for inbetween WBCs
analysis_plan_between <- data.table(
    outcome = c("overall", "male", "female", "overall", "male", "female"),
    method = c(rep("multivariable_between", 6)),
    model = c(rep("Model 1", 3), rep("Model 2", 3))
)

gecco_obs_between <- function (x) {
    my_outcome <- analysis_plan_between[[x,1]]
    my_method <- analysis_plan_between[[x,2]]
    my_model <- analysis_plan_between[[x,3]]
    if(my_outcome == "overall") {
        if(my_model == "Model 1") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ bas_residuals_rin + eos_residuals_rin + lym_residuals_rin + mon_residuals_rin + neu_residuals_rin"))
            my_glm <- glm(my_formula, data = covars, family = "binomial")
        } else if (my_model == "Model 2") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ bas_residuals_rin2 + eos_residuals_rin2 + lym_residuals_rin2 + mon_residuals_rin2 + neu_residuals_rin2"))
            my_glm <- glm(my_formula, data = covars, family = "binomial")
        }
    } else if (my_outcome == "male") {
        if(my_model == "Model 1") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ bas_male_residuals_rin + eos_male_residuals_rin + lym_male_residuals_rin + mon_male_residuals_rin + neu_male_residuals_rin"))
            my_glm <- glm(my_formula, data = covars_male, family = "binomial")
        } else if (my_model == "Model 2") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ bas_male_residuals_rin2 + eos_male_residuals_rin2 + lym_male_residuals_rin2 + mon_male_residuals_rin2 + neu_male_residuals_rin2"))
            my_glm <- glm(my_formula, data = covars_male, family = "binomial")
        }
    } else if (my_outcome == "female") {
        if(my_model == "Model 1") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ bas_female_residuals_rin + eos_female_residuals_rin + lym_female_residuals_rin + mon_female_residuals_rin + neu_female_residuals_rin"))
            my_glm <- glm(my_formula, data = covars_female, family = "binomial")
        } else if (my_model == "Model 2") {
            my_formula <- formula(paste0("incident_colorectal_cancer ~ bas_female_residuals_rin2 + eos_female_residuals_rin2 + lym_female_residuals_rin2 + mon_female_residuals_rin2 + neu_female_residuals_rin2"))
            my_glm <- glm(my_formula, data = covars_female, family = "binomial")
        }
    }
    my_glm_summmary <- as.data.table(summary(my_glm)[[12]], keep.rownames=TRUE)
    my_result_dt <- data.table(
        exposure = c("bas", "eos", "lym", "mon", "neu"),
        outcome = my_outcome,
        model = my_model,
        type = my_method,
        method = "Observational",
        b = my_glm_summmary$Estimate[2:6],
        se = my_glm_summmary$`Std. Error`[2:6],
        pval = my_glm_summmary$`Pr(>|z|)`[2:6],
        exp_var_name = my_glm_summmary$rn[2:6]
    )
    return(my_result_dt)
}

myCluster <- makeCluster(6, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

analyses_list3 <- foreach(x = 1:6) %dopar% {
  gecco_obs_between(
    x = x
  )
}

analyses_list3 <- rbindlist(analyses_list3)

stopCluster(myCluster)

analyses_list3$or <- exp(analyses_list3$b)
analyses_list3$or_lci <- exp(analyses_list3$b - 1.96*analyses_list3$se)
analyses_list3$or_uci <- exp(analyses_list3$b + 1.96*analyses_list3$se)


# Combine everything
analyses_list_all <- rbind(analyses_list, analyses_list2, analyses_list3)
fwrite(analyses_list_all, file="#.txt", quote=F, row.names=F, sep="\t", na=NA)

# Format as table and export
analyses_list_all <- analyses_list_all[, -c("b", "se", "exp_var_name")]
names(analyses_list_all) <- c("Exposure", "Outcome", "Model", "Analysis type", "Method", "P-value", "OR", "LCI", "UCI")
analyses_list_all <- analyses_list_all[, c("Exposure", "Outcome", "Analysis type", "Model", "OR", "LCI", "UCI", "P-value")]
analyses_list_all[, Exposure := fcase(
    Exposure == "bas", "Basophil", 
    Exposure == "eos", "Eosinophil", 
    Exposure == "lym", "Lymphocyte", 
    Exposure == "mon", "Monocyte", 
    Exposure == "neu", "Neutrophil", 
    Exposure == "wbc", "Overall WBC"
)]
analyses_list_all[, Outcome := str_to_title(Outcome)]
analyses_list_all[, `Analysis type` := str_to_title(`Analysis type`)]
analyses_list_all <- analyses_list_all %>% mutate_at(vars(OR, LCI, UCI), list(~ round(., 2)))
analyses_list_all[, `P-value` := format.pval(`P-value`)]
obs_ft <- flextable(analyses_list_all) %>% 
    theme_zebra() %>% 
    color(~ `P-value` < 0.05, ~ `P-value`, color = "red") %>%
    autofit()
save_as_docx(obs_ft, path = "#.docx")

