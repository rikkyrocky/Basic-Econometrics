#-----------------------------------------------
# Load libraries
#-----------------------------------------------
library(tidyverse)
library(naniar)
library(stargazer)
library(lmtest)
library(dplyr)

#-----------------------------------------------
# Import raw data
#-----------------------------------------------
bw <- read.csv("births_tidy_2019.csv")

# Summary stats of raw data
stargazer(bw, type = "text", title = "Raw Data Summary Statistics")

#-----------------------------------------------
# Replace "unknown" codes with NA
#-----------------------------------------------
bw_na <- bw %>%
  replace_with_na(replace = list(
    dob_tt      = 9999,
    bfacil3     = 3,
    dmar        = 3,
    mbstate_rec = 3,
    meduc       = 9,
    fagecomb    = 99,
    fagerec11   = 11,
    feduc       = 9,
    priordead   = 99,
    priorlive   = 99,
    priorterm   = 99,
    lbo_rec     = 9,
    tbo_rec     = 9,
    illb_r      = c(888, 999),
    ilp_r       = c(888, 999),
    precare     = 99,
    previs      = 99,
    cig_0       = 99,
    cig_1       = 99,
    cig_2       = 99,
    cig_3       = 99,
    m_ht_in     = 99,
    bmi         = 99.9,
    bmi_r       = 9,
    pwgt_r      = 999,
    dwgt_r      = 999,
    rf_pdiab    = 9,
    rf_gdiab    = 9,
    rf_phype    = 9,
    rf_ghype    = 9,
    rf_ehype    = 9,
    rf_ppterm   = 9,
    no_risks    = 9,
    dbwt        = 9999,
    bwgt_5      = 9,
    mm_aicu     = 9,
    ab_nicu     = 9
  ))

# Remove leftover outliers manually
bw_na <- bw_na %>%
  mutate(
    illb_r = replace(illb_r, which(illb_r > 300), NA),
    ilp_r  = replace(ilp_r,  which(ilp_r  > 300), NA),
    bmi    = replace(bmi,    which(bmi    > 70),  NA)
  )

# Summary stats after cleaning
stargazer(bw_na, type = "text",
          title = "Cleaned Data Summary Statistics (with NAs)")

#-----------------------------------------------
# Drop all NAs
#-----------------------------------------------
bw2 <- na.omit(bw_na)

# Summary stats after dropping NAs
stargazer(bw2, type = "text",
          title = "Edited Data Summary Statistics (NAs Removed)")



#-----------------------------------------------
# Continue with analysis using bw2
#-----------------------------------------------
# Additional key variables
bw2$avg_cigs <- (bw2$cig_1 + bw2$cig_2 + bw2$cig_3) / 3
bw2$log_bw   <- log(bw2$dbwt)

# Drop rows with missing in regression vars
bw2_clean <- na.omit(bw2[, c("log_bw", "avg_cigs", "rf_pdiab",
                             "rf_phype", "rf_ppterm", "sex", "bmi", "dbwt")])

# Run OLS regression
model <- lm(log_bw ~ avg_cigs + rf_pdiab + rf_phype +
              rf_ppterm + sex + bmi,
            data = bw2)

#-----------------------------------------------
# Table 2: Summary stats for regression variables
#-----------------------------------------------
bw2_clean_pretty <- bw2_clean %>%
  rename(
    "Avg cigarettes per day"     = avg_cigs,
    "Pre-pregnancy diabetes"     = rf_pdiab,
    "Pre-pregnancy hypertension" = rf_phype,
    "Previous preterm birth"     = rf_ppterm,
    "Infant's sex"               = sex,
    "Mother's BMI"               = bmi,
    "Birthweight"           = dbwt
  )

stargazer(bw2_clean_pretty, type = "text",
          title        = "Summary Statistics for Regression Variables",
          digits       = 3,
          summary.stat = c("mean", "sd", "min", "max", "median", "n"))

#-----------------------------------------------
# Table 3: Example regression results
#-----------------------------------------------

# OLS regression
model <- lm(log_bw ~ avg_cigs + rf_pdiab + rf_phype + rf_ppterm + sex + bmi,
            data = bw2_clean)

stargazer(model, type = "text",
          title = "OLS Regression Results: Birthweight and Smoking",
          dep.var.labels = "Log(Birthweight)",
          covariate.labels = c("Avg. Cigarettes/Day",
                               "Pre-pregnancy Diabetes",
                               "Pre-pregnancy Hypertension",
                               "Previous Preterm Birth",
                               "Infant Sex",
                               "Mother's BMI"),
          omit.stat = c("f", "ser"),
          no.space = TRUE)

#-----------------------------------------------
# Table 4: RESET test
#-----------------------------------------------
library(lmtest)
resettest(model)

#-----------------------------------------------
# Table 5: Joint hypothesis tests
#-----------------------------------------------
# install.packages("car")
library(car)

linearHypothesis(model, c("rf_pdiab = 0",
                          "rf_phype = 0",
                          "rf_ppterm = 0",
                          "sex = 0",
                          "bmi = 0"))

linearHypothesis(model, c("rf_pdiab = 0",
                          "rf_phype = 0"))

linearHypothesis(model, c("avg_cigs = 0", 
                          "sex = 0", "bmi=0"))

#-----------------------------------------------
# Table 6: Variance inflation factors
#-----------------------------------------------
vif(model)

#-----------------------------------------------
# Table 7: Breusch–Pagan test
#-----------------------------------------------
bp_test <- bptest(model)
print(bp_test)


# Required packages
library(dplyr)
library(broom)     # tidy
library(car)       # linearHypothesis
library(knitr)     # kable

base_formula <- log_bw ~ avg_cigs + rf_pdiab + rf_phype + rf_ppterm + sex + bmi
base_model <- lm(base_formula, data = bw2_clean)

# Terms to try: quadratics and interactions (one-at-a-time)
candidate_terms <- c(
  "I(avg_cigs^2)",          # quadratic in smoking
  "I(bmi^2)",               # quadratic in BMI
  "avg_cigs:bmi",           # interaction smoking * BMI
  "avg_cigs:sex",           # interaction smoking * sex
  "bmi:sex",                # interaction BMI * sex
  "I(avg_cigs^2):sex",      # quadratic smoking * sex (if you want)
  "I(bmi^2):sex"            # quadratic bmi * sex (if you want)
)

# A helper to safely make readable labels for output
nice_label <- function(term) {
  term %>% gsub("I\\(", "", .) %>% gsub("\\)", "", .) %>% gsub(":", " x ", .)
}

# Container for results
results <- data.frame(
  term = character(),
  term_type = character(),
  p_anova = numeric(),
  adj_r_squared = numeric(),
  delta_adj_r2 = numeric(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

# Fit each candidate model (base + term) and collect stats
base_adj_r2 <- summary(base_model)$adj.r.squared
base_aic <- AIC(base_model)

for(term in candidate_terms){
  # new formula
  new_formula <- as.formula(paste(". ~ . +", term))
  new_model <- update(base_model, new_formula)
  
  # Compare with ANOVA to test joint significance of added term(s)
  an <- tryCatch(anova(base_model, new_model), error = function(e) NULL)
  if(is.null(an)) {
    p_val <- NA
  } else {
    # anova returns a table; p-value in last row's "Pr(>F)"
    p_val <- an$`Pr(>F)`[nrow(an)]
  }
  
  # Adjusted R-squared and delta
  new_adj_r2 <- summary(new_model)$adj.r.squared
  delta_adj <- new_adj_r2 - base_adj_r2
  new_aic <- AIC(new_model)
  
  # Term type (quadratic vs interaction)
  type <- ifelse(grepl("\\^2|I\\(", term) & !grepl(":", term), "quadratic",
                 ifelse(grepl(":", term), "interaction", "other"))
  
  results <- rbind(results, data.frame(
    term        = nice_label(term),
    term_type   = type,
    p_anova     = p_val,
    adj_r_squared = round(new_adj_r2, 4),
    delta_adj_r2  = round(delta_adj, 4),
    AIC           = round(new_aic, 2),
    stringsAsFactors = FALSE
  ))
}

# Add the base model row for reference
base_row <- data.frame(
  term = "Base (no additions)",
  term_type = "base",
  p_anova = NA,
  adj_r_squared = round(base_adj_r2, 4),
  delta_adj_r2 = 0,
  AIC = round(base_aic, 2),
  stringsAsFactors = FALSE
)

summary_table <- bind_rows(base_row, results)

# Order by delta_adj_r2 (largest positive improvement first)
summary_table <- summary_table %>% arrange(desc(delta_adj_r2))

# Print nice table
kable(summary_table, caption = "One-at-a-time additions: significance and fit statistics")





best_terms <- c("I(avg_cigs^2)", "avg_cigs:bmi")  
best_formula <- as.formula(paste("log_bw ~ avg_cigs + rf_pdiab + rf_phype + rf_ppterm + sex + bmi +",
                                 paste(best_terms, collapse = " + ")))
best_model <- lm(best_formula, data = bw2_clean)

# Summary and RESET test
cat("\n\n--- Best combined model summary ---\n")
print(summary(best_model))

cat("\n\n--- RESET test on best_model ---\n")
library(lmtest)
print(resettest(best_model))

# Also give VIFs for the best model to check collinearity
cat("\n\n--- VIFs ---\n")
print(vif(best_model))

# ANOVA comparing base -> best combined model
cat("\n\n--- ANOVA base vs best combined model ---\n")
print(anova(base_model, best_model))

# If you prefer to save the summary table to CSV:
# write.csv(summary_table, "model_addition_summary.csv", row.names = FALSE)


#-----------------------------------------------
# Load required libraries
#-----------------------------------------------
library(dplyr)
library(lmtest)
library(car)
library(broom)
library(knitr)

#-----------------------------------------------
# Define 6 candidate models
#-----------------------------------------------
formulas <- list(
  # Model 1: Baseline linear model
  m1 = log_bw ~ avg_cigs + rf_pdiab + rf_phype + rf_ppterm + sex + bmi,
  
  # Model 2: Add quadratic terms for BMI and cigarettes
  m2 = log_bw ~ avg_cigs + I(avg_cigs^2) + rf_pdiab + rf_phype + rf_ppterm + sex + bmi + I(bmi^2),
  
  # Model 3: Add interaction between smoking and BMI (weight moderates smoking impact)
  m3 = log_bw ~ avg_cigs * bmi + rf_pdiab + rf_phype + rf_ppterm + sex + I(bmi^2),
  
  # Model 4: Allow sex to moderate BMI and smoking effects
  m4 = log_bw ~ avg_cigs * sex + bmi * sex + rf_pdiab + rf_phype + rf_ppterm + I(avg_cigs^2) + I(bmi^2),
  
  # Model 5: Full flexible model – combines key nonlinearities & interactions
  m5 = log_bw ~ avg_cigs * bmi + avg_cigs * sex + bmi * sex +
    rf_pdiab + rf_phype + rf_ppterm + I(avg_cigs^2) + I(bmi^2),

 # best expected model based on anova test
  m6 = log_bw ~ avg_cigs * bmi + rf_pdiab + rf_phype + rf_ppterm + bmi + I(bmi^2) + avg_cigs + I(avg_cigs^2) + sex +
   I(bmi^2):sex + avg_cigs:bmi, data = bw2_clean
)
#-----------------------------------------------
# Estimate all models
#-----------------------------------------------



models <- lapply(formulas, lm, data = bw2_clean)



stargazer(model, type = "text",
          title = "OLS Regression Results: Birthweight and Smoking",
          dep.var.labels = "Log(Birthweight)",
          covariate.labels = c("Avg. Cigarettes/Day",
                               "Pre-pregnancy Diabetes",
                               "Pre-pregnancy Hypertension",
                               "Previous Preterm Birth",
                               "Infant Sex",
                               "Mother's BMI"),
          omit.stat = c("f", "ser"),
          no.space = TRUE)

#-----------------------------------------------
# Create summary comparison table
#-----------------------------------------------
comparison <- data.frame(
  Model = names(models),
  Adj_R2 = sapply(models, function(x) summary(x)$adj.r.squared),
  AIC = sapply(models, AIC),
  RESET_pval = sapply(models, function(x) resettest(x)$p.value),
  BP_pval = sapply(models, function(x) bptest(x)$p.value)
)
print(comparison)

# Add a brief qualitative interpretation flag for RESET
comparison <- comparison %>%
  mutate(
    RESET_Signif = ifelse(RESET_pval < 0.05, "Fail (nonlinearities remain)", "Pass"),
    BP_Signif    = ifelse(BP_pval < 0.05, "Heteroskedasticity", "OK")
  )

#-----------------------------------------------
#Print comparison table
#-----------------------------------------------
kable(
  comparison %>%
    mutate(across(where(is.numeric), ~round(., 4))),
  caption = "Model Comparison: Adjusted R², AIC, RESET and BP Tests"
)

#-----------------------------------------------
# Detailed summaries & VIF checks
#-----------------------------------------------
for (i in seq_along(models)) {
  cat("\n\n==========================")
  cat(paste("\nModel", names(models)[i], "summary\n"))
  print(summary(models[[i]]))
  cat("\nRESET test:\n")
  print(resettest(models[[i]]))
  cat("\nBreusch-Pagan test:\n")
  print(bptest(models[[i]]))
  cat("\nVIFs:\n")
  print(vif(models[[i]]))
}

#-----------------------------------------------
# select models for regression results table
#-----------------------------------------------

model5 <- lm(log_bw ~ avg_cigs * bmi + avg_cigs * sex + bmi * sex +
               rf_pdiab + rf_phype + rf_ppterm + I(avg_cigs^2) + I(bmi^2), 
             data = bw2_clean)


model3 <- lm(log_bw ~ avg_cigs * bmi + rf_pdiab + rf_phype + rf_ppterm + sex + I(bmi^2), 
   data = bw2_clean)

model6 <- lm(log_bw ~ avg_cigs * bmi + rf_pdiab + rf_phype + rf_ppterm + bmi + I(bmi^2) + avg_cigs + I(avg_cigs^2) + sex +
                     I(bmi^2):sex + avg_cigs:bmi, data = bw2_clean)


stargazer(model, model3, model5, model6, type = "text",
          title = "OLS Regression Results: Birthweight and Smoking",
          dep.var.labels = "Log(Birthweight)",
          covariate.labels = c("Avg. Cigarettes/Day",
                               "Pre-pregnancy Diabetes",
                               "Pre-pregnancy Hypertension",
                               "Previous Preterm Birth",
                               "Infant Sex",
                               "Mother's BMI"),
          omit.stat = c("f", "ser"),
          no.space = TRUE)

ulmodel6 <- lm(dbwt ~ avg_cigs * bmi + rf_pdiab + rf_phype + rf_ppterm + bmi + I(bmi^2) + avg_cigs + I(avg_cigs^2) + sex +
                 I(bmi^2):sex + avg_cigs:bmi, data = bw2_clean)

ulmodel5 <- lm(dbwt ~ avg_cigs * bmi + avg_cigs * sex + bmi * sex +
               rf_pdiab + rf_phype + rf_ppterm + I(avg_cigs^2) + I(bmi^2), 
             data = bw2_clean)
ulmodel3 <- lm(dbwt ~ avg_cigs * bmi + rf_pdiab + rf_phype + rf_ppterm + sex + I(bmi^2), 
             data = bw2_clean)
ulmodel <- lm(dbwt ~ avg_cigs + rf_pdiab + rf_phype + rf_ppterm + sex
                       + bmi,
                       data = bw2_clean)


stargazer(ulmodel, ulmodel3, ulmodel5, ulmodel6, type = "text",
          title = "OLS Regression Results: Birthweight and Smoking",
          dep.var.labels = "Birthweight",
          covariate.labels = c("Avg. Cigarettes/Day",
                               "Pre-pregnancy Diabetes",
                               "Pre-pregnancy Hypertension",
                               "Previous Preterm Birth",
                               "Infant Sex",
                               "Mother's BMI"),
          omit.stat = c("f", "ser"),
          no.space = TRUE)









# -----------------------------
# Robust SEs + Stargazer tables
# -----------------------------
library(sandwich)   
library(lmtest) 

# Helper to extract robust SE vector for a model
rob_se <- function(fit, type = "HC1"){
  se_vec <- sqrt(diag(vcovHC(fit, type = type)))
  # Match order of coefficients returned by stargazer (intercept included)
  return(se_vec)
}

# Robust SEs for log models
rob_se_base_HC1 <- rob_se(model, type = "HC1")
rob_se_m3_HC1   <- rob_se(model3, type = "HC1")
rob_se_m5_HC1   <- rob_se(model5, type = "HC1")
rob_se_m6_HC1   <- rob_se(model6, type = "HC1")


# Robust SEs for level (dbwt) models
rob_se_ulbase_HC1 <- rob_se(ulmodel, type = "HC1")
rob_se_ul3_HC1    <- rob_se(ulmodel3, type = "HC1")
rob_se_ul5_HC1    <- rob_se(ulmodel5, type = "HC1")
rob_se_ul6_HC1    <- rob_se(ulmodel6, type = "HC1")


# -----------------------------
# Stargazer tables: LOG(Y) models
# -----------------------------
# 1) Conventional SEs (baseline)
stargazer(model, model3, model5, model6,
          type = "text",
          title = "OLS: Log(Birthweight) — Conventional SEs",
          dep.var.labels = "Log(Birthweight)",
          omit.stat = c("f", "ser"),
          no.space = TRUE)

# 2) Robust SEs (HC1)
stargazer(model, model3, model5, model6,
          se = list(rob_se_base_HC1, rob_se_m3_HC1, rob_se_m5_HC1, rob_se_m6_HC1),
          type = "text",
          title = "OLS: Log(Birthweight) — HC1 robust SEs",
          dep.var.labels = "Log(Birthweight)",
          omit.stat = c("f", "ser"),
          no.space = TRUE,
          notes = "Robust SEs: HC1")


# -----------------------------
# Stargazer tables: LEVEL(Y) models (dbwt)
# -----------------------------
stargazer(ulmodel, ulmodel3, ulmodel5, ulmodel6,
          type = "text",
          title = "OLS: Birthweight (grams) — Conventional SEs",
          dep.var.labels = "Birthweight (grams)",
          omit.stat = c("f", "ser"),
          no.space = TRUE)

stargazer(ulmodel, ulmodel3, ulmodel5, ulmodel6,
          se = list(rob_se_ulbase_HC1, rob_se_ul3_HC1, rob_se_ul5_HC1, rob_se_ul6_HC1),
          type = "text",
          title = "OLS: Birthweight (grams) — HC1 robust SEs",
          dep.var.labels = "Birthweight (grams)",
          omit.stat = c("f", "ser"),
          no.space = TRUE,
          notes = "Robust SEs: HC1")

# -----------------------------
# Stargazer tables: LOG(Y) models with HC1 robust SEs and 95% CIs
# -----------------------------
stargazer(model, model3, model5, model6,
          se = list(rob_se_base_HC1, rob_se_m3_HC1, rob_se_m5_HC1, rob_se_m6_HC1),
          type = "text",
          title = "OLS: Log(Birthweight) — HC1 Robust SEs with 95% Confidence Intervals",
          dep.var.labels = "Log(Birthweight)",
          omit.stat = c("f", "ser"),
          no.space = TRUE,
          ci = TRUE, ci.level = 0.95,
          notes = "Robust SEs: HC1. 95% confidence intervals shown."
)

# -----------------------------
# Stargazer tables: LEVEL(Y) models (dbwt) with HC1 robust SEs and 95% CIs
# -----------------------------
stargazer(ulmodel, ulmodel3, ulmodel5, ulmodel6,
          se = list(rob_se_ulbase_HC1, rob_se_ul3_HC1, rob_se_ul5_HC1, rob_se_ul6_HC1),
          type = "text",
          title = "OLS: Birthweight (grams) — HC1 Robust SEs with 95% Confidence Intervals",
          dep.var.labels = "Birthweight (grams)",
          omit.stat = c("f", "ser"),
          no.space = TRUE,
          ci = TRUE, ci.level = 0.95,
          notes = "Robust SEs: HC1. 95% confidence intervals shown."
)













