#-----------------------------------------------
# Load libraries
#-----------------------------------------------
library(tidyverse)
library(naniar)
library(stargazer)
library(lmtest)

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
                             "rf_phype", "rf_ppterm", "sex", "bmi")])

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
    "Log(Birthweight)"           = log_bw
  )

stargazer(bw2_clean_pretty, type = "text",
          title        = "Summary Statistics for Regression Variables",
          digits       = 3,
          summary.stat = c("mean", "sd", "min", "max", "median", "n"))

#-----------------------------------------------
# Table 3: Example summary statistics with kableExtra
#-----------------------------------------------
# install.packages("knitr")
# install.packages("kableExtra")

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

#-----------------------------------------------
# Table 6: Variance inflation factors
#-----------------------------------------------
vif(model)

#-----------------------------------------------
# Table 7: Breusch–Pagan test
#-----------------------------------------------
bp_test <- bptest(model)
print(bp_test)

