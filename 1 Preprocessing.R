# =========================================================
# Preprocessing - Eurobond DRC project
# =========================================================

# Tasks:
# 1) Import Excel data
# 2) Remove observations with missing yas_spread
# 3) Remove observations with simultaneous missing values on
#    debt, reserves, growth, inflation, current_account, fiscal_balance
# 4) Create regional dummies: africa, america, asia
# 5) Coerce variable types to their intended formats
# =========================================================

library(readxl)
library(dplyr)
library(stringr)
library(lubridate)
library(readr)

# ---------------------------------------------------------
# 1) Import data
# ---------------------------------------------------------

file_path <- "eurobond_dataset.xlsx"
sheet_name <- "cross sectional data"

raw_df <- read_excel(file_path, sheet = sheet_name)

# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------
# Convert numeric-like values safely, even if imported as text with commas
as_num <- function(x) {
  parse_number(
    as.character(x),
    locale = locale(decimal_mark = ".", grouping_mark = ",")
  )
}

# Parse dates robustly from common formats used in the dataset
as_dt <- function(x) {
  suppressWarnings(
    as.Date(parse_date_time(as.character(x), orders = c("d-b-Y", "d-b-y", "d/m/Y", "d-m-Y", "Y-m-d")))
  )
}

# ---------------------------------------------------------
# 2) Missing-value treatment
# ---------------------------------------------------------
# Rule A: drop rows with missing yas_spread
removed_yas <- raw_df %>%
  filter(is.na(yas_spread)) %>%
  mutate(drop_reason = "yas_spread missing")

step1_df <- raw_df %>%
  filter(!is.na(yas_spread))

# Rule B: drop rows where all core macro variables are missing simultaneously
macro_vars <- c("debt", "reserves", "growth", "inflation", "current_account", "fiscal_balance")

removed_macro <- step1_df %>%
  filter(if_all(all_of(macro_vars), is.na)) %>%
  mutate(drop_reason = "all macro variables missing")

clean_df <- step1_df %>%
  filter(!if_all(all_of(macro_vars), is.na))

# Combine removed rows and keep the reason(s)
removed_df <- bind_rows(removed_yas, removed_macro) %>%
  arrange(row_number()) %>%
  dplyr::select(country, region, installment, issue_date, maturity_date, yas_spread, all_of(macro_vars), drop_reason)

# NAs summary
na_summary <- clean_df %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "nb_NA") %>%
  arrange(desc(nb_NA)) %>%
  mutate(pct_NA = round(100 * nb_NA / nrow(clean_df), 2))
View(na_summary)

# ---------------------------------------------------------
# 3) Create regional dummies
# ---------------------------------------------------------
# We use the region field to build three dummies:
#   africa  = 1 if Africa
#   america = 1 if Latin America / America
#   asia    = 1 if Asia / Middle East

clean_df <- clean_df %>%
  mutate(
    africa  = if_else(region == "Afrique", 1, 0),
    america = if_else(region == "Amérique latine", 1, 0),
    asia    = if_else(region == "Asie & MO", 1, 0)
  )

# ---------------------------------------------------------
# 4) Correct variable formats
# ---------------------------------------------------------

# String variable
clean_df <- clean_df %>%
  mutate(country = as.character(country))

# Categorical variables (factors)
categorical_vars <- c("region", "installment", "currency", "maturity_type",
                      "bbg_composite", "sp_rating", "moodys_rating", "fitch_rating")

clean_df <- clean_df %>%
  mutate(across(all_of(categorical_vars), as.factor))

# Date variables
clean_df <- clean_df %>%
  mutate(
    issue_date = as.Date(issue_date, format = "%d-%b-%Y"),
    maturity_date = as.Date(maturity_date, format = "%d-%b-%Y")
  )

# Numerical variables
numeric_vars <- c("amount", "amount_usd", "maturity", "coupon_rate", 
                  "yas_spread", "asw_spread", "ytm", "debt", "reserves", 
                  "growth", "inflation", "current_account", "fiscal_balance",
                  "corruption_est", "corruption_scr", "effectiveness_est", 
                  "effectiveness_scr", "stability_est", "stability_scr", 
                  "regulatory_est", "regulatory_scr", "law_est", "law_scr", 
                  "voice_est", "voice_scr", "sp_note", "moodys_note",
                  "fitch_note", "volatility")

clean_df <- clean_df %>%
  mutate(across(all_of(numeric_vars), as.numeric))

# Integer variables
int_vars <- c("africa", "america", "asia", "hipc")

clean_df <- clean_df %>%
  mutate(across(all_of(int_vars), as.integer))

# Sanity check on types
str(clean_df)

# ---------------------------------------------------------
# 5) Optional quick checks
# ---------------------------------------------------------
# View removal summary
print(table(removed_df$drop_reason, useNA = "ifany"))
print(unique(removed_df$country))

# ---------------------------------------------------------
# 6) Save the cleaned dataset
# ---------------------------------------------------------

write.csv(clean_df, "cleaned_dataset.csv")


# The objects to use later are:
#   - clean_df   : cleaned dataset
#   - removed_df : observations removed and their reasons

# =========================================================
# End of script
# =========================================================