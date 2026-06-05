# =========================================================
# Simulations - Eurobond DRC project
# =========================================================
# This script implements:
#   1) PPE significance tests using resampling of the PPE(h) curve
#   2) Country-level PPE comparison vs. the Top-10 most similar countries
#   3) Learning-curve simulations for future RDC spreads and PPE paths
#
# Expected objects in the workspace:
#   - clean_df
#   - estimation_store
#
# Expected helper functions already available (from the estimation script):
#   - prepare_subset()
#   - complete_case_filter()
#   - varying_vars()
#   - mahalanobis_distance()
#   - gower_distance()
#   - kernel_weights()
#   - build_formula()
# =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(readxl)
library(tibble)
library(tidyr)

# ---------------------------------------------------------
# 0) Input checks and folders
# ---------------------------------------------------------
if (!exists("clean_df")) {
  stop("Object `clean_df` not found. Run the preprocessing and estimation scripts first.")
}
if (!exists("estimation_store")) {
  stop("Object `estimation_store` not found. Run the estimation script first.")
}
if (!exists("prepare_subset") || !exists("complete_case_filter") || !exists("varying_vars") ||
    !exists("mahalanobis_distance") || !exists("gower_distance") || !exists("kernel_weights") ||
    !exists("build_formula")) {
  stop("Some helper functions from the estimation script are missing. Source the estimation script first.")
}

dir.create("simulation_outputs", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("simulation_outputs", "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("simulation_outputs", "tables"), showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# ---------------------------------------------------------
# 1) Model specification (must match estimation script)
# ---------------------------------------------------------
base_rhs_vars <- c(
  "debt", "reserves", "growth", "inflation",
  "current_account", "fiscal_balance",
  "volatility", "sp_note", "law_est", "hipc"
)

outcome_vars <- c("yas_spread", "asw_spread")
distance_methods <- c("mahalanobis", "gower")
kernels <- c("gaussian", "epanechnikov")
installment_specs <- list(
  A = c("A", "Unique"),
  B = c("B", "Unique")
)

anchor_year <- 2025
future_years <- 2025:2031
lambda_scenarios <- c(0, 10, 25)

# ---------------------------------------------------------
# 2) Utility functions
# ---------------------------------------------------------

safe_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  parse_number(as.character(x))
}

country_weighted_mean <- function(df, value_var, weight_var = "amount_usd") {
  if (!value_var %in% names(df)) return(NA_real_)
  x <- safe_numeric(df[[value_var]])
  x[!is.finite(x)] <- NA_real_
  
  if (weight_var %in% names(df)) {
    w <- safe_numeric(df[[weight_var]])
    w[!is.finite(w)] <- NA_real_
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (sum(ok) > 0 && sum(w[ok]) > 0) {
      return(weighted.mean(x[ok], w[ok], na.rm = TRUE))
    }
  }
  
  mean(x, na.rm = TRUE)
}

country_total_amount <- function(df, weight_var = "amount_usd") {
  if (!weight_var %in% names(df)) return(NA_real_)
  w <- safe_numeric(df[[weight_var]])
  sum(w[is.finite(w) & w > 0], na.rm = TRUE)
}

select_target_rows <- function(df, target_country, installment_keep) {
  df %>%
    filter(country == target_country, installment %in% installment_keep)
}

get_target_curve_estimate <- function(data_all,
                                      target_country,
                                      installment_keep,
                                      outcome,
                                      rhs_vars,
                                      distance_vars,
                                      distance_method = c("mahalanobis", "gower"),
                                      kernel = c("gaussian", "epanechnikov"),
                                      h_grid,
                                      n_min_effective = 2) {
  distance_method <- match.arg(distance_method)
  kernel <- match.arg(kernel)
  
  ref_df <- prepare_subset(data_all, installment_keep) %>%
    filter(country != target_country)
  target_df <- select_target_rows(data_all, target_country, installment_keep)
  
  if (nrow(target_df) == 0) {
    stop(paste0("No target observation found for country: ", target_country))
  }
  
  model_formula <- build_formula(outcome, rhs_vars)
  model_vars <- setdiff(all.vars(model_formula), outcome)
  
  ref_keep_vars <- unique(c(outcome, model_vars, distance_vars))
  tar_keep_vars <- unique(c(outcome, model_vars, distance_vars, "amount_usd"))
  
  ref_df <- complete_case_filter(ref_df, ref_keep_vars)
  target_df <- complete_case_filter(target_df, tar_keep_vars)
  
  if (nrow(ref_df) < length(model_vars) + 2) {
    stop("Too few complete reference observations after filtering.")
  }
  
  distance_vec <- switch(
    distance_method,
    mahalanobis = mahalanobis_distance(ref_df, target_df, distance_vars),
    gower = gower_distance(ref_df, target_df, distance_vars)
  )
  
  distance_vec <- as.numeric(distance_vec)
  keep_idx <- is.finite(distance_vec)
  ref_df <- ref_df[keep_idx, , drop = FALSE]
  distance_vec <- distance_vec[keep_idx]
  
  if (length(distance_vec) == 0 || nrow(ref_df) == 0) {
    stop("No usable reference observations after distance filtering.")
  }
  
  rmse_tbl <- vector("list", length(h_grid))
  weight_tbl <- vector("list", length(h_grid))
  
  obs_spread_country <- country_weighted_mean(target_df, outcome, weight_var = "amount_usd")
  best_rmse <- Inf
  best_h <- NA_real_
  best_model <- NULL
  best_pred_spread <- NA_real_
  best_ppe <- NA_real_
  
  for (j in seq_along(h_grid)) {
    h <- h_grid[j]
    if (!is.finite(h) || h <= 0) {
      rmse_tbl[[j]] <- tibble(h = h, rmse = Inf, n_effective = 0, r2 = NA_real_, ppe = NA_real_, observed_spread = obs_spread_country, predicted_spread = NA_real_)
      next
    }
    
    w <- kernel_weights(distance_vec, h, kernel = kernel)
    weight_tbl[[j]] <- tibble(
      h = h,
      country = ref_df$country,
      distance = distance_vec,
      weight = w
    )
    
    reg_df <- ref_df %>%
      mutate(.w = w) %>%
      filter(is.finite(.w), .w > 0) %>%
      filter(if_all(all_of(unique(c(outcome, model_vars))), ~ !is.na(.x)))
    
    if (nrow(reg_df) < n_min_effective) {
      rmse_tbl[[j]] <- tibble(h = h, rmse = Inf, n_effective = nrow(reg_df), r2 = NA_real_, ppe = NA_real_, observed_spread = obs_spread_country, predicted_spread = NA_real_)
      next
    }
    
    fit <- tryCatch(
      lm(model_formula, data = reg_df, weights = .w),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      rmse_tbl[[j]] <- tibble(h = h, rmse = Inf, n_effective = nrow(reg_df), r2 = NA_real_, ppe = NA_real_, observed_spread = obs_spread_country, predicted_spread = NA_real_)
      next
    }
    
    pred_target_rows <- tryCatch(
      as.numeric(predict(fit, newdata = target_df)),
      error = function(e) rep(NA_real_, nrow(target_df))
    )
    
    pred_spread_country <- country_weighted_mean(
      tibble(pred = pred_target_rows, amount_usd = target_df$amount_usd),
      value_var = "pred",
      weight_var = "amount_usd"
    )
    
    resid <- residuals(fit)
    rmse <- sqrt(sum(reg_df$.w * resid^2, na.rm = TRUE) / sum(reg_df$.w, na.rm = TRUE))
    r2 <- summary(fit)$r.squared
    ppe_h <- obs_spread_country - pred_spread_country
    
    rmse_tbl[[j]] <- tibble(
      h = h,
      rmse = rmse,
      n_effective = nrow(reg_df),
      r2 = r2,
      ppe = ppe_h,
      observed_spread = obs_spread_country,
      predicted_spread = pred_spread_country
    )
    
    if (rmse < best_rmse) {
      best_rmse <- rmse
      best_h <- h
      best_model <- fit
      best_pred_spread <- pred_spread_country
      best_ppe <- ppe_h
    }
  }
  
  rmse_tbl <- bind_rows(rmse_tbl)
  weight_tbl <- bind_rows(weight_tbl)
  
  list(
    ref_df = ref_df,
    target_df = target_df,
    model_formula = model_formula,
    h_grid = h_grid,
    distance_method = distance_method,
    kernel = kernel,
    rmse_table = rmse_tbl,
    weight_table = weight_tbl,
    best_h = best_h,
    best_rmse = best_rmse,
    best_model = best_model,
    best_obs_spread = obs_spread_country,
    best_pred_spread = best_pred_spread,
    best_ppe = best_ppe
  )
}

bootstrap_curve_median <- function(ppe_values, B = 2000, conf = 0.95) {
  ppe_values <- ppe_values[is.finite(ppe_values)]
  if (length(ppe_values) == 0) {
    return(tibble(
      median = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      share_positive = NA_real_,
      ppe_min = NA_real_,
      ppe_max = NA_real_
    ))
  }
  
  n <- length(ppe_values)
  boot_meds <- replicate(B, median(sample(ppe_values, size = n, replace = TRUE)))
  
  tibble(
    median = median(ppe_values),
    ci_low = as.numeric(quantile(boot_meds, probs = (1 - conf) / 2, na.rm = TRUE)),
    ci_high = as.numeric(quantile(boot_meds, probs = 1 - (1 - conf) / 2, na.rm = TRUE)),
    share_positive = mean(ppe_values > 0, na.rm = TRUE),
    ppe_min = min(ppe_values, na.rm = TRUE),
    ppe_max = max(ppe_values, na.rm = TRUE)
  )
}

bootstrap_curve_band <- function(curve_df, B = 2000, conf = 0.95) {
  curve_df <- curve_df %>%
    arrange(h) %>%
    filter(is.finite(ppe))
  
  ppe_vals <- curve_df$ppe
  n <- length(ppe_vals)
  
  if (n == 0) {
    return(curve_df %>%
             mutate(lower = NA_real_, upper = NA_real_))
  }
  
  boot_mat <- replicate(B, sample(ppe_vals, size = n, replace = TRUE))
  
  curve_df %>%
    mutate(
      lower = apply(boot_mat, 1, quantile, probs = (1 - conf) / 2, na.rm = TRUE),
      upper = apply(boot_mat, 1, quantile, probs = 1 - (1 - conf) / 2, na.rm = TRUE)
    )
}

plot_ppe_curve <- function(curve_df, boot_stats, title_prefix) {
  band_df <- bootstrap_curve_band(curve_df)
  
  ggplot(band_df, aes(x = h, y = ppe)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey75", alpha = 0.55) +
    #geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_hline(yintercept = boot_stats$median, linetype = "dashed", color = "gray20") +
    geom_line(linewidth = 1.2, color = "black") +
    geom_point(size = 1.6, color = "black") +
    labs(
      title = paste0("MEP curve | ", title_prefix),
      subtitle = paste0(
        "Median = ", round(boot_stats$median, 2),
        "; 95% CI [", round(boot_stats$ci_low, 2), ", ", round(boot_stats$ci_high, 2), "]",
        "; Share positive = ", round(100 * boot_stats$share_positive, 1), "%"
      ),
      x = "Bandwidth h",
      y = "MEP(h)"
    ) +
    scale_x_continuous(limits = c(2, 13.5), n.breaks = 7) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 24, face = "bold"),
          plot.subtitle = element_text(size = 20),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
}

get_top_similar_countries <- function(estimation_object, top_n = 10) {
  best_h <- estimation_object$best_h
  weight_tbl <- estimation_object$weight_table %>%
    filter(dplyr::near(h, best_h)) %>%
    group_by(country) %>%
    summarise(
      weight = max(weight, na.rm = TRUE),
      distance = min(distance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(weight))
  
  weight_tbl %>% slice_head(n = top_n)
}

normalize_future_names <- function(df) {
  names(df) <- tolower(trimws(names(df)))
  names(df) <- gsub("\\s+", "_", names(df))
  if ("10d_volatility" %in% names(df) && !"volatility" %in% names(df)) {
    names(df)[names(df) == "10d_volatility"] <- "volatility"
  }
  df
}

get_anchor_spread_2025 <- function(future_raw, outcome) {
  future_raw <- normalize_future_names(future_raw)
  
  year_col <- intersect(c("year", "forecast_year", "date"), names(future_raw))[1]
  if (is.na(year_col)) {
    stop("No year column found in the future DRC sheet.")
  }
  
  df_2025 <- future_raw %>%
    mutate(.year = safe_numeric(.data[[year_col]])) %>%
    filter(.year == 2025)
  
  if (nrow(df_2025) == 0) {
    stop("No 2025 row found in the future DRC sheet.")
  }
  
  if (outcome %in% names(df_2025)) {
    vals <- safe_numeric(df_2025[[outcome]])
    vals[!is.finite(vals)] <- NA_real_
    
    if ("amount_usd" %in% names(df_2025)) {
      w <- safe_numeric(df_2025$amount_usd)
      ok <- is.finite(vals) & is.finite(w) & w > 0
      if (sum(ok) > 0) {
        return(weighted.mean(vals[ok], w[ok], na.rm = TRUE))
      }
    }
    
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0) return(vals[1])
  }
  
  stop(paste0("Outcome column `", outcome, "` not found or empty in 2025 future data."))
}

build_future_projection_rows <- function(future_raw, years_vec) {
  future_raw <- normalize_future_names(future_raw)
  
  if ("country" %in% names(future_raw)) {
    future_raw <- future_raw %>% filter(is.na(country) | country == "RDC")
  }
  
  year_col <- intersect(c("year", "forecast_year", "date"), names(future_raw))[1]
  if (is.na(year_col)) {
    stop("No year column found in the future DRC sheet.")
  }
  
  future_raw[[year_col]] <- safe_numeric(future_raw[[year_col]])
  future_raw <- future_raw %>% arrange(.data[[year_col]])
  
  rows <- lapply(years_vec, function(y) {
    row <- future_raw %>% filter(.data[[year_col]] == y)
    if (nrow(row) == 0) {
      row <- future_raw %>% filter(.data[[year_col]] <= y)
      if (nrow(row) == 0) {
        row <- future_raw %>% slice_tail(n = 1)
      } else {
        row <- row %>% slice_tail(n = 1)
      }
    }
    row <- row %>% slice(1)
    row$forecast_year <- y
    row
  })
  
  bind_rows(rows)
}

get_anchor_values_2025 <- function(data_all, installment_keep, vars_to_fill) {
  anchor_df <- data_all %>%
    filter(country == "RDC", installment %in% installment_keep)
  
  if (nrow(anchor_df) == 0) {
    stop("No RDC anchor row found in clean_df for the selected installment group.")
  }
  
  anchor_values <- sapply(vars_to_fill, function(v) {
    if (!v %in% names(anchor_df)) return(0)
    
    x <- safe_numeric(anchor_df[[v]])
    x[!is.finite(x)] <- NA_real_
    x <- x[is.finite(x)]
    
    if (length(x) == 0) return(0)
    mean(x, na.rm = TRUE)
  })
  
  names(anchor_values) <- vars_to_fill
  anchor_values
}

fill_missing_model_vars_naive <- function(df, model_vars, anchor_values) {
  out <- df
  
  for (v in model_vars) {
    if (!v %in% names(out)) out[[v]] <- anchor_values[[v]]
    out[[v]] <- safe_numeric(out[[v]])
    
    na_idx <- !is.finite(out[[v]]) | is.na(out[[v]])
    if (any(na_idx)) {
      out[[v]][na_idx] <- anchor_values[[v]]
    }
  }
  
  out
}

simulate_learning_curve <- function(best_model,
                                    future_df,
                                    anchor_spread_2025,
                                    anchor_values_2025,
                                    anchor_year = 2025,
                                    years_vec,
                                    lambda_vec = c(0, 10, 25)) {
  model_formula <- formula(best_model)
  model_vars <- attr(terms(model_formula), "term.labels")
  future_df <- fill_missing_model_vars_naive(future_df, model_vars, anchor_values_2025)
  
  pred_spread <- tryCatch(
    as.numeric(predict(best_model, newdata = future_df)),
    error = function(e) rep(NA_real_, nrow(future_df))
  )
  
  future_tbl <- tibble(
    forecast_year = future_df$forecast_year,
    predicted_spread = pred_spread
  )
  
  scenario_tbl <- expand.grid(
    forecast_year = years_vec,
    lambda = lambda_vec,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    arrange(lambda, forecast_year) %>%
    group_by(lambda) %>%
    mutate(
      observed_spread = pmax(
        anchor_spread_2025 - lambda * pmax(forecast_year - anchor_year, 0),
        0
      )
    ) %>%
    ungroup() %>%
    left_join(future_tbl, by = "forecast_year") %>%
    mutate(
      ppe = observed_spread - predicted_spread
    )
  
  scenario_tbl
}

plot_learning_curve <- function(learning_tbl, title_prefix) {
  label_df <- learning_tbl %>%
    group_by(lambda) %>%
    slice_max(order_by = forecast_year, n = 1, with_ties = FALSE) %>%
    mutate(
      scenario_label = case_when(
        lambda == 0  ~ "Scénario 0",
        lambda == 10 ~ "Scénario 1",
        lambda == 25 ~ "Scénario 2",
        TRUE ~ paste0("Scénario : lambda=", lambda)
      )
    ) %>%
    ungroup()
  
  ggplot(learning_tbl, aes(x = forecast_year, y = ppe, group = lambda)) +
    #geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_line(color = "black", linewidth = 1.2) +
    geom_point(color = "black", size = 1.8) +
    geom_text(
      data = label_df,
      aes(label = scenario_label),
      fontface = "italic",
      hjust = 0,
      nudge_x = 0.15,
      size = 8,
      color = "black",
      show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = sort(unique(learning_tbl$forecast_year)),
      expand = expansion(mult = c(0.02, 0.22))
    ) +
    scale_y_continuous(limits = c(-50,150), n.breaks = 5,
                       expand = expansion(mult = c(0.08, 0.10))) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("Learning curve | ", title_prefix),
      x = " ",
      y = "MEP"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 35, 5.5, 5.5),
      plot.title = element_text(size = 24, face = "bold"),
      axis.title.x = element_text(size = 22),
      axis.title.y = element_text(size = 22),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20)
    )
}

# ---------------------------------------------------------
# 3) Load future RDC projections
# ---------------------------------------------------------
future_file_path <- "eurobond_dataset.xlsx"
if (!file.exists(future_file_path)) {
  stop("Future data file not found: eurobond_dataset.xlsx")
}
future_raw <- read_excel(future_file_path, sheet = "drc data")

# ---------------------------------------------------------
# 4) Master loops over all estimation combinations
# ---------------------------------------------------------

ppe_significance_results <- list()
country_comparison_tables <- list()
learning_curve_tables <- list()
curve_plots <- list()
learning_plots <- list()

for (inst_name in names(installment_specs)) {
  for (outcome in outcome_vars) {
    for (dist_m in distance_methods) {
      for (kernel in kernels) {
        
        combo_key <- paste(inst_name, outcome, dist_m, kernel, sep = "|")
        message("Processing combination: ", combo_key)
        
        combo_est <- estimation_store[[inst_name]][[outcome]][[dist_m]][[kernel]]
        h_grid_use <- combo_est$h_grid
        
        # -----------------------------------------------------------------
        # 4.1 RDC-specific PPE curve and bootstrap significance test
        # -----------------------------------------------------------------
        rdc_curve <- get_target_curve_estimate(
          data_all = clean_df,
          target_country = "RDC",
          installment_keep = installment_specs[[inst_name]],
          outcome = outcome,
          rhs_vars = base_rhs_vars,
          distance_vars = setdiff(attr(terms(combo_est$model_formula), "term.labels"), NULL),
          distance_method = dist_m,
          kernel = kernel,
          h_grid = h_grid_use
        )
        
        boot_stats <- bootstrap_curve_median(rdc_curve$rmse_table$ppe, B = 2000, conf = 0.95)
        significant_positive <- is.finite(boot_stats$ci_low) && boot_stats$ci_low > 0
        
        ppe_significance_results[[combo_key]] <- tibble(
          installment_group = inst_name,
          outcome = outcome,
          distance_method = dist_m,
          kernel = kernel,
          best_h = rdc_curve$best_h,
          best_rmse = rdc_curve$best_rmse,
          median_ppe = boot_stats$median,
          ci_low = boot_stats$ci_low,
          ci_high = boot_stats$ci_high,
          share_positive = boot_stats$share_positive,
          ppe_min = boot_stats$ppe_min,
          ppe_max = boot_stats$ppe_max,
          significant_positive = significant_positive
        )
        
        curve_plot <- plot_ppe_curve(
          curve_df = rdc_curve$rmse_table %>% filter(is.finite(ppe)),
          boot_stats = boot_stats,
          title_prefix = combo_key
        )
        curve_plots[[combo_key]] <- curve_plot
        ggsave(
          filename = file.path("simulation_outputs", "plots", paste0("ppe_curve_", gsub("[|]", "_", combo_key), ".pdf")),
          plot = curve_plot,
          width = 16,
          height = 8,
          units = "in"
        )
        
        # -----------------------------------------------------------------
        # 4.2 Country-level PPE comparison vs Top-10 similar countries
        # -----------------------------------------------------------------
        top10_tbl <- get_top_similar_countries(combo_est, top_n = 10)
        compare_countries <- c("RDC", top10_tbl$country)
        comp_rows <- vector("list", length(compare_countries))
        
        for (k in seq_along(compare_countries)) {
          tgt_country <- compare_countries[k]
          sim_rank <- if (tgt_country == "RDC") 0L else which(top10_tbl$country == tgt_country)
          
          tgt_curve <- tryCatch(
            get_target_curve_estimate(
              data_all = clean_df,
              target_country = tgt_country,
              installment_keep = installment_specs[[inst_name]],
              outcome = outcome,
              rhs_vars = base_rhs_vars,
              distance_vars = distance_vars,
              distance_method = dist_m,
              kernel = kernel,
              h_grid = h_grid_use
            ),
            error = function(e) NULL
          )
          
          if (is.null(tgt_curve)) {
            comp_rows[[k]] <- tibble(
              installment_group = inst_name,
              outcome = outcome,
              distance_method = dist_m,
              kernel = kernel,
              target_country = tgt_country,
              similarity_rank = sim_rank,
              best_h = NA_real_,
              best_rmse = NA_real_,
              observed_spread = NA_real_,
              predicted_spread = NA_real_,
              ppe = NA_real_,
              r_squared = NA_real_
            )
            next
          }
          
          comp_rows[[k]] <- tibble(
            installment_group = inst_name,
            outcome = outcome,
            distance_method = dist_m,
            kernel = kernel,
            target_country = tgt_country,
            similarity_rank = sim_rank,
            best_h = tgt_curve$best_h,
            best_rmse = tgt_curve$best_rmse,
            observed_spread = tgt_curve$best_obs_spread,
            predicted_spread = tgt_curve$best_pred_spread,
            ppe = tgt_curve$best_ppe,
            r_squared = if (!is.null(tgt_curve$best_model)) summary(tgt_curve$best_model)$r.squared else NA_real_
          )
        }
        
        comp_tbl <- bind_rows(comp_rows) %>%
          arrange(similarity_rank, desc(target_country == "RDC"), target_country)
        
        country_comparison_tables[[combo_key]] <- comp_tbl
        write.csv(
          comp_tbl,
          file = file.path("simulation_outputs", "tables", paste0("ppe_comparison_", gsub("[|]", "_", combo_key), ".csv")),
          row.names = FALSE
        )
        
        # -----------------------------------------------------------------
        # 4.3 Learning-curve simulation for RDC only
        # -----------------------------------------------------------------
        anchor_spread_2025 <- get_anchor_spread_2025(future_raw, outcome)
        future_rows <- build_future_projection_rows(future_raw, years_vec = future_years)
        model_vars_future <- attr(terms(combo_est$model_formula), "term.labels")
        vars_to_fill <- setdiff(model_vars_future, names(future_rows))
        anchor_values_2025 <- get_anchor_values_2025(clean_df, installment_specs[[inst_name]], vars_to_fill)
        future_rows <- fill_missing_model_vars_naive(future_rows, vars_to_fill, anchor_values_2025)
        
        learning_tbl <- simulate_learning_curve(
          best_model = rdc_curve$best_model,
          future_df = future_rows,
          anchor_spread_2025 = anchor_spread_2025,
          anchor_values_2025 = anchor_values_2025,
          anchor_year = anchor_year,
          years_vec = future_years,
          lambda_vec = lambda_scenarios
        )
        
        learning_tbl <- learning_tbl %>%
          mutate(
            installment_group = inst_name,
            outcome = outcome,
            distance_method = dist_m,
            kernel = kernel,
            combo_key = combo_key
          )
        
        learning_curve_tables[[combo_key]] <- learning_tbl
        learning_plot <- plot_learning_curve(learning_tbl, title_prefix = combo_key)
        learning_plots[[combo_key]] <- learning_plot
        ggsave(
          filename = file.path("simulation_outputs", "plots", paste0("learning_curve_", gsub("[|]", "_", combo_key), ".pdf")),
          plot = learning_plot,
          width = 16,
          height = 8,
          units = "in"
        )
      }
    }
  }
}

# ---------------------------------------------------------
# 5) Save master tables
# ---------------------------------------------------------
ppe_significance_tbl <- bind_rows(ppe_significance_results)
write.csv(
  ppe_significance_tbl,
  file = file.path("simulation_outputs", "tables", "ppe_significance_summary.csv"),
  row.names = FALSE
)

country_comparison_tbl <- bind_rows(country_comparison_tables)
write.csv(
  country_comparison_tbl,
  file = file.path("simulation_outputs", "tables", "ppe_country_comparison_all_combos.csv"),
  row.names = FALSE
)

learning_curve_tbl <- bind_rows(learning_curve_tables)
write.csv(
  learning_curve_tbl,
  file = file.path("simulation_outputs", "tables", "learning_curve_all_combos.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# 6) Console summary
# ---------------------------------------------------------
print(ppe_significance_tbl)
print(head(country_comparison_tbl, 20))
print(head(learning_curve_tbl, 20))

# Main objects created:
#   - ppe_significance_tbl
#   - country_comparison_tbl
#   - learning_curve_tbl
#   - curve_plots
#   - learning_plots
#   - ppe_significance_results
#   - country_comparison_tables
#   - learning_curve_tables
#
# =========================================================
# End of script
# =========================================================
