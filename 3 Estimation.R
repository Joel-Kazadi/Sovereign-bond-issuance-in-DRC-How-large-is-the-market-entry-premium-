# =========================================================
# Estimation - Eurobond DRC project
# =========================================================
# This script implements four stages:
#   1) Similarity measures: Mahalanobis and Gower distances
#   2) Kernel weights: Gaussian and Epanechnikov
#   3) Weighted least squares over a grid of bandwidths (h)
#   4) Similarity barplots at the optimal h (RMSE-minimizing)
#
# Notes:
# - The script expects `clean_df` to be available in the workspace.
# - The two RDC tranches (A and B) are handled separately via installment subsets:
#       A-subset: installment in c("A", "Unique")
#       B-subset: installment in c("B", "Unique")
# =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(cluster)
library(MASS)
library(readr)
library(sandwich)
library(lmtest)
library(tibble)
library(tidyr)

# ---------------------------------------------------------
# 0) Input checks, output folder and outliers removing
# ---------------------------------------------------------
if (!exists("clean_df")) {
  stop("Object `clean_df` not found. Run the preprocessing script first.")
}

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

clean_df <- clean_df %>%
  filter(!country %in% c("Canada", "Grenade", "Kazakhstan"))

# ---------------------------------------------------------
# 1) Model specification
# ---------------------------------------------------------
# Baseline controls (RHS).
# sp_note and effectiveness_est are now included directly in the baseline.
base_rhs_vars <- c(
  "debt", "reserves", "growth", "inflation",
  "current_account", "fiscal_balance",
  "volatility", "sp_note", "law_est", "hipc"
)

# Dependent variables.
outcome_vars <- c("yas_spread", "asw_spread")

# Variables used for similarity calculations.
distance_vars <- c(
  "debt", "reserves", "growth", "inflation",
  "current_account", "fiscal_balance",
  "law_est", "sp_note", "volatility"
)

# Installment subsets.
installment_specs <- list(
  A = c("A", "Unique"),
  B = c("B", "Unique")
)

# Distances and kernels to be compared.
distance_methods <- c("mahalanobis", "gower")
kernels <- c("gaussian", "epanechnikov")

# ---------------------------------------------------------
# 2) Helper functions
# ---------------------------------------------------------

quote_var <- function(x) paste0("`", x, "`")

build_formula <- function(outcome, rhs_vars) {
  outcome_q <- quote_var(outcome)
  rhs_q <- paste(quote_var(rhs_vars), collapse = " + ")
  as.formula(paste(outcome_q, "~", rhs_q, "- 1"))
}

prepare_subset <- function(df, installment_keep) {
  df %>%
    filter(country != "RDC", installment %in% installment_keep)
}

prepare_target <- function(df, installment_keep) {
  df %>%
    filter(country == "RDC", installment %in% installment_keep)
}

complete_case_filter <- function(df, vars) {
  vars <- intersect(vars, names(df))
  df %>% filter(if_all(all_of(vars), ~ !is.na(.x)))
}

varying_vars <- function(df, vars) {
  vars <- intersect(vars, names(df))
  keep <- vapply(vars, function(v) {
    x <- df[[v]]
    x <- x[is.finite(x)]
    length(unique(x)) > 1
  }, logical(1))
  vars[keep]
}

mahalanobis_distance <- function(ref_df, target_df, vars) {
  vars <- varying_vars(ref_df, vars)
  if (length(vars) < 2) stop("Mahalanobis distance requires at least two varying variables.")
  
  x_ref <- ref_df %>% dplyr::select(all_of(vars)) %>% as.matrix()
  x_tar <- target_df %>% dplyr::select(all_of(vars)) %>% as.matrix()
  
  x_ref_sc <- scale(x_ref)
  x_tar_sc <- scale(x_tar, center = attr(x_ref_sc, "scaled:center"), scale = attr(x_ref_sc, "scaled:scale"))
  
  cov_mat <- stats::cov(x_ref_sc, use = "pairwise.complete.obs")
  cov_mat <- cov_mat + diag(1e-8, ncol(cov_mat))
  
  d2 <- stats::mahalanobis(x_ref_sc, center = as.numeric(x_tar_sc[1, ]), cov = cov_mat)
  sqrt(pmax(d2, 0))
}

gower_distance <- function(ref_df, target_df, vars) {
  vars <- varying_vars(ref_df, vars)
  if (length(vars) < 2) stop("Gower distance requires at least two varying variables.")
  
  x_ref <- ref_df %>% dplyr::select(all_of(vars))
  x_tar <- target_df %>% dplyr::select(all_of(vars))
  combined <- bind_rows(x_tar, x_ref)
  
  dmat <- as.matrix(cluster::daisy(combined, metric = "gower"))
  as.numeric(dmat[1, -1])
}

kernel_weights <- function(distance, bandwidth, kernel = c("gaussian", "epanechnikov")) {
  kernel <- match.arg(kernel)
  u <- distance / bandwidth
  
  w <- if (kernel == "gaussian") {
    exp(-0.5 * u^2)
  } else {
    pmax(0, 1 - u^2)
  }
  
  w[!is.finite(w)] <- 0
  w
}

fit_wls_grid <- function(data_all,
                         installment_keep,
                         outcome,
                         rhs_vars,
                         distance_vars,
                         distance_method = c("mahalanobis", "gower"),
                         kernel = c("gaussian", "epanechnikov"),
                         n_grid = 40) {
  distance_method <- match.arg(distance_method)
  kernel <- match.arg(kernel)
  
  # Split into reference sample and RDC target.
  ref_df <- prepare_subset(data_all, installment_keep)
  target_df <- prepare_target(data_all, installment_keep)
  
  if (nrow(target_df) == 0) {
    stop("No RDC observation found for the selected installment subset.")
  }
  
  model_formula <- build_formula(outcome, rhs_vars)
  model_vars <- setdiff(all.vars(model_formula), outcome)
  
  # Keep only complete cases for the variables needed in this specification.
  ref_keep_vars <- unique(c(outcome, model_vars, distance_vars))
  tar_keep_vars <- unique(c(model_vars, distance_vars))
  ref_df <- complete_case_filter(ref_df, ref_keep_vars)
  target_df <- complete_case_filter(target_df, tar_keep_vars)
  
  if (nrow(ref_df) < length(model_vars) + 2) {
    stop("Too few complete observations after filtering for the selected specification.")
  }
  
  # Distances from each reference country to RDC.
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
  
  max_d <- max(distance_vec, na.rm = TRUE)
  if (!is.finite(max_d) || max_d <= 0) max_d <- 1
  
  # Reasonable grid around the observed distance scale.
  h_grid <- seq(from = 0.30 * max_d, to = 2.00 * max_d, length.out = n_grid)
  h_grid <- unique(pmax(h_grid, 1e-6))
  
  rmse_tbl <- vector("list", length(h_grid))
  weight_tbl <- vector("list", length(h_grid))
  
  for (j in seq_along(h_grid)) {
    h <- h_grid[j]
    if (!is.finite(h) || h <= 0) {
      rmse_tbl[[j]] <- tibble(h = h, rmse = Inf, n_effective = 0, r2 = NA_real_)
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
    
    if (nrow(reg_df) == 0) {
      rmse_tbl[[j]] <- tibble(h = h, rmse = Inf, n_effective = 0, r2 = NA_real_)
      next
    }
    
    fit <- tryCatch(
      lm(model_formula, data = reg_df, weights = .w),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      rmse_tbl[[j]] <- tibble(h = h, rmse = Inf, n_effective = nrow(reg_df), r2 = NA_real_)
      next
    }
    
    resid <- residuals(fit)
    rmse <- sqrt(sum(reg_df$.w * resid^2, na.rm = TRUE) / sum(reg_df$.w, na.rm = TRUE))
    r2 <- summary(fit)$r.squared
    
    rmse_tbl[[j]] <- tibble(h = h, rmse = rmse, n_effective = nrow(reg_df), r2 = r2)
  }
  
  rmse_tbl <- bind_rows(rmse_tbl)
  weight_tbl <- bind_rows(weight_tbl)
  best_idx <- which.min(rmse_tbl$rmse)
  best_h <- rmse_tbl$h[best_idx]
  
  ref_df_w <- ref_df %>% mutate(.w = kernel_weights(distance_vec, best_h, kernel = kernel))
  best_model <- lm(model_formula, data = ref_df_w, weights = .w)
  
  rdc_pred <- tryCatch(
    as.numeric(predict(best_model, newdata = target_df)),
    error = function(e) NA_real_
  )
  
  list(
    ref_df = ref_df,
    target_df = target_df,
    model_formula = model_formula,
    distance_method = distance_method,
    kernel = kernel,
    h_grid = h_grid,
    rmse_table = rmse_tbl,
    weight_table = weight_tbl,
    best_h = best_h,
    best_rmse = rmse_tbl$rmse[best_idx],
    best_model = best_model,
    best_weights = ref_df_w$.w,
    distance_vec = distance_vec,
    rdc_pred = rdc_pred
  )
}

print_optimal_model <- function(estimation_object, title) {
  cat("
=====================================================
")
  cat(title, "
")
  cat("Best h:", round(estimation_object$best_h, 6),
      " | RMSE:", round(estimation_object$best_rmse, 6),
      " | R-squared:", round(summary(estimation_object$best_model)$r.squared, 6), "
")
  print(estimation_object$model_formula)
  cat("
Coefficient table (Estimate, Robust Std. Error, t value, Pr(>|t|)):
")
  robust_vcov <- sandwich::vcovHC(estimation_object$best_model, type = "HC1")
  print(lmtest::coeftest(estimation_object$best_model, vcov. = robust_vcov))
  cat("=====================================================
")
}

make_similarity_barplot <- function(estimation_object,
                                    title_prefix = "Similarity to RDC",
                                    top_n = 10) {
  weights <- kernel_weights(
    distance = estimation_object$distance_vec,
    bandwidth = estimation_object$best_h,
    kernel = estimation_object$kernel
  )
  
  weights_df <- data.frame(
    country = estimation_object$ref_df$country,
    weight = weights,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(weight))
  
  top_df <- weights_df %>%
    slice_head(n = top_n) %>%
    mutate(country = factor(country, levels = rev(country)))
  
  bottom_df <- weights_df %>%
    slice_tail(n = top_n) %>%
    arrange(weight) %>%
    mutate(country = factor(country, levels = rev(country)))
  
  p_top <- ggplot(top_df, aes(x = country, y = weight)) +
    geom_col(fill = "gray80", color = "black", width = 0.7) +
    coord_flip() +
    labs(
      title = "Top-10 (pays les plus similaires)",
      subtitle = title_prefix,
      x = NULL,
      y = "Poids"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  
  p_bottom <- ggplot(bottom_df, aes(x = country, y = weight)) +
    geom_col(fill = "gray80", color = "black", width = 0.7) +
    coord_flip() +
    labs(
      title = "Bottom-10 (pays les moins similaires)",
      subtitle = title_prefix,
      x = NULL,
      y = "Poids"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  
  p_top + p_bottom + plot_layout(ncol = 2)
}

# ---------------------------------------------------------
# 3) Run the estimations
# ---------------------------------------------------------
# The store is nested as:
# estimation_store[[installment_group]][[outcome]][[distance_method]][[kernel]]
estimation_store <- list()
estimation_summary <- list()

for (inst_name in names(installment_specs)) {
  estimation_store[[inst_name]] <- list()
  
  for (outcome in outcome_vars) {
    estimation_store[[inst_name]][[outcome]] <- list()
    
    for (dist_m in distance_methods) {
      estimation_store[[inst_name]][[outcome]][[dist_m]] <- list()
      
      for (kernel in kernels) {
        message(
          "Estimating | installment = ", inst_name,
          " | outcome = ", outcome,
          " | distance = ", dist_m,
          " | kernel = ", kernel
        )
        
        est_obj <- fit_wls_grid(
          data_all = clean_df,
          installment_keep = installment_specs[[inst_name]],
          outcome = outcome,
          rhs_vars = base_rhs_vars,
          distance_vars = distance_vars,
          distance_method = dist_m,
          kernel = kernel,
          n_grid = 40
        )
        
        estimation_store[[inst_name]][[outcome]][[dist_m]][[kernel]] <- est_obj
        
        estimation_summary[[length(estimation_summary) + 1]] <- tibble(
          installment_group = inst_name,
          outcome = outcome,
          distance_method = dist_m,
          kernel = kernel,
          best_h = est_obj$best_h,
          best_rmse = est_obj$best_rmse,
          r_squared = summary(est_obj$best_model)$r.squared,
          adj_r_squared = summary(est_obj$best_model)$adj.r.squared,
          rdc_pred = est_obj$rdc_pred,
          n_reference = nrow(est_obj$ref_df),
          n_target = nrow(est_obj$target_df)
        )
      }
    }
  }
}

estimation_summary <- bind_rows(estimation_summary)
if (interactive()) View(estimation_summary)

# ---------------------------------------------------------
# 4) Display optimal equations before plotting
# ---------------------------------------------------------
for (inst_name in names(estimation_store)) {
  for (outcome in names(estimation_store[[inst_name]])) {
    for (dist_m in names(estimation_store[[inst_name]][[outcome]])) {
      for (kernel in names(estimation_store[[inst_name]][[outcome]][[dist_m]])) {
        print_optimal_model(
          estimation_store[[inst_name]][[outcome]][[dist_m]][[kernel]],
          title = paste(
            "Optimal model | installment =", inst_name,
            "| outcome =", outcome,
            "| distance =", dist_m,
            "| kernel =", kernel
          )
        )
      }
    }
  }
}

# ---------------------------------------------------------
# 5) Similarity plots (all combinations)
# ---------------------------------------------------------
# Wide weight table per bandwidth
weights_df <- estimation_store[["A"]][["yas_spread"]][["mahalanobis"]][["gaussian"]]$weight_table
weights_wide <- weights_df %>%
  mutate(h = paste0("h_", round(h, 4))) %>%
  tidyr::pivot_wider(names_from = h, values_from = weight)
View(weights_wide)

similarity_plots <- list()

for (inst_name in names(estimation_store)) {
  for (outcome in names(estimation_store[[inst_name]])) {
    for (dist_m in names(estimation_store[[inst_name]][[outcome]])) {
      for (kernel in names(estimation_store[[inst_name]][[outcome]][[dist_m]])) {
        est_obj <- estimation_store[[inst_name]][[outcome]][[dist_m]][[kernel]]
        
        plot_title <- paste0(
          tools::toTitleCase(inst_name),
          " | ", outcome,
          " | ", dist_m,
          " | ", kernel
        )
        
        p <- make_similarity_barplot(
          est_obj,
          title_prefix = plot_title
        )
        
        plot_name <- paste0(
          "similarity_bars_",
          inst_name, "_", outcome, "_", dist_m, "_", kernel, ".pdf"
        )
        
        ggsave(
          filename = file.path("plots", plot_name),
          plot = p,
          width = 16,
          height = 8,
          units = "in"
        )
        
        similarity_plots[[paste(inst_name, outcome, dist_m, kernel, sep = "|")]] <- p
      }
    }
  }
}

# ---------------------------------------------------------
# 6) Optional outputs for downstream use
# ---------------------------------------------------------
# Objects available in the workspace:
#   - estimation_store
#   - estimation_summary
#   - similarity_plots
# Save estimation summary as CSV
write.csv(estimation_summary, "estimation_summary.csv")

# =========================================================
# End of script
# =========================================================