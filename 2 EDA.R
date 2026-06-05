# =========================================================
# EDA - Eurobond DRC project
# =========================================================
# Outputs saved in: ./plots
# Graphs:
#   1) Correlation plot (lower triangle: scatter, upper triangle: correlations)
#   2) Coupon vs maturity (2 panels: A / B)
#   3) Spread vs S&P rating (2 panels: A / B)
#   4) Debt by region / hipc (2 panels)
# =========================================================

library(dplyr)
library(ggplot2)
library(GGally)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(scales)
library(here)

# ---------------------------------------------------------
# 0) Expected input
# ---------------------------------------------------------
# This script assumes the raw dataset is already available
# in the R environment under the object name `raw_df`.
# If not, stop with a clear message.
if (!exists("raw_df")) {
  stop("Object `raw_df` not found. Run the preprocessing script first.")
}

# ---------------------------------------------------------
# 1) Output folder
# ---------------------------------------------------------
plots_dir <- here("plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# -----------------------------------------------------------
# 2) Light factor recoding for plotting / removing outliers
# -----------------------------------------------------------
raw_df <- raw_df %>%
  filter(!country %in% c("Liban", "Venezuela")) %>%
  mutate(
    region = factor(region),
    installment = factor(installment),
    hipc = factor(hipc, levels = c(0, 1), labels = c("Non-IPPTE", "Post-IPPTE"))
  )

# ---------------------------------------------------------
# 3) Correlation plot
# ---------------------------------------------------------
# Use numeric columns only. This includes continuous variables,
# ordinal scores, and dummies, but excludes pure categoricals.
corr_df <- raw_df %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-any_of(c("N°", ".row_id", "amount", "amount_usd", "maturity", "coupon_rate",
                   "corruption_est", "corruption_scr", "effectiveness_est",
                   "effectiveness_scr", "stability_est", "stability_scr",
                   "regulatory_est", "regulatory_scr", "law_scr",
                   "voice_est", "voice_scr")))

corr_plot <- ggpairs(
  corr_df,
  lower = list(continuous = wrap("points", alpha = 0.35, size = 0.7)),
  upper = list(continuous = wrap("cor", size = 3)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  theme_bw(base_size = 8)

ggsave(
  filename = file.path(plots_dir, "correlation_plot.pdf"),
  plot = corr_plot,
  width = 16,
  height = 16,
  units = "in"
)

# ---------------------------------------------------------
# 4) Scatter plot: coupon vs maturity
# ---------------------------------------------------------
coupon_maturity_base <- raw_df %>%
  dplyr::select(country, installment, maturity, coupon_rate) %>%
  filter(!is.na(maturity), !is.na(coupon_rate)) %>%
  filter(!country %in% c("Chine"))

# Subsets of the data
coupon_maturity_subset_A <- subset(coupon_maturity_base, installment %in% c("A", "Unique"))
coupon_maturity_subset_B <- subset(coupon_maturity_base, installment %in% c("B", "Unique"))

coupon_maturity_A <- ggplot(coupon_maturity_subset_A, aes(x = maturity, y = coupon_rate)) +
  geom_point(color = "black", size = 2.6, alpha = 0.85) +
  geom_text_repel(aes(label = country), size = 5, max.overlaps = 40, show.legend = FALSE) +
  labs(
    title = "(a) Tranche A",
    x = "Maturité (années)",
    y = "Coupon (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

coupon_maturity_B <- ggplot(coupon_maturity_subset_B, aes(x = maturity, y = coupon_rate)) +
  geom_point(color = "black", size = 2.6, alpha = 0.85) +
  geom_text_repel(aes(label = country), size = 5, max.overlaps = 40, show.legend = FALSE) +
  labs(
    title = "(b) Tranche B",
    x = "Maturité (années)",
    y = "Coupon (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

coupon_maturity_plot <- coupon_maturity_A + coupon_maturity_B

ggsave(
  filename = file.path(plots_dir, "coupon_vs_maturity.pdf"),
  plot = coupon_maturity_plot,
  width = 16,
  height = 8,
  units = "in"
)

# ---------------------------------------------------------
# 5) Scatter plot: spread vs rating
# ---------------------------------------------------------
# Use S&P numeric ordinal score: smaller = better rating.
spread_rating_base <- raw_df %>%
  dplyr::select(country, installment, sp_note, asw_spread) %>%
  filter(!is.na(sp_note), !is.na(asw_spread)) %>%
  filter(!country %in% c("Sénégal"))

# Subsets of the data
spread_rating_subset_A <- subset(spread_rating_base, installment %in% c("A", "Unique"))
spread_rating_subset_B <- subset(spread_rating_base, installment %in% c("B", "Unique"))

spread_rating_A <- ggplot(spread_rating_subset_A, aes(x = sp_note, y = asw_spread)) +
  geom_point(color = "black", size = 2.6, alpha = 0.85) +
  geom_text_repel(aes(label = country), size = 5, max.overlaps = 40, show.legend = FALSE) +
  labs(
    title = "(a) Tranche A",
    x = "S&P note (score ordinal numérique)",
    y = "Spread (bps)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

spread_rating_B <- ggplot(spread_rating_subset_B, aes(x = sp_note, y = asw_spread)) +
  geom_point(color = "black", size = 2.6, alpha = 0.85) +
  geom_text_repel(aes(label = country), size = 5, max.overlaps = 40, show.legend = FALSE) +
  labs(
    title = "(b) Tranche B",
    x = "S&P note (score ordinal numérique)",
    y = "Spread (bps)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

spread_rating_plot <- spread_rating_A + spread_rating_B

ggsave(
  filename = file.path(plots_dir, "spread_vs_rating.pdf"),
  plot = spread_rating_plot,
  width = 16,
  height = 8,
  units = "in"
)

# ---------------------------------------------------------
# 6) Box plot: debt by region / hipc
# ---------------------------------------------------------
debt_base <- raw_df %>%
  dplyr::select(country, region, hipc, debt) %>%
  filter(!is.na(debt))

p_debt_region <- ggplot(debt_base, aes(x = region, y = debt)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.6) +
  labs(
    title = "(a) Distribution par région",
    x = "",
    y = "Dette brute (% PIB)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)
  ) + theme(legend.position = "none")

p_debt_hipc <- ggplot(debt_base, aes(x = hipc, y = debt)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.6) +
  labs(
    title = "(b) Distribution par statut IPPTE",
    x = "",
    y = "Dette brute (% PIB)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)
  ) + theme(legend.position = "none")

p_debt <- p_debt_region + p_debt_hipc + plot_layout(ncol = 2)

ggsave(
  filename = file.path(plots_dir, "debt_distribution.pdf"),
  plot = p_debt,
  width = 16,
  height = 8,
  units = "in"
)

# ---------------------------------------------------------
# 7) Console confirmation
# ---------------------------------------------------------
message("EDA plots saved in: ", plots_dir)

# =========================================================
# End of script
# =========================================================
