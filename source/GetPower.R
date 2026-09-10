source("./source/AnalysisFunctions.R")

vars_df <- readRDS(
  "./data/VogelModelParameters.rds"
)

strata <- c(
  "Bos",
  "Heide.en.duin",
  "Landbouw",
  "Moeras.en.water",
  "Suburbaan",
  "Urbaan"
)
average_est <- mean(colMeans(vars_df[, strata], na.rm = TRUE))
mean_sigma <- mean(vars_df$sigma, na.rm = TRUE)
summary(rowMeans(vars_df[, strata], na.rm = TRUE))
plot(vars_df[, strata])
plot(sort(unlist(vars_df[, strata])))

design <- list(
  n_locations_total = 21,
  n_years = 24,
  beta_0 = average_est,
  beta_1 = 0.009352669,
  sigma = mean_sigma
)

object_testing <- designpower::find_power(
  design = design,
  design_digits = c(
    n_locations_total = 0,
    n_years = 0,
    beta_0 = 0,
    beta_1 = 4,
    sigma = 0
  ),
  opti = "n_locations_total",
  sim_power = power_analysis_design,
  extra_args = list(),
  power = 0.9,
  alpha = 0.1,
  n_sim = 10,
  max_sim = 500,
  filename = "testing_averagebeta_averagesigma_correct.duckdb"
)
