library(n2kanalysis)
library(future)
library(furrr)
library(n2kanalysis)
library(odbc)
library(aws.s3)
library(inbodb)
library(FactoMineR)
library(designpower)
library(rvest)
library(dplyr)
library(purrr)
library(stringr)
library(factoextra)
library(tidyverse)
library(dplyr)
library(INLA)

variability_comparison <- function(model) {
  prec_marg_square <-
    model@Model$marginals.hyperpar[["Precision for square"]]

  prec_marg_point <-
    model@Model$marginals.hyperpar[["Precision for point"]]

  sd_square <- inla.emarginal(
    function(tau) 1 / sqrt(tau),
    prec_marg_square
  )

  sd_point <- inla.emarginal(
    function(tau) 1 / sqrt(tau),
    prec_marg_point
  )

  var_square <- inla.emarginal(
    function(tau) 1 / tau,
    prec_marg_square
  )

  var_point <- inla.emarginal(
    function(tau) 1 / tau,
    prec_marg_point
  )
  # Ratio of expected values != expected value of ratio of variances, just icc approximation
  ICC <- var_square / (var_square + var_point)

  return(
    c(
      sd_square = sd_square,
      sd_point = sd_point,
      var_square = var_square,
      var_point = var_point,
      ICC = ICC
    )
  )
}

variability_comparison_new_strata <- function(model) {
  prec_marg_square <-
    model$marginals.hyperpar[["Precision for square"]]

  prec_marg_point <-
    model$marginals.hyperpar[["Precision for point"]]

  sd_square <- inla.emarginal(
    function(tau) 1 / sqrt(tau),
    prec_marg_square
  )

  sd_point <- inla.emarginal(
    function(tau) 1 / sqrt(tau),
    prec_marg_point
  )

  var_square <- inla.emarginal(
    function(tau) 1 / tau,
    prec_marg_square
  )

  var_point <- inla.emarginal(
    function(tau) 1 / tau,
    prec_marg_point
  )
  # Ratio of expected values != expected value of ratio of variances, just icc approximation
  ICC <- var_square / (var_square + var_point)

  return(
    c(
      sd_square = sd_square,
      sd_point = sd_point,
      var_square = var_square,
      var_point = var_point,
      ICC = ICC
    )
  )
}

get_scenario_ests <- function(model) {
  stratum <- c(
    "Bos",
    "Heide en duin",
    "Landbouw",
    "Moeras en water",
    "Suburbaan",
    "Urbaan"
  )
  fix <- setNames(
    model@Model$summary.fixed$mean,
    rownames(model@Model$summary.fixed)
  )
  return(setNames(fix[paste0("stratum", stratum)], stratum))
}

get_scenario_ests_new_strata <- function(model) {
  stratum <- c(
    "akker",
    "bebouwd",
    "bebouwed_groen",
    "bos",
    "duinen/heide",
    "estuarium",
    "gemengd",
    "getijdensgebied/overgangswater",
    "grasland",
    "moeras",
    "water"
  )
  fix <- setNames(
    model$summary.fixed$mean,
    rownames(model$summary.fixed)
  )
  return(setNames(fix[paste0("new_stratum", stratum)], stratum))
}


data_generation <- function(counts_per_year, n_years, beta_0, beta_1, sigma) {
  df <- expand.grid(
    telpunt = seq_len(counts_per_year),
    jaar = seq_len(n_years) - 1
  )

  rf_telpunt <- rnorm(
    counts_per_year,
    mean = 0,
    sd = sigma
  )

  eta <- beta_0 + beta_1 * df$jaar + rf_telpunt
  df$y <- rpois(length(eta), lambda = exp(eta))
  return(df)
}

data_generation_rotational <- function(
  n_locations_total,
  n_years,
  beta_0,
  beta_1,
  sigma
) {
  df <- data.frame(
    telpunt = seq_len(n_locations_total),
    rotation = rep(1:3, length.out = n_locations_total),
    jaar = unlist(sapply(
      (3 * (1:(n_years / 3))),
      function(x) {
        rep(1:3, length.out = n_locations_total) - 3 + x
      },
      simplify = FALSE
    ))
  )

  rf_telpunt <- rnorm(n_locations_total, mean = 0, sd = sigma)
  df$rf <- rf_telpunt[df$telpunt]
  df$jaar <- df$jaar - 1
  # add the random effects in the log scale and then exponentiate to lambda for the rpois function
  eta <- beta_0 + beta_1 * df$jaar + df$rf
  df$y <- rpois(nrow(df), lambda = exp(eta))

  return(df)
}


fit_model_extract_params <- function(df, n_locations_total) {
  model_fit <- glmer(
    y ~ jaar + (1 | telpunt),
    data = df,
    family = poisson
  )

  beta_jaar <- fixef(model_fit)["jaar"]
  se_jaar <- sqrt(vcov(model_fit)["jaar", "jaar"])
  pval <- 2 * pnorm(-abs(beta_jaar / se_jaar))

  return(c(
    beta = beta_jaar,
    pval = pval,
    n_locations_total = n_locations_total
  ))
}


power_analysis_design <- function(design, n_sim, ...) {
  p_vals <- replicate(n_sim, {
    df <- data_generation_rotational(
      n_locations_total = design$n_locations_total,
      n_years = design$n_years,
      beta_0 = design$beta_0,
      beta_1 = design$beta_1,
      sigma = design$sigma
    )

    res <- fit_model_extract_params(df, design$n_locations_total)
    return(res["pval.jaar"])
  })
  return(list(p = unname(p_vals)))
}

signorini_marginal <- function(
  phi,
  beta_0,
  beta_1,
  sigma,
  sampled_years,
  alpha,
  power
) {
  Z_alpha <- qnorm(1 - alpha / 2)
  Z_power <- qnorm(power)
  V_x <- var(sampled_years)
  n_obs_per_point <- length(sampled_years)
  lambda_bar <- mean(exp(beta_0 + beta_1 * sampled_years + sigma^2 / 2))
  N <- (phi * (Z_alpha + Z_power)^2) /
    (n_obs_per_point * lambda_bar * (beta_1^2) * V_x)
  return(N)
}
