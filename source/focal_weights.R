library(mvtnorm)
library(tidyverse)

sigma <- 7.77
window_size <- 41
pixel_size <- 10
radius <- 200
weight_file <- "data/gaussian_weights.txt"

# check sigma and radius are consistent with each other
quant <- qmvnorm(
  p = 0.99,
  tail = "lower.tail",
  sigma = diag(2) * sigma^2
)$quantile
stopifnot(
  "mismatch between `sigma` and `radius`" = abs(quant * pixel_size - radius) <=
    0.1
)

# calculate the focal weights
-floor(window_size / 2) |>
  seq(floor(window_size / 2)) -> offsets
list(offsets) |>
  rep(2) |>
  expand.grid() |>
  mutate(
    dens = cbind(.data$Var1, .data$Var2) |>
      dmvnorm(sigma = sigma^2 * diag(2)),
    across(c("Var1", "Var2"), ~ .x * pixel_size),
  ) |>
  pivot_wider(names_from = "Var2", values_from = "dens") |>
  arrange(.data$Var1) |>
  select(-"Var1") |>
  as.matrix() |>
  write.table(
    file = weight_file,
    row.names = FALSE,
    col.names = FALSE
  )
