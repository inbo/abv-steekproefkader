library(terra)
library(sf)
library(dplyr)

# 1. Load the data
grts <- rast("data/grts.tiff")
strata <- rast("data/strata.tif")
buffer <- vect("data/road_path_small_buffer.gpkg")

# Ensure CRS matches and combine into one stack
crop(grts, strata) |>
  c(strata) -> stacked
names(stacked) <- c("grts", "stratum")

# 2. Extract values AND coordinates simultaneously
# This guarantees row alignment and only extracts cells inside the polygon.
# It returns a dataframe with columns: ID, grts, stratum, x, y
project(buffer, crs(grts)) |>
  extract(x = stacked, xy = TRUE) |>
  mutate(
    stratum = factor(
      .data$stratum,
      levels = c(1, 2, 4, 5, 6, 7, 8, 11),
      labels = c(
        "urbaan",
        "suburbaan",
        "akker",
        "grasland",
        "bos",
        "heide",
        "moeras",
        "mix"
      )
    )
  ) |>
  filter(!is.na(.data$grts), !is.na(.data$stratum)) |>
  slice_min(.data$grts, by = "stratum", n = 3 * 100) |>
  st_as_sf(coords = c("x", "y"), crs = crs(grts)) |>
  group_by(.data$stratum) |>
  arrange(.data$grts) |>
  mutate(
    volgorde = row_number() - 1,
    set = .data$volgorde %/% 50
  ) |>
  ungroup() -> steekproef
st_write(steekproef, "data/steekproef.gpkg", delete_dsn = TRUE)
st_write(steekproef, "data/steekproef.shp", delete_dsn = TRUE)


library(terra)
library(sf)
library(dplyr)
grts <- rast("data/grts.tiff")
strata <- ifel(rast("data/strata.tif") == 8, 1, NA)
crop(grts, strata) |>
  c(strata) -> stacked
names(stacked) <- c("grts", "stratum")

stacked |>
  as.data.frame(xy = TRUE, na.rm = TRUE) |>
  slice_min(.data$grts, n = 3 * 100) |>
  st_as_sf(coords = c("x", "y"), crs = crs(stacked)) -> moeras
st_write(moeras, "data/moeras.gpkg", delete_dsn = TRUE)
st_write(moeras, "data/moeras.shp", delete_dsn = TRUE)
