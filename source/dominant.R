library(qgisprocess)
qgis_configure(use_cached_data = TRUE)

type <- "niv1_estat"
year <- 2022
version <- "v12"

weight_file <- "data/gaussian_weights.txt"

working_dir <- "data/temp"
dir.create(working_dir, showWarnings = FALSE)
output_dir <- "data/dominant"
dir.create(output_dir, showWarnings = FALSE)

source_raster <- sprintf("data/ecosystem/%s_%i_%s.tiff", type, year, version)
stopifnot(file.exists(source_raster))

# definitions for reclassification and focal processing
definitions <- list(
  settlement = rbind(c(0, 0, 1), c(0, 1, 32500), c(1, 12, 1)),
  cropland = rbind(c(0, 1, 1), c(1, 2, 32500), c(2, 12, 1)),
  grassland = rbind(c(0, 2, 1), c(2, 3, 32500), c(3, 12, 1)),
  forest = rbind(c(0, 3, 1), c(3, 4, 32500), c(4, 12, 1)),
  heathland = rbind(c(0, 4, 1), c(4, 5, 32500), c(5, 12, 1)),
  sparse = rbind(c(0, 5, 1), c(5, 6, 32500), c(6, 12, 1)),
  wetland = rbind(c(0, 6, 1), c(6, 7, 32500), c(7, 12, 1)),
  river = rbind(c(0, 7, 1), c(7, 8, 32500), c(8, 12, 1)),
  lake = rbind(c(0, 8, 1), c(8, 9, 32500), c(9, 12, 1)),
  inlet = rbind(c(0, 9, 1), c(9, 10, 32500), c(10, 12, 1)),
  dune = rbind(c(0, 10, 1), c(10, 11, 32500), c(11, 12, 1))
)

read.table(weight_file) |>
  nrow() -> window_size

for (landuse in names(definitions)) {
  # convert the land use raster to a binary raster for the land use of interest
  landuse_subset <- sprintf(
    "%s/%s_%s_%s_%i_subset.tiff",
    working_dir,
    landuse,
    type,
    version,
    year
  )
  if (!file.exists(landuse_subset)) {
    message("Processing land use: ", landuse)
    qgis_run_algorithm_p(
      source_raster,
      algorithm = "native:reclassifybytable",
      OUTPUT = landuse_subset,
      RASTER_BAND = 1,
      TABLE = definitions[[landuse]],
      NO_DATA = 0,
      DATA_TYPE = 1,
      RANGE_BOUNDARIES = 0,
      NODATA_FOR_MISSING = TRUE,
      .quiet = TRUE
    )
  }
  # calculate the focal average for the land use of interest
  landuse_focal <- sprintf(
    "%s/%s_%s_%s_%i_focal.tiff",
    working_dir,
    landuse,
    type,
    version,
    year
  )
  if (!file.exists(landuse_focal)) {
    message("Processing focal for land use: ", landuse)
    qgis_run_algorithm_p(
      landuse_subset,
      algorithm = "grass:r.neighbors",
      selection = landuse_subset,
      method = 0,
      size = window_size,
      weight = weight_file,
      output = landuse_focal,
      .quiet = TRUE
    )
  }
}

# calculate dominant landuse
dominant_amount <- sprintf(
  "%s/dominant_amount_%s_%s_%i.tiff",
  output_dir,
  type,
  version,
  year
)
dominant_which <- sprintf(
  "%s/dominant_which_%s_%s_%i.tiff",
  output_dir,
  type,
  version,
  year
)
qgis_run_algorithm(
  algorithm = "grass:r.mapcalc.simple",
  output = dominant_amount,
  a = sprintf(
    "%s/%s_%s_%s_%i_focal.tiff",
    working_dir,
    names(definitions)[1],
    type,
    version,
    year
  ),
  b = sprintf(
    "%s/%s_%s_%s_%i_focal.tiff",
    working_dir,
    names(definitions)[2],
    type,
    version,
    year
  ),
  expression = "A * (A > B) + B * (B >= A)"
)
qgis_run_algorithm(
  algorithm = "grass:r.mapcalc.simple",
  output = dominant_which,
  a = sprintf(
    "%s/%s_%s_%s_%i_focal.tiff",
    working_dir,
    names(definitions)[1],
    type,
    version,
    year
  ),
  b = sprintf(
    "%s/%s_%s_%s_%i_focal.tiff",
    working_dir,
    names(definitions)[2],
    type,
    version,
    year
  ),
  expression = "(A > B) + 2 * (B >= A)"
)
for (i in tail(seq_along(definitions), -2)) {
  landuse <- names(definitions)[i]
  message("Processing dominant for land use: ", landuse)
  qgis_run_algorithm(
    algorithm = "grass:r.mapcalc.simple",
    output = dominant_which,
    a = sprintf(
      "%s/%s_%s_%s_%i_focal.tiff",
      working_dir,
      landuse,
      type,
      version,
      year
    ),
    b = dominant_amount,
    c = dominant_which,
    expression = paste(i, "* (A > B) + C * (B >= A)")
  )
  qgis_run_algorithm(
    algorithm = "grass:r.mapcalc.simple",
    output = dominant_amount,
    a = dominant_amount,
    b = sprintf(
      "%s/%s_%s_%s_%i_focal.tiff",
      working_dir,
      landuse,
      type,
      version,
      year
    ),
    expression = "A * (A > B) + B * (B >= A)"
  )
}
dominant_strata <- sprintf(
  "%s/dominant_strata_%s_%s_%i.tiff",
  output_dir,
  type,
  version,
  year
)
# calculate dominant stratum
qgis_run_algorithm(
  algorithm = "grass:r.mapcalc.simple",
  output = dominant_strata,
  a = dominant_amount,
  b = dominant_which,
  expression = "B * 10 + (A > 21667) + (A > 16250) + (A > 10833)"
)
