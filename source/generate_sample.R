library(terra)
library(DBI)
library(grtsdb)

type <- "niv1-estat"
year <- 2022
version <- "v1-2"
target <- "data/grts.tiff"
source_raster <- sprintf("data/ecosystem/%s-%s-%i.tiff", type, version, year)
stopifnot(file.exists(source_raster))


# Use the source raster as the geometry template so that the output has exactly
# the same extent, resolution and CRS.
template <- rast(source_raster)
n_col <- ncol(template)
n_row <- nrow(template)
xmin(template) |>
  c(xmax(template)) |>
  mean() -> mean_x
ymin(template) |>
  c(ymax(template)) |>
  mean() -> mean_y

# The middle `x1` value (2^14) maps to the middle x coordinate and every unit
# increase in `x1` shifts one pixel to the east.  The same holds for `x2` and
# the y coordinate, but there the axis is reversed: y decreases with the row
# number, so a unit increase in `x2` shifts one pixel to the north.
# The half cell offset on x follows the convention used for the CSV export
# above.
mid_index <- 2^14
col_mid <- colFromX(template, mean_x)
row_mid <- rowFromY(template, mean_y)

# Inverse mappings between raster indices and the `level15` coordinates.
x1_from_col <- function(col) col - col_mid + mid_index
x2_from_row <- function(row) row_mid - row + mid_index
col_from_x1 <- function(x1) x1 + col_mid - mid_index
row_from_x2 <- function(x2) row_mid + mid_index - x2

# `ranking` has a maximum of 4^15 which still fits in a signed 32 bit integer.
output <- rast(template)
chunk_rows <- 512

writeStart(
  output,
  filename = target,
  overwrite = TRUE,
  datatype = "INT4S",
  NAflag = -1,
  gdal = c("COMPRESS=DEFLATE", "BIGTIFF=IF_NEEDED")
)

db <- connect_db("data/grts.sqlite")
c(xmin(template), xmax(template)) |>
  rbind(
    c(ymin(template), ymax(template))
  ) |>
  add_level(grtsdb = db, level = 15, cellsize = res(template)[1])
compact_db(db)

for (row_start in seq(1, n_row, by = chunk_rows)) {
  message(row_start, " ")
  row_end <- min(row_start + chunk_rows - 1, n_row)
  chunk <- rep(NA_integer_, (row_end - row_start + 1) * n_col)
  # `x2` decreases with increasing row number, hence the reversed bounds.
  sprintf(
    paste(
      "SELECT x1, x2, ranking FROM level15",
      "WHERE x2 BETWEEN %i AND %i AND x1 BETWEEN %i AND %i"
    ),
    x2_from_row(row_end),
    x2_from_row(row_start),
    x1_from_col(1),
    x1_from_col(n_col)
  ) |>
    dbGetQuery(conn = db) -> ranking
  if (nrow(ranking) > 0) {
    # Row major index within the chunk; `x1` and `x2` values outside the raster
    # are already excluded by the query bounds.
    chunk[
      (row_from_x2(ranking$x2) - row_start) * n_col + col_from_x1(ranking$x1)
    ] <- as.integer(ranking$ranking)
  }
  writeValues(output, chunk, start = row_start, nrows = row_end - row_start + 1)
}
writeStop(output)
