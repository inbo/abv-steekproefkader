library(httr)
library(osmextract)
library(qgisprocess)
library(sf)

# download border of Flanders from the WFS service of the Flemish government
# will be used to clip the Belgian OSM data
url <- parse_url("https://geo.api.vlaanderen.be/VRBG/wfs")
url$query <- list(
  service = "wfs",
  version = "2.0.0",
  request = "GetFeature",
  typename = "Refgew",
  crsName = "EPSG:31370"
)
build_url(url) |>
  read_sf() -> flanders

# check if we have the most recent Belgium OSM data, if not download it
oe_download_directory() |>
  file.path("geofabrik_belgium-latest.osm.pbf") -> osm_target
difftime(Sys.time(), file.info(osm_target)$mtime, units = "hours") |>
  as.numeric() -> delta
force_download <- delta > 24
osm_source <- oe_match("Belgium")
osm_pbf <- oe_download(
  file_url = osm_source$url,
  file_size = osm_source$file_size,
  force = force_download
)

# select water-like features from the OSM data
# to be used to clip the Belgian OSM data
# only relevant for detecting the large region without roads
waterlike <- c(
  "natural" = "water",
  "natural" = "wetland",
  "landuse" = "aquaculture",
  "landuse" = "bassin",
  "landuse" = "harbour",
  "landuse" = "lock_gate",
  "landuse" = "reservoir",
  "landuse" = "water_storage",
  "leisure" = "marina",
  "leisure" = "swimming_pool",
  "leisure" = "slipway",
  "leisure" = "fishing",
  "leisure" = "water_park"
)
# extract the relevant water polygons from OSM
oe_vectortranslate(
  osm_target,
  layer = "multipolygons",
  force_vectortranslate = force_download
) |>
  paste0("|layername=multipolygons") |>
  setNames("INPUT") |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyexpression",
    EXPRESSION = sprintf(
      "(\"%s\" = '%s')", # nolint: quotes_lintr
      names(waterlike),
      waterlike
    ) |>
      paste(collapse = " OR "),
    OUTPUT = "data/waterlike.gpkg",
    FAIL_OUTPUT = NULL
  ) -> water_selection
# extract the relevant road lines from OSM
oe_vectortranslate(
  osm_target,
  layer = "lines",
  force_vectortranslate = force_download
) |>
  paste0("|layername=lines") |>
  setNames("INPUT") |>
  qgis_run_algorithm_p(
    algorithm = "native:retainfields",
    OUTPUT = qgis_tmp_vector(),
    FIELDS = c("osm_id", "highway")
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    VALUE = NULL,
    OPERATOR = "is not null",
    OUTPUT = qgis_tmp_vector(),
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "proposed",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "planned",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "construction",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "destroyed",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "razed",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "motorway",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "motorway_link",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "emergency_bay",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "services",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "restarea",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "corridor",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "no",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:extractbyattribute",
    FIELD = "highway",
    OPERATOR = "≠",
    OUTPUT = qgis_tmp_vector(),
    VALUE = "raceway",
    FAIL_OUTPUT = NULL
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:reprojectlayer",
    OUTPUT = qgis_tmp_vector(),
    TARGET_CRS = paste(
      "PROJ4:+proj=lcc +lat_0=90 +lon_0=4.36748666666667",
      "+lat_1=51.1666672333333 +lat_2=49.8333339 +x_0=150000.013",
      "+y_0=5400088.438 +ellps=intl",
      "+towgs84=-99.059,53.322,-112.486,0.419,-0.83,1.885,-1",
      "+units=m +no_defs"
    )
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:clip",
    OVERLAY = flanders,
    OUTPUT = "data/road_path.gpkg"
  ) -> road_path
# create a small buffer around the road paths to be used for clipping the GRTS
# sample
road_path |>
  qgis_run_algorithm_p(
    algorithm = "native:buffer",
    DISTANCE = 10,
    DISSOLVE = FALSE,
    END_CAP_STYLE = 0,
    JOIN_STYLE = 0,
    MITER_LIMIT = 2,
    SEGMENTS = 5,
    OUTPUT = "data/road_path_small_buffer.gpkg",
    SEPARATE_DISJOINT = TRUE
  ) -> small_buffer
qgis_run_algorithm_p(small_buffer, algorithm = "native:createspatialindex")
qgis_run_algorithm_p(
  algorithm = "grass:r.mapcalc.simple",
  a = "data/strata.tif",
  expression = paste("A ==", c(1, 2, 4, 5, 6, 7, 8, 11), collapse = " || ") |>
    sprintf(fmt = "if(%s, A)"),
  OUTPUT = qgis_tmp_vector()
)

# create a wider buffer (200 meters) to detect locations that are not observable
# from the roads
road_path |>
  qgis_run_algorithm_p(
    algorithm = "native:buffer",
    DISTANCE = 200,
    DISSOLVE = TRUE,
    END_CAP_STYLE = 0,
    JOIN_STYLE = 0,
    MITER_LIMIT = 2,
    SEGMENTS = 5,
    OUTPUT = "data/road_path_buffer.gpkg",
    SEPARATE_DISJOINT = TRUE
  ) -> road_buffer
# the erosion and dilution removes small regions between the roads that are not
# observable from the roads
flanders |>
  # clip the Flanders border with the road buffer
  qgis_run_algorithm_p(
    algorithm = "native:difference",
    OVERLAY = road_buffer$OUTPUT,
    OUTPUT = qgis_tmp_vector()
  ) |>
  # ignore water polygons
  qgis_run_algorithm_p(
    algorithm = "native:difference",
    OVERLAY = water_selection$OUTPUT,
    OUTPUT = qgis_tmp_vector()
  ) |>
  # erosion
  qgis_run_algorithm_p(
    algorithm = "native:buffer",
    DISTANCE = -50,
    END_CAP_STYLE = 0,
    JOIN_STYLE = 0,
    MITER_LIMIT = 2,
    SEGMENTS = 5,
    OUTPUT = qgis_tmp_vector(),
    SEPARATE_DISJOINT = TRUE
  ) |>
  # dilation
  qgis_run_algorithm_p(
    algorithm = "native:buffer",
    DISTANCE = 50,
    DISSOLVE = TRUE,
    END_CAP_STYLE = 0,
    JOIN_STYLE = 0,
    MITER_LIMIT = 2,
    SEGMENTS = 5,
    OUTPUT = qgis_tmp_vector(),
    SEPARATE_DISJOINT = TRUE
  ) |>
  qgis_run_algorithm_p(
    algorithm = "native:multiparttosingleparts",
    OUTPUT = qgis_tmp_vector()
  ) |>
  # calculate the area
  qgis_run_algorithm_p(
    algorithm = "native:fieldcalculator",
    FIELD_NAME = "ha",
    FORMULA = "$area / 10000",
    FIELD_TYPE = "Decimal (double)", # nolint
    OUTPUT = "data/not_covered.gpkg"
  ) -> road_buffer_clipped
