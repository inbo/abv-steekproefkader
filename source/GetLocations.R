# Script to check the locations of the birds from ABV
# What is the occupancy of each bird
# What would the strata distribution be if we use ABV2.0 but with the original points

# Read in the location data from the database
library(ggplot2)
library(rvest)
library(dplyr)
library(purrr)
library(stringr)
library(odbc)
library(n2kanalysis)
library(inbodb)
project <- "abv"
#base <- keyring::key_get("n2kbucket") |>
#  aws.s3::get_bucket(, prefix = project, max = 1)

con <- connect_inbo_dbase("S0008_00_Meetnetten")

library(dplyr)
library(DBI)

# 1. Query only necessary fields and extract Lon/Lat directly
query_coords <- "
SELECT 
    l.name AS point,
    lp.name AS square,
    l.geom.STX AS longitude,
    l.geom.STY AS latitude
FROM [staging_Meetnetten].[projects_projectlocation] pl
INNER JOIN [staging_Meetnetten].[locations_location] l ON l.id = pl.location_id
INNER JOIN [staging_Meetnetten].[locations_location] lp ON lp.id = l.parent_id
INNER JOIN [staging_Meetnetten].[projects_project] p ON p.id = pl.project_id
WHERE p.name = 'Algemene Broedvogelmonitoring (ABV)'
"

coords_df <- DBI::dbGetQuery(con, query_coords) %>%
  distinct()


# Get the coordinates from location name
# read through all the models and create a dataframe from those as they have the strata
# link both together and then i will have the original strata linked to the locations.
base_url <- "https://www.vlaanderen.be/inbo/rapporten/abv-trends-2007-2025/"

overview <- read_html(paste0(base_url, "soort/trends.html"))
species_links <- overview |>
  html_elements("a[href*='/soort/']") |>
  {
    \(x) {
      tibble(
        naam = html_text2(x),
        url = html_attr(x, "href")
      )
    }
  }() |>
  filter(!naam %in% c("", "Overzicht")) |>
  mutate(
    url = ifelse(
      str_starts(url, "http"),
      url,
      paste0(base_url, "soort/", basename(url))
    )
  ) |>
  distinct(url, .keep_all = TRUE)
get_hash_table <- function(url, naam) {
  tables <- url |> read_html() |> html_table()
  hash_tbl <- keep(tables, ~ all(c("analyse", "status") %in% names(.x)))
  if (length(hash_tbl) == 0) {
    return(NULL)
  }
  hash_tbl[[1]] |>
    mutate(
      soort = naam,
      analyse = str_remove_all(analyse, "\\s"),
      status = str_remove_all(status, "\\s")
    )
}

hash_data <- map2_dfr(
  species_links$url,
  species_links$naam,
  possibly(get_hash_table, otherwise = NULL),
  .progress = TRUE
)

all_hashes <- hash_data |>
  dplyr::filter(
    frequentie == "jaarlijks" & model == "lineair"
  )

# loop over the analysis
# download each one if it has not been downloaded yet
failed_soorten <- c()
vars_list <- vector(mode = "list", length = nrow(all_hashes))
names(vars_list) <- all_hashes$soort

for (i in 1:nrow(all_hashes)) {
  soort_name <- all_hashes$soort[i]
  print(soort_name)

  file_name <- paste0(
    "./source/Nieuwe map/PowerAnalyse/Data/",
    soort_name,
    ".rds"
  )
  if (!file.exists(file_name)) {
    tmp <- read_model(x = all_hashes$analyse[i], base = base, project = project)
    saveRDS(tmp, file = file_name)
  }

  model <- readRDS(file_name)
  vars_list[[soort_name]] <- model@Data[, c(
    "stratum",
    "square",
    "point",
    "weight"
  )] |>
    unique() |>
    mutate(species = soort_name)
}
print(failed_soorten)
vars_df <- do.call(rbind, vars_list)
#TODO: check where the nas come frome
vars_df <- na.omit(vars_df)

vars_df_updated <- vars_df %>%
  left_join(coords_df, by = c("square", "point"))
library(sf)

pt_wgs84 <- st_as_sf(
  vars_df_updated[, c("longitude", "latitude")] |> unique(),
  coords = c("longitude", "latitude"),
  crs = 4326
)

pt_lambert72 <- st_transform(pt_wgs84, 31370)

st_coordinates(pt_lambert72)


library(terra)
map_file <- "./data/gauss100m.tif"
flanders_raster <- rast(map_file)
nieuwe_strata <- c(
  "bebouwd",
  "bebouwed_groen",
  "water",
  "akker",
  "grasland",
  "bos",
  "duinen/heide",
  "moeras",
  "estuarium",
  "getijdensgebied/overgangswater",
  "gemengd"
)

extracted_vals <- terra::extract(flanders_raster, pt_lambert72, xy = TRUE)
table(nieuwe_strata[extracted_vals$max])
freq(flanders_raster)


unique_pts_sf <- vars_df_updated %>%
  distinct(longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(31370)

extracted_vals <- terra::extract(flanders_raster, unique_pts_sf, xy = TRUE)

raster_lookup <- unique_pts_sf %>%
  st_drop_geometry() %>%
  mutate(
    x_lambert = extracted_vals$x,
    y_lambert = extracted_vals$y,
    new_stratum_id = extracted_vals$max,
    new_stratum = nieuwe_strata[extracted_vals$max]
  )

# 4. Join back to full dataset and select desired columns
final_df <- vars_df_updated %>%
  left_join(raster_lookup, by = c("longitude", "latitude")) %>%
  select(
    stratum_old = stratum,
    square,
    point,
    species,
    longitude,
    latitude,
    x_lambert,
    y_lambert,
    new_stratum_id,
    new_stratum
  )


strata_stats_flanders_df <- freq(flanders_raster)
strata_stats_flanders_df$new_stratum <- nieuwe_strata[
  strata_stats_flanders_df$value
]

strata_stats_flanders_df$A_s <- strata_stats_flanders_df$count /
  sum(strata_stats_flanders_df$count)

final_df <- merge(
  final_df,
  strata_stats_flanders_df[, c("count", "new_stratum", "A_s")]
)

## ---- draai dit éénmalig na je bestaande script -------------------------------
## abv2_point_classes.rds: 1 rij per uniek punt
abv2_point_classes <- final_df %>%
  distinct(point, square, stratum_old, longitude, latitude, new_stratum)

## abv2_relevance.rds: relevante punt-soort combinaties (analyse-universum)
abv2_relevance <- vars_df %>% # model@Data-gebaseerd zoals in je script
  distinct(species, point)

## abv2_raster_freq.rds: cellen per nieuwe klasse over heel Vlaanderen
abv2_raster_freq <- freq(flanders_raster) %>%
  as.data.frame() %>%
  transmute(new_stratum = nieuwe_strata[value], cellen = count)


saveRDS(
  abv2_point_classes,
  "C:/Users/stijn_vandenbulcke/Documents/abv/abv-steekproefkader/source/Nieuwe map/PowerAnalyse/Data/abv2_point_classes.rds"
)
saveRDS(
  abv2_relevance,
  "C:/Users/stijn_vandenbulcke/Documents/abv/abv-steekproefkader/source/Nieuwe map/PowerAnalyse/Data/abv2_relevance.rds"
)
saveRDS(
  abv2_raster_freq,
  "C:/Users/stijn_vandenbulcke/Documents/abv/abv-steekproefkader/source/Nieuwe map/PowerAnalyse/Data/abv2_raster_freq.rds"
)

####################################
##### SAME THING BUT FOR ROADS #####
####################################
map_file <- "C:/Users/stijn_vandenbulcke/Documents/abv/abv-steekproefkader/source/Nieuwe map/PowerAnalyse/RoadBufferBias/buffer_mask_100m.tif"
buffer_mask_100m <- rast(map_file)
flanders_raster <- mask(flanders_raster, buffer_mask_100m)
# Find each point on the new map (ignore out of bounds points)

extracted_vals <- terra::extract(flanders_raster, pt_lambert72, xy = TRUE)
table(nieuwe_strata[extracted_vals$max])
# count the strata distributions and occurance of each bird for each strata
freq(flanders_raster)

unique_pts_sf <- vars_df_updated %>%
  distinct(longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(31370)

# 2. Extract raster values and Lambert X/Y coordinates
extracted_vals <- terra::extract(flanders_raster, unique_pts_sf, xy = TRUE)

# 3. Map extracted values to unique coordinates lookup dataframe
raster_lookup <- unique_pts_sf %>%
  st_drop_geometry() %>%
  mutate(
    x_lambert = extracted_vals$x,
    y_lambert = extracted_vals$y,
    new_stratum_id = extracted_vals$max,
    new_stratum = nieuwe_strata[extracted_vals$max]
  )

# 4. Join back to full dataset and select desired columns
final_df <- vars_df_updated %>%
  left_join(raster_lookup, by = c("longitude", "latitude")) %>%
  select(
    stratum_old = stratum,
    square,
    point,
    species,
    longitude,
    latitude,
    x_lambert,
    y_lambert,
    new_stratum_id,
    new_stratum
  )


strata_stats_flanders_df <- freq(flanders_raster)
strata_stats_flanders_df$new_stratum <- nieuwe_strata[
  strata_stats_flanders_df$value
]

strata_stats_flanders_df$A_s <- strata_stats_flanders_df$count /
  sum(strata_stats_flanders_df$count)

final_df <- merge(
  final_df,
  strata_stats_flanders_df[, c("count", "new_stratum", "A_s")]
)

## ---- draai dit éénmalig na je bestaande script -------------------------------
## abv2_point_classes.rds: 1 rij per uniek punt
abv2_point_classes <- final_df %>%
  distinct(point, square, stratum_old, longitude, latitude, new_stratum)

## abv2_relevance.rds: relevante punt-soort combinaties (analyse-universum)
abv2_relevance <- vars_df %>% # model@Data-gebaseerd zoals in je script
  distinct(species, point)

## abv2_raster_freq.rds: cellen per nieuwe klasse over heel Vlaanderen
abv2_raster_freq <- freq(flanders_raster) %>%
  as.data.frame() %>%
  transmute(new_stratum = nieuwe_strata[value], cellen = count)


saveRDS(
  abv2_point_classes,
  "C:/Users/stijn_vandenbulcke/Documents/abv/abv-steekproefkader/source/Nieuwe map/PowerAnalyse/Data/abv2_point_classes_road.rds"
)
saveRDS(
  abv2_relevance,
  "C:/Users/stijn_vandenbulcke/Documents/abv/abv-steekproefkader/source/Nieuwe map/PowerAnalyse/Data/abv2_relevance_road.rds"
)
saveRDS(
  abv2_raster_freq,
  "C:/Users/stijn_vandenbulcke/Documents/abv/abv-steekproefkader/source/Nieuwe map/PowerAnalyse/Data/abv2_raster_freq_road.rds"
)
