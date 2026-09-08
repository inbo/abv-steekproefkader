# ==============================================================================
# FLANDERS NIV2 - gauss_100m - focal calculations
# ==============================================================================

library(tidyverse)
library(terra)
library(qgisprocess)
library(sf)

src_dir <- "./data/ecosystem/"
data_dir <- "./data/ecosystem/flanders/"

# Beperking geheugen gebruik, anders kon ik het script niet runnen
terraOptions(
  memfrac = 0.35,
  progress = 10,
  tempdir = "./data/ecosystem/terratemp/"
)

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

r_template <- rast(paste0(src_dir, "ecosysteemkaart_niv2_2022_v12.tiff"))

kernel_100m_B <- focalMat(
  r_template,
  d = 100 / 2.57,
  type = "Gauss",
  fillNA = TRUE
)

rm(r_template)
gc()

#todo: change the way kernel configs is read
#it was a list because i was testing different parameters and running all at once
kernel_configs <- list(
  list(
    label = "gauss_100m",
    focal_subdir = "Gauss100m/focal_intermediary/",
    kernel_mat = kernel_100m_B
  )
)

for (cfg in kernel_configs) {
  dir.create(
    paste0(data_dir, cfg$focal_subdir),
    showWarnings = FALSE,
    recursive = TRUE
  )
}

eco_dict_short <- c(
  "101" = "urbaan",
  "102" = "urbaan",
  "103" = "urbaan_groen",
  "104" = "urbaan_groen",
  "105" = "stilstaand_lopend_water",
  "106" = "urbaan",
  "201" = "landbouw",
  "202" = "landbouw",
  "203" = "landbouw",
  "204" = "landbouw",
  "205" = "landbouw",
  "301" = "grasland",
  "302" = "grasland",
  "400" = "bos",
  "401" = "bos",
  "402" = "bos",
  "403" = "bos",
  "404" = "bos",
  "501" = "duinen_heide",
  "602" = "duinen_heide",
  "701" = "moeras",
  "702" = "moeras",
  "801" = "stilstaand_lopend_water",
  "802" = "stilstaand_lopend_water",
  "901" = "stilstaand_lopend_water",
  "902" = "stilstaand_lopend_water",
  "1001" = "estuarium",
  "1002" = "getijdengebied_overgangswater",
  "1003" = "getijdengebied_overgangswater",
  "1102" = "duinen_heide",
  "1103" = "duinen_heide"
)

eco_order <- c(
  "urbaan",
  "urbaan_groen",
  "stilstaand_lopend_water",
  "landbouw",
  "grasland",
  "bos",
  "duinen_heide",
  "moeras",
  "estuarium",
  "getijdengebied_overgangswater"
)

# Segregate Niv2 into binary layers

r <- rast(paste0(src_dir, "ecosysteemkaart_niv2_2022_v12.tiff"))
seg_file <- paste0(data_dir, "segregated_temp.tif")

if (file.exists(seg_file)) {
  message(sprintf(
    "[%s] Segregated file exists — loading from disk.",
    Sys.time()
  ))
  flanders_layers <- rast(seg_file)
} else {
  message(sprintf("[%s] Segregating layers to disk...", Sys.time()))
  flanders_layers <- segregate(
    r,
    filename = seg_file,
    overwrite = TRUE,
    wopt = list(
      gdal = c(
        "COMPRESS=DEFLATE",
        "TILED=YES",
        "BLOCKXSIZE=256",
        "BLOCKYSIZE=256",
        "BIGTIFF=YES"
      )
    )
  )
  gc()
}

layer_names <- names(flanders_layers)

flanders_layers <- flanders_layers[[layer_names[
  layer_names %in% names(eco_dict_short)
]]]
layer_names <- names(flanders_layers)

# Grouping the ecosystems

eco_layers <- sapply(eco_dict_short[layer_names], function(x) {
  which(eco_order == x)
}) |>
  as.vector()
grouped_file <- paste0(data_dir, "grouped_ecosystems_temp.tif")
layer_names <- eco_order[unique(eco_layers)]

if (file.exists(grouped_file)) {
  message(sprintf("[%s] Grouped file exists — loading from disk.", Sys.time()))
  flanders_layers <- rast(grouped_file)
} else {
  message(sprintf("[%s] Grouping with tapp()...", Sys.time()))
  flanders_layers <- tapp(
    flanders_layers,
    index = eco_layers,
    fun = sum,
    filename = grouped_file,
    overwrite = TRUE,
    wopt = list(
      names = layer_names,
      gdal = c(
        "COMPRESS=DEFLATE",
        "TILED=YES",
        "BLOCKXSIZE=256",
        "BLOCKYSIZE=256",
        "BIGTIFF=YES"
      )
    )
  )
}

total_layers <- nlyr(flanders_layers)
rm(flanders_layers, r)
gc()

# Could not run the focal window on my laptop
# Below is gemini code in order to run it with tiling
# If i do it without tiling my R session crashes

MAX_ROWS <- 2000
MAX_COLS <- 2000

template <- rast(grouped_file, lyrs = 1)
total_cols <- ncol(template)
total_rows <- nrow(template)
rm(template)
gc()

TILE_NROWS <- ceiling(total_rows / MAX_ROWS)
TILE_NCOLS <- ceiling(total_cols / MAX_COLS)
tile_grid <- expand.grid(ti = seq_len(TILE_NROWS), tj = seq_len(TILE_NCOLS))
n_tiles <- nrow(tile_grid)

get_tile_window <- function(ti, tj, overlap) {
  row_start <- floor((ti - 1) * total_rows / TILE_NROWS) + 1
  row_end <- floor(ti * total_rows / TILE_NROWS)
  col_start <- floor((tj - 1) * total_cols / TILE_NCOLS) + 1
  col_end <- floor(tj * total_cols / TILE_NCOLS)

  read_row_start <- max(1, row_start - overlap)
  read_row_end <- min(total_rows, row_end + overlap)
  read_col_start <- max(1, col_start - overlap)
  read_col_end <- min(total_cols, col_end + overlap)

  crop_row_start <- row_start - read_row_start + 1
  crop_row_end <- crop_row_start + (row_end - row_start)
  crop_col_start <- col_start - read_col_start + 1
  crop_col_end <- crop_col_start + (col_end - col_start)

  list(
    read = c(read_row_start, read_row_end, read_col_start, read_col_end),
    crop = c(crop_row_start, crop_row_end, crop_col_start, crop_col_end),
    core = c(row_start, row_end, col_start, col_end)
  )
}

kernel_configs <- kernel_configs[1]
for (cfg in kernel_configs) {
  focal_dir <- paste0(data_dir, cfg$focal_subdir)
  kernel_mat <- cfg$kernel_mat
  kernel_label <- cfg$label

  # OVERLAP_CELLS derived directly from the kernel so tiling is always correct.
  OVERLAP_CELLS <- floor(nrow(kernel_mat) / 2)

  message(sprintf(
    "\n[%s] ══ Kernel: %s | overlap: %d cells | %d layers × %d tiles = %d ops ══",
    Sys.time(),
    kernel_label,
    OVERLAP_CELLS,
    total_layers,
    n_tiles,
    total_layers * n_tiles
  ))

  layer_result_files <- vector("list", total_layers)

  for (i in seq_len(total_layers)) {
    current_name <- layer_names[i]
    merged_file <- sprintf("%sfocal_%s.tif", focal_dir, current_name)

    if (file.exists(merged_file)) {
      message(sprintf(
        "\n[%s] [%s] Layer %d/%d (%s): already done — skipping.",
        Sys.time(),
        kernel_label,
        i,
        total_layers,
        current_name
      ))
      layer_result_files[[i]] <- merged_file
      next
    }

    message(sprintf(
      "\n[%s] [%s] Layer %d/%d: %s",
      Sys.time(),
      kernel_label,
      i,
      total_layers,
      current_name
    ))
    tile_result_files <- vector("list", n_tiles)

    for (t in seq_len(n_tiles)) {
      tile_crop_file <- sprintf(
        "%stemp_crop_%s_tile%02d.tif",
        focal_dir,
        current_name,
        t
      )

      if (file.exists(tile_crop_file)) {
        message(sprintf("  tile %d/%d already done — skipping.", t, n_tiles))
        tile_result_files[[t]] <- tile_crop_file
        next
      }

      ti <- tile_grid$ti[t]
      tj <- tile_grid$tj[t]
      w <- get_tile_window(ti, tj, OVERLAP_CELLS)

      message(sprintf(
        "  tile %d/%d (rows %d–%d, cols %d–%d)",
        t,
        n_tiles,
        w$read[1],
        w$read[2],
        w$read[3],
        w$read[4]
      ))

      single_layer <- rast(grouped_file, lyrs = i)
      padded_tile <- single_layer[
        w$read[1]:w$read[2],
        w$read[3]:w$read[4],
        drop = FALSE
      ]
      rm(single_layer)
      gc()

      tile_focal_file <- sprintf(
        "%stemp_focal_%s_tile%02d.tif",
        focal_dir,
        current_name,
        t
      )
      focal(
        padded_tile,
        w = kernel_mat,
        fun = "sum",
        na.rm = TRUE,
        filename = tile_focal_file,
        overwrite = TRUE,
        wopt = list(
          gdal = c(
            "COMPRESS=DEFLATE",
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "BIGTIFF=YES"
          )
        )
      )
      rm(padded_tile)
      gc()

      focal_result <- rast(tile_focal_file)
      core_tile <- focal_result[
        w$crop[1]:w$crop[2],
        w$crop[3]:w$crop[4],
        drop = FALSE
      ]

      writeRaster(
        core_tile,
        tile_crop_file,
        overwrite = TRUE,
        wopt = list(
          gdal = c(
            "COMPRESS=DEFLATE",
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "BIGTIFF=YES"
          )
        )
      )

      rm(focal_result, core_tile)
      gc()
      file.remove(tile_focal_file)
      tmpFiles(current = TRUE, orphan = TRUE, remove = TRUE)

      tile_result_files[[t]] <- tile_crop_file
    }

    message(sprintf(
      "  [%s] [%s] Merging %d tiles for: %s",
      Sys.time(),
      kernel_label,
      n_tiles,
      current_name
    ))
    tile_stack <- sprc(lapply(tile_result_files, rast))
    merge(
      tile_stack,
      filename = merged_file,
      overwrite = TRUE,
      wopt = list(
        gdal = c(
          "COMPRESS=DEFLATE",
          "TILED=YES",
          "BLOCKXSIZE=256",
          "BLOCKYSIZE=256",
          "BIGTIFF=YES"
        )
      )
    )

    rm(tile_stack)
    gc()
    lapply(tile_result_files, file.remove)
    tmpFiles(current = TRUE, orphan = TRUE, remove = TRUE)

    layer_result_files[[i]] <- merged_file
  }

  # stack layers and check

  message(sprintf(
    "[%s] [%s] Stacking all focal layers...",
    Sys.time(),
    kernel_label
  ))
  weighted_sums <- rast(unlist(layer_result_files))
  names(weighted_sums) <- layer_names

  message(sprintf(
    "[%s] [%s] Calculating total focal sum across all habitats...",
    Sys.time(),
    kernel_label
  ))
  original_mask <- rast(paste0(data_dir, "grouped_ecosystems_temp.tif"))[[1]]
  original_mask <- !is.na(original_mask)
  original_mask <- ifel(original_mask == 1, 1, NA)

  weighted_sums_clipped <- terra::mask(weighted_sums, original_mask)

  total_focal_sum <- sum(weighted_sums_clipped, na.rm = TRUE)
  message(sprintf(
    "[%s] [%s] Converting to edge-corrected proportions (0-1)...",
    Sys.time(),
    kernel_label
  ))
  proportions <- ifel(
    total_focal_sum > 0,
    weighted_sums_clipped / total_focal_sum,
    NA
  )
  names(proportions) <- layer_names

  message(sprintf(
    "[%s] [%s] Global min/max per layer (Proportions):",
    Sys.time(),
    kernel_label
  ))
  print(global(proportions, c("min", "max"), na.rm = TRUE))
  print(global(sum(proportions), c("min", "max"), na.rm = TRUE))

  # we moeten de randen van vlaanderen corrigeren aangezien er niets is voorbij de kust op de vlaamse kaart
  # anders zouden we geen strata krijgen op de randen.
  writeRaster(
    proportions,
    paste0(focal_dir, "focal_all_proportions_edgecorrected.tif"),
    overwrite = TRUE,
    wopt = list(
      gdal = c(
        "COMPRESS=DEFLATE",
        "TILED=YES",
        "BLOCKXSIZE=256",
        "BLOCKYSIZE=256",
        "BIGTIFF=YES"
      )
    )
  )

  plot_file <- paste0(focal_dir, "proportions_verification_edgecorrected.png")
  png(plot_file, width = 1920, height = 1080, res = 150)
  plot(
    proportions[[1]],
    main = paste0(
      "Proportional Habitat Density — ",
      names(proportions)[1],
      "  [",
      kernel_label,
      "]"
    ),
    col = map.pal("viridis", 100)
  )
  dev.off()

  message(sprintf(
    "[%s] [%s] Verification plot saved to: %s",
    Sys.time(),
    kernel_label,
    plot_file
  ))

  rm(
    weighted_sums,
    weighted_sums_clipped,
    total_focal_sum,
    proportions,
    original_mask
  )
  gc()
}

message(sprintf("[%s] All runs complete.", Sys.time()))


# ==============================================================================
# Niv2 Gauss 100m classification
# ==============================================================================

message(sprintf(
  "[%s] Starting Categorical Classification (>50%% rule)...",
  Sys.time()
))
data_dir <- "./data/"
proportions <- rast(paste0(
  data_dir,
  "focal_all_proportions_edgecorrected.tif"
))


max_prop_val <- max(proportions, na.rm = TRUE)
dom_class_idx <- which.max(proportions)

n_base_classes <- nlyr(proportions)
mixed_id <- n_base_classes + 1
landbouw_grasland_id <- n_base_classes + 2

landbouw_id <- which(names(proportions) == "landbouw")
grasland_id <- which(names(proportions) == "grasland")

# OPTION1: simplest method
# If the highest proportion is > 0.5, it's the dominanting class otherwise it's mixed
final_class <- ifel(max_prop_val > 0.5, dom_class_idx, mixed_id)
focal_dir <- paste0(data_dir, "Gauss100m")

writeRaster(
  final_class,
  paste0(data_dir, "gauss100m.tif"),
  overwrite = TRUE,
  datatype = "INT1U",
  wopt = list(gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
)

#  OPTION 2
p_landbouw <- proportions[[landbouw_id]]
p_grasland <- proportions[[grasland_id]]
p_comb <- p_landbouw + p_grasland

# criteria:
is_agri_mix <- (p_comb > 0.5) &
  (p_landbouw >= 0.15) &
  (p_grasland >= 0.15) &
  (max_prop_val <= 0.5)

final_class <- ifel(
  max_prop_val > 0.5,
  dom_class_idx,
  ifel(is_agri_mix, landbouw_grasland_id, mixed_id)
)

writeRaster(
  final_class,
  paste0(data_dir, "gauss100m_landbouwgraslandmix.tif"),
  overwrite = TRUE,
  datatype = "INT1U",
  wopt = list(gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
)


# OPTION 3
p_landbouw <- proportions[[landbouw_id]]
p_grasland <- proportions[[grasland_id]]
p_comb <- p_landbouw + p_grasland

# criteria (in total three strata):
# - agriculture  if above 0.5
# - akker if akker is more than 2/3 of p_comb
# - grasland if grasland is more than 2/3 of p_comb
is_agri <- p_comb > 0.5
is_akker <- is_agri & (p_landbouw > (2 / 3) * p_comb)
is_grasland <- is_agri & (p_grasland > (2 / 3) * p_comb)
is_agri_mix <- is_agri & !is_akker & !is_grasland


final_class <- ifel(
  is_akker,
  landbouw_id,
  ifel(
    is_grasland,
    grasland_id,
    ifel(
      is_agri_mix,
      landbouw_grasland_id,
      ifel(max_prop_val > 0.5, dom_class_idx, mixed_id)
    )
  )
)

writeRaster(
  final_class,
  paste0(data_dir, "gauss100m_landbouwgraslandmix_option3.tif"),
  overwrite = TRUE,
  datatype = "INT1U",
  wopt = list(gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
)


# OPTION 4
p_landbouw <- proportions[[landbouw_id]]
p_grasland <- proportions[[grasland_id]]
p_comb <- p_landbouw + p_grasland

# criteria (in total two strata)
# combined above 50%
# assign the akker or landbouw strata to whichever is dominant (no mixed stratum)
is_agri <- (p_comb > 0.5)

final_class <- ifel(
  is_agri,
  ifel(p_landbouw >= p_grasland, landbouw_id, grasland_id),
  ifel(max_prop_val > 0.5, dom_class_idx, mixed_id)
)

writeRaster(
  final_class,
  paste0(data_dir, "gauss100m_landbouwgraslandmix_option4.tif"),
  overwrite = TRUE,
  datatype = "INT1U",
  wopt = list(gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
)


#  OPTION 5 LOGIC
p_landbouw <- proportions[[landbouw_id]]
p_grasland <- proportions[[grasland_id]]
p_comb <- p_landbouw + p_grasland

# criteria (in total three strata):
# - agriculture  if above 0.5
# - akker if akker is more than 1/3
# - grasland if grasland is more than 1/3
is_agri <- p_comb > 0.5
is_akker <- is_agri & (p_landbouw > (1 / 3))
is_grasland <- is_agri & (p_grasland > (1 / 3))
is_agri_mix <- is_agri & !is_akker & !is_grasland


final_class <- ifel(
  is_akker,
  landbouw_id,
  ifel(
    is_grasland,
    grasland_id,
    ifel(
      is_agri_mix,
      landbouw_grasland_id,
      ifel(max_prop_val > 0.5, dom_class_idx, mixed_id)
    )
  )
)

writeRaster(
  final_class,
  paste0(data_dir, "gauss100m_landbouwgraslandmix_option5.tif"),
  overwrite = TRUE,
  datatype = "INT1U",
  wopt = list(gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
)

# Compare the n_cells for each strata
strata <- c(
  "bebouwd",
  "bebouwd_groen",
  "water",
  "akker",
  "grasland",
  "bos",
  "duinen/heide",
  "moeras",
  "estuarium",
  "getijdensgebied/overgangswater",
  "gemengd",
  "akker_grasland_mix"
)

option3 <- freq(rast(paste0(
  data_dir,
  "gauss100m_landbouwgraslandmix_option3.tif"
)))
option5 <- freq(rast(paste0(
  data_dir,
  "gauss100m_landbouwgraslandmix_option5.tif"
)))


option4 <- freq(rast(paste0(
  data_dir,
  "gauss100m_landbouwgraslandmix_option4.tif"
)))
original <- freq(rast(paste0(data_dir, "gauss100m.tif")))


option3$strata <- strata
option4$strata <- strata[1:11]
option5$strata <- strata
original$strata <- strata[1:11]

colnames(option3) <- c("layer", "value", "count_3", "strata")
colnames(option4) <- c("layer", "value", "count_4", "strata")
colnames(option5) <- c("layer", "value", "count_5", "strata")
colnames(original) <- c("layer", "value", "count_original", "strata")

left_join(option3, option4, by = c("strata", "layer", "value")) |>
  left_join(option5, by = c("strata", "layer", "value")) |>
  left_join(original, by = c("strata", "layer", "value"))
