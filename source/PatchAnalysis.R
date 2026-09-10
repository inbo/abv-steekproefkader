# script for patch analysis
# to see how many patches of swamp there are so that not all locations are in the same patch.

library(terra)
library(sf)

map <- rast("data/gauss100m_landbouwgraslandmix_option4.tif")
strata <- 1:11
names(strata) <- c(
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
  "gemengd"
)
freq(map)
small_strata <- c("duinen/heide", "moeras")
patch_freq_list <- vector(mode = "list", length = length(small_strata))
names(patch_freq_list) <- names(small_strata)
for (stratum in strata[small_strata]) {
  print(stratum)
  stratum_raster <- map == stratum
  patches_raster <- patches(stratum_raster, directions = 8, zeroAsNA = TRUE)
  patch_freq_list[[stratum]] <- mutate(
    freq(patches_raster),
    stratum = names(strata)[stratum]
  )[, c(2, 3, 4)]
  colnames(patch_freq_list[[stratum]]) <- c("id", "count", "stratum")
}
patch_freq_df <- do.call(rbind, patch_freq_list)
patch_freq_df |>
  group_by(stratum) |>
  summarise(
    n_patches = n(),
    mean_size = mean(count)
  )

# save for use in the steekproefgrootte.qmd:
saveRDS(patch_freq_df, "data/patch_freq_df.rds")
# compare these values then to the needed sample size.
