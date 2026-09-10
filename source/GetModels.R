library(rvest)
library(dplyr)
library(purrr)
library(stringr)
library(n2kanalysis)
library(abvanalysis)
source("./source/AnalysisFunctions.R")
project <- "abv"

con <- connect_inbo_dbase("S0008_00_Meetnetten")

query <- "
select lp.name as Hok
, l.name as Telpunt
, l.geom.STX as Lon
, l.geom.STY as Lat
from [staging_Meetnetten].[projects_projectlocation] pl
inner join [staging_Meetnetten].[locations_location] l on l.id = pl.location_id
inner join [staging_Meetnetten].[locations_location] lp on lp.id = l.parent_id
inner join [staging_Meetnetten].[projects_project] p on p.id = pl.project_id
where 1=1
and p.name = 'Algemene Broedvogelmonitoring (ABV)'
"

abv_coordinates <- DBI::dbGetQuery(con, query)

head(abv_coordinates)

colnames(abv_coordinates) <- c("square", "point", "Lon", "Lat")

# merge the data from the models above with the lcoation data

#claude code to scrape all the hashes
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
  res <- variability_comparison(model)
  ests <- get_scenario_ests(model)

  if (inherits(res, "error")) {
    message("--> Error encountered for ", soort_name, ": ", res$message)
    failed_soorten <- c(failed_soorten, soort_name)
    vars_list[[soort_name]] <- NA
  } else {
    vars_list[[soort_name]] <- c(res, ests)
  }
}
print(failed_soorten)

vars_df <- data.frame(do.call(rbind, vars_list))
vars_df <- vars_df |>
  mutate(
    species = rownames(vars_df),
    sigma2 = var_square + var_point,
    sigma = sqrt(var_square + var_point)
  )


saveRDS(
  vars_df,
  "./source/Nieuwe map/PowerAnalyse/Data/VogelModelParameters.rds"
)

# quick visual check of the betas and lambdas per stratum

betas <- vars_df[, c(
  "species",
  "Bos",
  "Heide.en.duin",
  "Landbouw",
  "Moeras.en.water",
  "Suburbaan",
  "Urbaan"
)]
betas |>
  pivot_longer(
    cols = c(
      "Bos",
      "Heide.en.duin",
      "Landbouw",
      "Moeras.en.water",
      "Suburbaan",
      "Urbaan"
    ),
    names_to = "stratum",
    values_to = "betas"
  ) |>
  mutate(lambdas = exp(betas)) |>
  ggplot(aes(x = betas, col = stratum)) +
  stat_ecdf() +
  xlim(c(-2.8, 2.8))
