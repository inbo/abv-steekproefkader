library(rvest)
library(tidyverse)
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
  ) |>
  as.data.frame()

write.csv(all_hashes, "./data/hashes_list.csv")
