library(magrittr)
library(here)
library(data.table)
#library(dplyr)

#--- importing and cleaning OpenSky flight numbers data
source(here::here("R/data_helper_functions.R"))

opensky_raw_data_path <- "data/flight_data/flightlist_20210401_20210430.csv.gz"

#opensky_raw_data <- readr::read_csv(here::here(opensky_raw_data_path))
opensky_raw_data <- fread(here(opensky_raw_data_path))

airport_lookup <- read.delim(here("data/airports.dat"), sep = ",") %>%
  data.table() %>%
  .[, .SD, .SDcols = c("Country", "three_letter_code", "four_letter_code")] %>% 
  setnames(., old = "Country", new = "country")

opensky_cleaned_data  <- clean_flight_data(opensky_raw_data)

opensky_total_flights <- opensky_cleaned_data %>%
  dplyr::filter(origin_country != destination_country) %>%
  dplyr::group_by(origin_country, destination_country) %>%
  dplyr::summarise(total_flights = dplyr::n()) %>%
  dplyr::mutate(origin_country_iso_code = countrycode::countrycode(origin_country, "country.name", "iso3c"),
                destination_country_iso_code = countrycode::countrycode(destination_country, "country.name", "iso3c")) %>%
  dplyr::ungroup() %>%
  dplyr::select(origin_country_iso_code, destination_country_iso_code, total_flights)

#--- its quite slow, so I save the output as a .csv as import it in
#--- the main script

write.csv(opensky_total_flights, here("data/flight_data/opensky_april_clean.csv"))
