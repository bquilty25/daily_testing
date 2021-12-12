library(dplyr)
library(purrr)
library(here)
library(tidyr)

#--- strange option in R to remove scientific notation from small floating point
#--- numbers
options(scipen = 999)

source(here::here("R/data_helper_functions.R"))

case_data <- case_data_function()

#--- computing and joining ascertainment_estimates
ascertainment_estimates <- under_reporting_data_function(case_data, 
                                                         here("data/underascertainment_estimates/"))

dates_of_interest <- lubridate::ymd(seq(as.Date("2021-04-01"),
                                        as.Date("2021-04-30"), 1))
  
num_of_days <- dates_of_interest %>% length()
#--- computing incidence estimates
incidence_estimates <- get_adjusted_case_data_national(ascertainment_estimates)

incidence_recent <- incidence_estimates %>% 
  group_by(country, iso_code) %>%
  filter(date %in% dates_of_interest) %>%
  summarise(incidence_mid  = mean(new_cases_adjusted_mid),
            incidence_low  = mean(new_cases_adjusted_low),
            incidence_high = mean(new_cases_adjusted_high))

#--- computing prevalence estimates. Today's date - 8 is the earliest we can estimate for
#--- given reporting delays etc


prevalence_estimates <- dates_of_interest %>% 
  map_df(~global_prevalence_estimates_function(case_data,
                                               ascertainment_estimates,
                                               .)) %>%
  select(iso_code, country, date, everything()) %>%
  drop_na()

prevalence_estimates_mean <- prevalence_estimates %>% 
  group_by(iso_code) %>% 
  summarise(prevalence_mid  = mean(prevalence_mid),
            prevalence_low  = mean(prevalence_low),
            prevalence_high = mean(prevalence_high))

prevalence_estimates_reduced <- prevalence_estimates_mean %>%
  select(iso_code, prevalence_mid, prevalence_low, prevalence_high) %>%
  arrange()

#--- importing clean flight data, cleaned using the clean_flight_data script
travel_data <- read.csv(here("data/flight_data/opensky_april_clean.csv"),
                        row.names=1)
#travel_data_december <- read.csv(here::here("data/flight_data/opensky_december_clean.csv"), row.names=1)
#travel_data_january <- read.csv(here::here("data/flight_data/opensky_january_clean.csv"), row.names=1)

#-------------------- CALCULATING EXPECTED IMPORTS --------------------#
#--- asymptomatic proportions taken from Buitrago-Garcia et al. (2020)
#--- living review, using studies where individuals were followed up only

asymptomatic_prop_mid     <- 0.31
asymptomatic_prop_low     <- 0.26
asymptomatic_prop_high    <- 0.37

#--- number of passengers on each plane, taken from OAG publicly available
#--- figures for 2020
number_of_passegers_per_plane <- 142

imported_cases_pre_sum <- travel_data %>%
  left_join(prevalence_estimates %>% rename(origin_country_iso_code = iso_code)) %>%
  # dplyr::select(origin_country, destination_country, origin_country_iso_code, 
  #               destination_country_iso_code, total_passengers, scaled_travellers, 
  #               prevalence_mid, prevalence_low, prevalence_high) %>%
  group_by(destination_country_iso_code) %>%
  mutate(
    imported_cases_mid   = total_flights*number_of_passegers_per_plane*prevalence_mid,
    imported_cases_low   = total_flights*number_of_passegers_per_plane*prevalence_low,
    imported_cases_high  = total_flights*number_of_passegers_per_plane*prevalence_high)

imported_cases <- imported_cases_pre_sum %>%
  group_by(origin_country_iso_code, destination_country_iso_code) %>%
  summarise_at(.vars = vars(starts_with("imported_cases")),
               .funs = function(x){sum(x, na.rm = T)/num_of_days}) %>%
  filter(origin_country_iso_code == "GBR") %>% 
  left_join(incidence_recent %>% rename(destination_country_iso_code = iso_code))


write.csv(imported_cases, here("outputs/uk_exports_april.csv"))

