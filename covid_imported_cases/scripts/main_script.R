library(dplyr)
library(ggplot2)
library(here)
library(tidyr)
library(data.table)
library(countrycode)
options(scipen = 999)

source(here("R/data_helper_functions.R"))

case_data <- case_data_function()

#--- computing and joining ascertainment_estimates
ascertainment_estimates <- under_reporting_data_function(case_data, here("data/underascertainment_estimates/"))

#--- computing incidence estimates
incidence_estimates <- get_adjusted_case_data_national(ascertainment_estimates)

#--- given that the flight volume data is at a monthly temporal resolution,
#--- we calculate local incidence and prevalence over the same month and average
#--- for the final estimates  
dates_of_interest <- lubridate::ymd(seq(as.Date("2021-04-01"),
                                        as.Date("2021-04-30"), 1))

num_of_days <- dates_of_interest %>% length()

incidence_recent <- incidence_estimates %>% 
  group_by(country, iso_code) %>%
  filter(date %in% dates_of_interest) %>%
  summarise(incidence_mid  = mean(new_cases_adjusted_mid),
            incidence_low  = mean(new_cases_adjusted_low),
            incidence_high = mean(new_cases_adjusted_high))

#--- plotting distribution of all incidence estimates
incidence_recent %>% 
  ggplot(aes(x = incidence_mid)) + 
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 25) + 
  geom_density(alpha = .2, fill="#FF6666") +
  scale_x_continuous(trans = "log10") +
  expand_limits(x = 0.055) + 
  xlab("Incidence") + 
  ylab("Count")

#--- computing prevalence estimates. Today's date - 8 is the earliest we can estimate for
#--- given reporting delays etc
prevalence_estimates <- dates_of_interest %>% 
  map_df(~global_prevalence_estimates_function(case_data,
                                               ascertainment_estimates,
                                               .)) %>%
  select(iso_code, country, date, everything()) %>%
  drop_na() %>% 
  group_by(country, iso_code) %>% 
  summarise(prevalence_mid  = mean(prevalence_mid),
            prevalence_low  = mean(prevalence_low),
            prevalence_high = mean(prevalence_high))


prevalence_estimates_reduced <- prevalence_estimates %>%
  select(iso_code, country, prevalence_mid, prevalence_low, prevalence_high) %>%
  arrange()

#--- plotting distribution of all prevalence estimates
prevalence_estimates %>% 
  ggplot(aes(x = prevalence_mid)) + 
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 25) + 
  geom_density(alpha = .2, fill="#FF6666") +
  scale_x_continuous(trans = "log10", labels = scales::percent) +
  expand_limits(x = 0.055) + 
  xlab("Prevalence (%)") + 
  ylab("Count")

#--- importing clean flight data, cleaned using the clean_flight_data script
travel_data <- read.csv(here("data/flight_data/opensky_april_clean.csv"), row.names=1)

#-------------------- CALCULATING EXPECTED IMPORTS --------------------#
#--- asymptomatic proportions taken from Buitrago-Garcia et al. (2020)
#--- living review, using studies where individuals were followed up only

asymptomatic_prop_mid  <- 0.31
asymptomatic_prop_low  <- 0.26
asymptomatic_prop_high <- 0.37

#--- number of passengers on each plane, taken from OAG publicly available
#--- figures for 2020
number_of_passegers_per_plane <- 142

imported_cases_pre_sum <- travel_data %>%
  left_join(prevalence_estimates %>% rename(origin_country_iso_code = iso_code)) %>%
  group_by(destination_country_iso_code) %>%
  mutate(
    imported_cases_mid  = total_flights*number_of_passegers_per_plane*prevalence_mid,
    imported_cases_low  = total_flights*number_of_passegers_per_plane*prevalence_low,
    imported_cases_high = total_flights*number_of_passegers_per_plane*prevalence_high)

imported_cases <- imported_cases_pre_sum %>%
  dplyr::summarise_at(.vars = vars(starts_with("imported_cases")),
                      .funs = function(x){sum(x, na.rm = T)/num_of_days}) %>%
  dplyr::arrange(destination_country_iso_code)


# calculating the incidence in each destination country in May
incidence_recent_asymptomatic <-  incidence_recent %>%
  group_by(iso_code) %>%
  mutate(
    incidence_total_mid  = incidence_mid/(1 - asymptomatic_prop_mid),
    incidence_total_low  = incidence_low/(1 - asymptomatic_prop_low),
    incidence_total_high = incidence_high/(1 - asymptomatic_prop_high)) %>%
  mutate(destination_country_iso_code = iso_code) %>%
  select(destination_country_iso_code,
         incidence_total_mid, 
         incidence_total_low, 
         incidence_total_high)

#--- putting all the estimate together
all_estimates_together <- incidence_recent_asymptomatic %>%
  select(-destination_country_iso_code) %>%
  left_join(prevalence_estimates_reduced) %>%
  left_join(imported_cases %>% rename(iso_code = destination_country_iso_code)) %>%
  mutate(country = countrycode(iso_code, "iso3c", "iso.name.en")) %>% 
  select(iso_code, country, everything())

write.csv(all_estimates_together, here("outputs/april_estimates.csv"))

