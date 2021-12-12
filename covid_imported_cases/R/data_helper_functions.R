get_adjusted_case_data_national <- function(under_reporting_data_arg)
{
  
  asymptomatic_mid <- 0.5
  asymptomatic_low <- 0.1
  asymptomatic_high <- 0.7
  
  under_reporting_raw_data <- under_reporting_data_arg
  
  under_reporting_and_case_data <- case_data %>%
    dplyr::left_join(under_reporting_raw_data) %>%
    dplyr::group_by(country) %>%
    dplyr::arrange(country, date) %>%
    tidyr::drop_na() %>%
    dplyr::select(date, iso_code, country, new_cases, new_deaths, estimate, lower, upper)
  
  
  data_out <- under_reporting_and_case_data %>%
    dplyr::group_by(country) %>%
    dplyr::mutate(new_cases_smoothed             = zoo::rollmean(new_cases, k = 7, fill = NA),
                  new_cases_adjusted_mid         = new_cases/estimate,
                  new_cases_adjusted_low         = new_cases/upper,
                  new_cases_adjusted_high        = new_cases/lower,
                  new_cases_adjusted_smooth_mid  = zoo::rollmean(new_cases_adjusted_mid, k = 7, fill = NA),
                  new_cases_adjusted_smooth_low  = zoo::rollmean(new_cases_adjusted_low, k = 7, fill = NA),
                  new_cases_adjusted_smooth_high = zoo::rollmean(new_cases_adjusted_high, k = 7, fill = NA)) %>%
    # dplyr::mutate(cumulative_incidence_mid  = cumsum(new_cases_adjusted_mid)/(popData2019*(1 - asymptomatic_mid)),
    #               cumulative_incidence_low  = cumsum(new_cases_adjusted_low)/(popData2019*(1 - asymptomatic_low)),
    #               cumulative_incidence_high = cumsum(new_cases_adjusted_high)/(popData2019*(1 - asymptomatic_high))) %>%
    dplyr::mutate(date_infection  = date - 9) %>%
    # dplyr::mutate(cumulative_incidence_mid = dplyr::case_when(cumulative_incidence_mid >= 1 ~ 1,
    #                                                           cumulative_incidence_mid <= 0 ~ 0,
    #                                                           cumulative_incidence_mid > 0 & cumulative_incidence_mid < 1 ~ cumulative_incidence_mid)) %>%
    # dplyr::mutate(cumulative_incidence_low = dplyr::case_when(cumulative_incidence_low > 1 ~ 1,
    #                                                           cumulative_incidence_low < 0 ~ 0,
    #                                                           cumulative_incidence_low > 0 & cumulative_incidence_low < 1 ~ cumulative_incidence_low)) %>%
    # dplyr::mutate(cumulative_incidence_high = dplyr::case_when(cumulative_incidence_high > 1 ~ 1,
    #                                                            cumulative_incidence_high < 0 ~ 0,
    #                                                            cumulative_incidence_high > 0 & cumulative_incidence_high < 1 ~ cumulative_incidence_high)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(country = stringr::str_replace_all(country, "_", " ")) %>%
    dplyr::mutate(country = dplyr::case_when(country == "United States of America" ~ "USA",
                                             country != "United States of America" ~ country)) %>%
    dplyr::mutate(country = dplyr::case_when(country == "United Kingdom" ~ "UK",
                                             country != "United Kingdom" ~ country))
  
  
  return(data_out)
  
}


case_data_function <- function()
{
  
  case_data <- readr::read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv") %>%
    dplyr::mutate(date = lubridate::ymd(date)) %>%
    dplyr::select(iso_code, country = location, date, new_cases, new_deaths, population)
  
  return(case_data)
}

under_reporting_data_function <- function(case_data_arg, ascertainment_estimates_path)
{
  
  data_path <- here::here(ascertainment_estimates_path)
  files <- dir(path = data_path,
               pattern = "*.rds")
  
  under_reporting_data_raw <- dplyr::tibble(iso_code = files) %>% 
    dplyr::mutate(file_contents = purrr::map(iso_code, 
                                             ~ readRDS(file.path(data_path, .)))
                  
    ) %>% 
    tidyr::unnest(cols = c(file_contents)) %>%
    dplyr::mutate(iso_code = stringr::str_remove(iso_code, "result_")) %>% 
    dplyr::mutate(iso_code = stringr::str_remove(iso_code, ".rds")) %>%
    dplyr::group_by(iso_code) %>%
    dplyr::select(date, everything()) %>%
    dplyr::left_join(case_data_arg) %>%
    dplyr::select(date, country, iso_code, everything()) %>%
    dplyr::group_by(iso_code) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(country = stringr::str_replace_all(country, "_", " "))
  
  return(under_reporting_data_raw)
}

global_prevalence_estimates_function <- function(case_data, under_reporting_data_arg, date_arg, recent = TRUE)
{
  
  date_arg <- as.Date(date_arg)
  if(date_arg - Sys.Date() > 13)
  {
    date_arg_new <- date_arg
  }
  else(date_arg - Sys.Date() <= 13)
  {
    date_arg_new <- Sys.Date() - 13
    
  }
  
  #--- small block of code that imports and cleans populations if needed
  #--- population also included in OWID dataset so not required at the moment
  # library(wpp2019)
  # data(pop)
  # population_raw_data <- pop
  # 
  # populations <- population_raw_data %>%
  #   dplyr::mutate(iso_code = countrycode::countrycode(country_code, "iso3n", "iso3c")) %>% 
  #   tidyr::drop_na() %>% 
  #   dplyr::select(iso_code, population = "2020") %>%
  #   dplyr::mutate(population = population*1000) %>%
  #   dplyr::add_row(iso_code = "AND", population = 77323) %>%
  #   dplyr::arrange(iso_code)
  
  #--- asymptomatic proportions taken from Buitrago-Garcia et al. (2020)
  #--- living review, using studies where individuals were followed up only
  asymptomatic_estimate_mid <- 0.31
  asymptomatic_estimate_low <- 0.26
  asymptomatic_estimate_high <- 0.37
  
  new_case_estimates_recent <- case_data %>%
    dplyr::filter(date_arg - 9 < date) %>% 
    dplyr::mutate(country = stringr::str_replace_all(country, "_", " ")) %>%
    dplyr::mutate(iso_code = countrycode::countrycode(country, "country.name", "iso3c", custom_match = c('Kosovo' = 'RKS'))) %>%
    dplyr::group_by(country, iso_code, population) %>%
    dplyr::summarise(total_new_cases = sum(new_cases))
  
  # turn off scientific notation
  options(scipen=999)
  
  if(recent == TRUE)
  {
    
    under_reporting_estimates_at_date <- under_reporting_data_arg %>%
      dplyr::group_by(iso_code) %>%
      dplyr::filter(date == max(date))
  } else
  {
    under_reporting_estimates_at_date <- under_reporting_data_arg %>%
      dplyr::group_by(iso_code) %>%
      dplyr::filter(date == date_arg - 13)
  }
  
  most_recent_estimates_together <-  new_case_estimates_recent %>% 
    dplyr::left_join(under_reporting_estimates_at_date) %>%
    dplyr::group_by(iso_code) %>%
    dplyr::mutate(prevalence_mid  = (total_new_cases/(estimate*(1 - asymptomatic_estimate_mid))/population),
                  prevalence_low  = (total_new_cases/(upper*(1 - asymptomatic_estimate_low))/population),
                  prevalence_high = (total_new_cases/(lower*(1 - asymptomatic_estimate_high))/population)) %>%
    dplyr::select(iso_code, country, total_new_cases, estimate, lower, upper, population, prevalence_mid, prevalence_low, prevalence_high) %>%
    dplyr::arrange(country) %>%
    dplyr::mutate(date = date_arg) %>% 
    tidyr::drop_na()
  
  return(most_recent_estimates_together)
  
}

oxford_restrictions_fun <- function(month_arg)
{
  #--- downloading and cleaning the oxford travel restriction ratings
  oxford_url <- "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"
  oxford_data <- oxford_url %>% readr::read_csv()
  
  #--- keeping Aprils data, cleaning and keeping just maximum rating in April
  oxford_travel_neat <- oxford_data %>%
    dplyr::select(CountryCode, Date, restriction_rating = "C8_International travel controls")  %>%
    dplyr::mutate(date = lubridate::ymd(Date)) %>%
    dplyr::select(-Date) %>%
    dplyr::filter(lubridate::month(date) == month_arg) %>%
    dplyr::group_by(CountryCode) %>%
    dplyr::rename(iso_code = CountryCode) %>%
    dplyr::group_by(iso_code) %>%    
    dplyr::summarise(travel_restrictions_rating = max(restriction_rating)) %>%
    dplyr::mutate(country = countrycode::countrycode(iso_code, "iso3c", "country.name", custom_match = c("RKS" = "Kosovo"))) %>%
    dplyr::select(country, iso_code, travel_restrictions_rating) %>%
    dplyr::arrange(country)
  
  return(oxford_travel_neat)
}

#--- function which reads in 2020 data and scales it down (for rating 4 countries)
#--- using the ratio of 2020 to 2019 flight volumes
oag_2020_april_data_function <- function()
{
  
  #--- reading in 2019 and 2020 OAG data
  oag_2019_data_raw <- readr::read_csv(here::here("data", "raw_data/oag_flight_data_2019.csv"))
  oag_2020_data_raw <- readr::read_csv(here::here("data", "raw_data/flight_data_all_2020_Feb_Apr_final_monthly.csv"))
  
  #--- cleaning the 2019 data
  oag_2019_data_neat <- oag_2019_data_raw %>% 
    dplyr::select(iso_code_dep = dep_country_code,
                  iso_code_arr = arr_country_code,
                  total_travellers = bookings,
                  date = month_rev)
  
  #--- filtering for just April 2019, moving from 2 letter to 3 letter codes and cleaning
  oag_2019_april <- oag_2019_data_neat %>%
    dplyr::filter(date == "2019-04") %>%
    dplyr::mutate(iso_code_dep = countrycode::countrycode(iso_code_dep, "iso2c", "iso3c"),
                  iso_code_arr = countrycode::countrycode(iso_code_arr, "iso2c", "iso3c")) %>%
    dplyr::arrange(iso_code_dep, iso_code_arr)
  
  #--- filtering for just April 2019, moving from 2 letter to 3 letter codes and cleaning
  oag_april_2020_neat <- oag_2020_data_raw %>% 
    dplyr::select(iso_code_dep = dep_country_code,
                  iso_code_arr = arr_country_code, 
                  total_travellers = flight_volume_total, 
                  date = year_month_rev) %>%
    dplyr::mutate(iso_code_dep = countrycode::countrycode(iso_code_dep, "iso2c", "iso3c"),
                  iso_code_arr = countrycode::countrycode(iso_code_arr, "iso2c", "iso3c")) %>%
    dplyr::filter(date == "2020-04") %>%
    dplyr::select(-date) %>%
    dplyr::arrange(iso_code_dep, iso_code_arr)
  
  #--- downloading the oxford travel restriction ratings
  oxford_travel_neat <- oxford_restrictions_fun(month = 4)
  
  #--- 2019 April data for countries with level 3 or lower restriction rating
  #--- to scale down the level 4 countries
  countries_level_3_or_lower <- oxford_travel_neat %>% 
    dplyr::filter(travel_restrictions_rating != 4) %>%
    dplyr::pull(iso_code)
  
  #--- calculating average scaling factor from 2019 to 2020
  oag_april_2019_tmp <- oag_2019_april %>%
    dplyr::filter(iso_code_dep %in% countries_level_3_or_lower) %>%
    dplyr::rename(traveller_volume_2019 = total_travellers) %>%
    dplyr::select(-date) %>%
    dplyr::filter(iso_code_dep != iso_code_arr) %>%
    dplyr::group_by(iso_code_dep, iso_code_arr) %>%
    dplyr::summarise(traveller_volume_2019 = sum(traveller_volume_2019))
  
  oag_april_2020_tmp <- oag_april_2020_neat %>%
    dplyr::filter(iso_code_dep %in% countries_level_3_or_lower) %>%
    dplyr::rename(traveller_volume_2020 = total_travellers) %>%
    dplyr::filter(iso_code_dep != iso_code_arr) %>%
    dplyr::group_by(iso_code_dep, iso_code_arr) %>%
    dplyr::summarise(traveller_volume_2020 = sum(traveller_volume_2020))
  
  scaling_factor_2019_to_2020 <- oag_april_2019_tmp %>%
    dplyr::left_join(oag_april_2020_tmp) %>%
    dplyr::filter(traveller_volume_2019 != 0) %>%
    dplyr::mutate(scaling_factor = traveller_volume_2020/traveller_volume_2019) %>%
    dplyr::ungroup() %>%
    tidyr::drop_na() %>%
    dplyr::summarise(scaling_factor_mean = mean(scaling_factor)) %>%
    dplyr::pull(scaling_factor_mean)
  
  
  scaling_factor_2019_to_2020 <- oag_april_2019_tmp %>%
    dplyr::left_join(oag_april_2020_tmp) %>%
    dplyr::filter(traveller_volume_2019 != 0) %>%
    dplyr::mutate(scaling_factor = traveller_volume_2020/traveller_volume_2019) %>%
    dplyr::ungroup() %>%
    tidyr::drop_na() %>%
    dplyr::summarise(scaling_factor_mean = mean(scaling_factor)) %>%
    dplyr::pull(scaling_factor_mean)
  
  #--- 2020 data set to be used in analysis
  #--- after scaling down for reduction between 2019 and 2020 
  #--- and cleaning removing intra-country travel
  
  oag_april_2020_scaled <- oag_april_2020_neat %>%
    dplyr::group_by(iso_code_dep, iso_code_arr) %>%
    dplyr::filter(iso_code_dep != iso_code_arr) %>%
    dplyr::group_by(iso_code_dep, iso_code_arr) %>%
    dplyr::summarise(total_travellers = sum(total_travellers)) %>%
    dplyr::mutate(scaled_travellers = dplyr::case_when(!(iso_code_dep %in% countries_level_3_or_lower) ~ total_travellers*scaling_factor_2019_to_2020,
                                                       iso_code_dep %in% countries_level_3_or_lower  ~ total_travellers)) %>%
    dplyr::mutate(origin_country = countrycode::countrycode(iso_code_dep, "iso3c", 'country.name', 
                                                            custom_match = c('RKS' = 'Kosovo')),
                  destination_country = countrycode::countrycode(iso_code_arr, "iso3c", 'country.name', 
                                                                 custom_match = c('RKS' = 'Kosovo')))
  
  return(oag_april_2020_scaled)
  
}

clean_flight_data <- function(open_sky_data)
{
  
  flights_by_origin <- open_sky_data %>% 
    dplyr::mutate(date = lubridate::ymd(day))%>%
    dplyr::select(date, four_letter_code = origin) %>%
    dplyr::left_join(airport_lookup) %>%
    #tidyr::drop_na() %>%
    dplyr::rename(origin_country = country,
                  origin_code = four_letter_code) %>%
    dplyr::select(-three_letter_code)
  
  flights_by_destination <- open_sky_data %>% 
    dplyr::mutate(date = lubridate::ymd(day))%>%
    dplyr::select(date, four_letter_code = destination) %>%
    dplyr::left_join(airport_lookup) %>%
    #tidyr::drop_na() %>%
    dplyr::rename(destination_country = country,
                  destination_code = four_letter_code) %>%
    dplyr::select(-three_letter_code)
  
  
  flights_by_origin$destination_code <- flights_by_destination$destination_code
  flights_by_origin$destination_country <- flights_by_destination$destination_country
  
  flight_data_out <- flights_by_origin %>%
    tidyr::drop_na() %>%
    dplyr::select(date, origin_code, destination_code, origin_country, destination_country) %>%
    dplyr::rowwise() %>%
    dplyr::filter(origin_country != destination_country)
  
  return(flight_data_out)
  
}

calculate_exports <- function(month_arg){
  
  file_path <- paste0("outputs/uk_exports_", month_arg, ".csv")
  
  uk_exports <- read.csv(file_path) %>%
    rename(iso_a3 = destination_country_iso_code) %>%
    select(iso_a3, 
           exports_mid = imported_cases_mid, 
           exports_low = imported_cases_low, 
           exports_high = imported_cases_high, 
           incidence_mid, 
           incidence_low,
           incidence_high) %>%
    mutate(month = month_arg) %>%
    add_row(iso_a3 = "GBR", 
            exports_mid = NA,
            exports_low = NA, 
            exports_high = NA, 
            month = month_arg) %>%
    arrange(iso_a3)
  
  return(uk_exports)
}

