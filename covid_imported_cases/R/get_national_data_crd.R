get_national_data_crd <- function()
{
  library(data.table)
  library(covidregionaldata)
  library(countrycode)
  crd_data_dt <- get_national_data() %>% 
    data.table() %>% 
    .[, c("date", "country", "iso_code", "cases_new", "deaths_new")] %>%
    .[order(country, date)] %>% 
    .[!is.na(country)] %>% 
    .[!is.na(cases_new)] %>% 
    .[!is.na(deaths_new)] %>% 
    .[, iso_code := countrycode(iso_code, "iso2c", "iso3c", 
                                custom_match = c("XK" = "XXK"))] %>% # for Kosovo
    setnames(., 
             old = c("iso_code", "cases_new", "deaths_new"),
             new = c("iso3c", "new_cases", "new_deaths")) 
  
  return(crd_data_dt)
}

tmp <- get_national_data_crd()

p_tmp_cases <- tmp %>% 
  ggplot(aes(x = date)) + 
  geom_line(aes(y = new_cases), colour = "blue") + 
  facet_wrap(~country, scales = "free")

p_tmp_deaths <- tmp %>% 
  ggplot(aes(x = date)) + 
  geom_line(aes(y = new_deaths), colour = "red") + 
  facet_wrap(~country, scales = "free")

tmp[, total_deaths := sum(new_deaths), by = "country"]
tmp[, .SD[total_deaths < 10 & date == max(date)], by = "country"]

ggsave("all_countries_death_data.png", p_tmp_deaths, width = 20, height = 20)

countrycode::countrycode(tmp[, unique(iso3c)], "iso3c", "iso.name.en")
