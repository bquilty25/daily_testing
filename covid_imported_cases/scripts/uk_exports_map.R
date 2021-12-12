library(ggplot2)
library(sf)
library(rnaturalearth)
library(tidyr)
library(dplyr)
library(data.table)
library(patchwork)
library(here)
options(scipen = 999)

source(here("R/data_helper_functions.R"))

# reading in s-gene dropout proportion timeseries in UK
s_gene_dropout <- fread("data/sgtfv_clean.csv", drop = 1) 

# taking three time points from the timeseries for plotting later
all_dates <- lubridate::ymd(seq(as.Date("2020-10-01"),
                                as.Date("2021-04-30"), 1)) %>% 
  data.table(date = .)

t1_dates <- all_dates[date <= "2020-10-31"]
t2_dates <- all_dates["2020-11-01" <= date & date <= "2020-11-30"]
t3_dates <- all_dates["2020-12-01" <= date & date <= "2020-12-31"]
t4_dates <- all_dates["2021-01-01" <= date & date <= "2021-01-31"]
t5_dates <- all_dates["2021-02-01" <= date & date <= "2021-02-28"]
t6_dates <- all_dates["2021-03-01" <= date & date <= "2021-03-31"]
t6_dates <- all_dates["2021-04-01" <= date & date <= "2021-04-30"]

s_gene_prop_t1 <- s_gene_dropout[date %in% t1_dates[, date]] %>%
  .[region == "national"] %>%
  .[, mean(sgtfv)]

s_gene_prop_t2 <- s_gene_dropout[date %in% t2_dates[, date]] %>%
  .[region == "national"] %>%
  .[, mean(sgtfv)]

s_gene_prop_t3 <- s_gene_dropout[date %in% t3_dates[, date]] %>%
  .[region == "national"] %>%
  .[, mean(sgtfv)]

s_gene_prop_t4 <- s_gene_dropout[date %in% t4_dates[, date]] %>%
  .[region == "national"] %>%
  .[, mean(sgtfv)]

#--- we assume (in the absence of easily available data) that the proportion of
#--- the new variant nationally has been fixed at 1 since the end of janurary,
#--- given that the data suggests it was very close to one before then and
#--- accelerating quickly
s_gene_prop_t5 <- 1
s_gene_prop_t6 <- 1
s_gene_prop_t7 <- 1

# calculating number of exported cases from the UK at the same three time points
uk_exports_october  <- calculate_exports("october")
uk_exports_november <- calculate_exports("november")
uk_exports_december <- calculate_exports("december")
uk_exports_january  <- calculate_exports("january")
uk_exports_february <- calculate_exports("february")
uk_exports_march    <- calculate_exports("march")
uk_exports_april    <- calculate_exports("april")

#--- calculating the number of "non-variant" exports
uk_exports_october_non_variant <- uk_exports_october %>%
  mutate(exports_mid  = exports_mid*(1 - s_gene_prop_t1), 
         exports_low  = exports_low*(1 - s_gene_prop_t1), 
         exports_high = exports_high*(1 - s_gene_prop_t1), 
         variant = FALSE)

uk_exports_november_non_variant <- uk_exports_november %>%
  mutate(exports_mid  = exports_mid*(1 - s_gene_prop_t2), 
         exports_low  = exports_low*(1 - s_gene_prop_t2), 
         exports_high = exports_high*(1 - s_gene_prop_t2), 
         variant = FALSE) 

uk_exports_december_non_variant <- uk_exports_december %>%
  mutate(exports_mid  = exports_mid*(1 - s_gene_prop_t3),
         exports_low  = exports_low*(1 - s_gene_prop_t3),
         exports_high = exports_high*(1 - s_gene_prop_t3),
         variant = FALSE)

uk_exports_january_non_variant <- uk_exports_january %>%
  mutate(exports_mid  = exports_mid*(1 - s_gene_prop_t4), 
         exports_low  = exports_low*(1 - s_gene_prop_t4), 
         exports_high = exports_high*(1 - s_gene_prop_t4), 
         variant = FALSE) 

uk_exports_february_non_variant <- uk_exports_february %>%
  mutate(exports_mid  = exports_mid*(1 - s_gene_prop_t5),
         exports_low  = exports_low*(1 - s_gene_prop_t5),
         exports_high = exports_high*(1 - s_gene_prop_t5),
         variant = FALSE) 

uk_exports_march_non_variant <- uk_exports_march %>%
  mutate(exports_mid = exports_mid*(1 - s_gene_prop_t6), 
         exports_low = exports_low*(1 - s_gene_prop_t6), 
         exports_high = exports_high*(1 - s_gene_prop_t6), 
         variant = FALSE) 

uk_exports_april_non_variant <- uk_exports_april %>%
  mutate(exports_mid = exports_mid*(1 - s_gene_prop_t7), 
         exports_low = exports_low*(1 - s_gene_prop_t7), 
         exports_high = exports_high*(1 - s_gene_prop_t7), 
         variant = FALSE) 

#--- calculating the number of "variant" exports
uk_exports_october_variant <- uk_exports_october %>% 
  mutate(exports_mid  = exports_mid*s_gene_prop_t1, 
         exports_low  = exports_low*s_gene_prop_t1, 
         exports_high = exports_high*s_gene_prop_t1, 
         variant = TRUE) 

uk_exports_november_variant <- uk_exports_november %>% 
  mutate(exports_mid  = exports_mid*s_gene_prop_t2, 
         exports_low  = exports_low*s_gene_prop_t2, 
         exports_high = exports_high*s_gene_prop_t2, 
         variant = TRUE) 

uk_exports_december_variant <- uk_exports_december %>% 
  mutate(exports_mid  = exports_mid*s_gene_prop_t3,
         exports_low  = exports_low*s_gene_prop_t3,
         exports_high = exports_high*s_gene_prop_t3,
         variant = TRUE)

uk_exports_january_variant <- uk_exports_january %>% 
  mutate(exports_mid  = exports_mid*s_gene_prop_t4,
         exports_low  = exports_low*s_gene_prop_t4,
         exports_high = exports_high*s_gene_prop_t4,
         variant = TRUE)

uk_exports_february_variant <- uk_exports_february %>% 
  mutate(exports_mid  = exports_mid*s_gene_prop_t5, 
         exports_low  = exports_low*s_gene_prop_t5, 
         exports_high = exports_high*s_gene_prop_t5, 
         variant = TRUE)

uk_exports_march_variant <- uk_exports_march %>% 
  mutate(exports_mid  = exports_mid*s_gene_prop_t6, 
         exports_low  = exports_low*s_gene_prop_t6, 
         exports_high = exports_high*s_gene_prop_t6, 
         variant = TRUE)

uk_exports_april_variant <- uk_exports_april %>% 
  mutate(exports_mid  = exports_mid*s_gene_prop_t7, 
         exports_low  = exports_low*s_gene_prop_t7, 
         exports_high = exports_high*s_gene_prop_t7, 
         variant = TRUE)


uk_exports <- uk_exports_april_non_variant %>% 
  bind_rows(uk_exports_march_non_variant) %>%
  bind_rows(uk_exports_february_non_variant) %>%
  bind_rows(uk_exports_january_non_variant) %>%
  bind_rows(uk_exports_december_non_variant) %>%
  bind_rows(uk_exports_november_non_variant) %>%
  bind_rows(uk_exports_october_non_variant) %>%
  bind_rows(uk_exports_october_variant) %>%
  bind_rows(uk_exports_november_variant) %>%
  bind_rows(uk_exports_december_variant) %>%
  bind_rows(uk_exports_january_variant) %>% 
  bind_rows(uk_exports_february_variant) %>% 
  bind_rows(uk_exports_march_variant) %>%
  bind_rows(uk_exports_april_variant) %>%
  mutate(month = factor(month, 
                        levels = c("october",
                                   "november",
                                   "december",
                                   "january",
                                   "february",
                                   "march", 
                                   "april"))) %>% 
  data.table() %>% 
  .[order(month, variant, iso_a3)]

write.csv(uk_exports, 
          "outputs/uk_exports_all.csv",
          row.names = FALSE)
