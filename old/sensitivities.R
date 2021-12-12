# look at all sensitivities

results_baseline <-
  input %>%
  filter(
    delay_scaling != 0.5, 
    waning        != "waning_canada_total") %>%
  pull(scenario) %>%
  get(results_name)[.] %>%
  make_days_plots(.,
                  faceting = index_test_delay ~ stringency,
                  y_labels = grep(value = T, pattern = "prior",
                                  x = infectivity_labels, invert = T),
                  dir = results_name,
                  base = "baseline",
                  sum = F)

results_scaling <-
  input %>%
  filter(waning != "waning_canada_total") %>%
  pull(scenario) %>%
  get(results_name)[.] %>%
  make_days_plots(.,
                  faceting = index_test_delay + delay_scaling ~ stringency,
                  y_labels = grep(value = T, pattern = "prior",
                                  x = infectivity_labels, invert = T),
                  dir = results_name,
                  base = "scaling",
                  sum = F)

results_waning <-
  input %>%
  filter(delay_scaling != 0.5) %>%
  pull(scenario) %>%
  get(results_name)[.] %>%
  make_days_plots(.,
                  faceting = index_test_delay + waning ~ stringency,
                  y_labels = grep(value = T, pattern = "prior",
                                  x = infectivity_labels, invert = T),
                  dir = results_name,
                  base = "waning",
                  sum = F)
