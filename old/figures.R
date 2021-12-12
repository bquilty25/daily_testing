# main figures

faceting      <- index_test_delay + delay_scaling + waning ~ stringency
faceting_wide <- index_test_delay + delay_scaling ~ waning + stringency

results_infectivity <- 
  get(results_name) %>%
  make_days_plots(.,
      faceting = faceting,
      y_labels = grep(value = T, pattern = "prior",
                      x = infectivity_labels, invert = T),
      dir = results_name,
      base = "all",# all
      sum = F)

infectivity_labels %>%
  map2(.x = ., .y = names(.),
       .f = ~set_names(.x, .y)) %>%
  map(
    ~make_days_plots(get(results_name),
                     input, 
                     faceting = faceting_wide,
                     y_labels = .x,
                     dir = results_name,
                     base = paste(results_name, names(.x),sep="_"),
                     sum = F))


