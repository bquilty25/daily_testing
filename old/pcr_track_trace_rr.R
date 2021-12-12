# RELATIVE RISKS

results <- mutate(results, 
                  infectivity_total = infectivity_pre + 
                    infectivity_post)

baseline_low <- data.frame(
  screening             = FALSE,
  first_test_delay      = 0,
  second_test_delay     = NA,
  stringency            = "low"
)

rr_low <- run_rr_analysis(results,
                          main_scenarios, 
                          y_var             = infectivity_total,
                          baseline_scenario = baseline_low,
                          faceting          = index_test_delay ~ stringency,
                          log_scale         = FALSE)

rr_low %>% 
  filter(stringency %in% c("low", "maximum")) %>%
  show_results

rr_low %>% 
  filter(stringency %in% c("high","maximum"), time_in_iso > 8,index_test_delay==2) %>%
  show_results


baseline_max <- data.frame(
  screening = FALSE,
  first_test_delay      = 14,
  second_test_delay     = NA,
  stringency            = "maximum"
)

rr_max <- run_rr_analysis(results,
                          main_scenarios, 
                          y_var = infectivity_post,
                          baseline_scenario = baseline_max,
                          faceting = index_test_delay ~ stringency,
                          log_scale=TRUE)


rr_max %>% 
  filter(stringency %in% c("high"), time_in_iso > 8,index_test_delay==2) %>%
  show_results(reduction = FALSE)
