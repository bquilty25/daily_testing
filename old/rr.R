run_rr_analysis <- function(
  released_times,
  main_scenarios,
  baseline_scenario,
  text_size = 2.5,
  log_scale=TRUE,
  y_var = infectivity_post,
  faceting = ~ stringency){
  set.seed(145)
  
  #browser()
  #Parameters
  
  y_var <- enquo(y_var)
  
  baseline <- inner_join(baseline_scenario, input )
  
  stringencies <- distinct(released_times, stringency, scenario)
  
  
  released_times_summaries <- 
    mutate(released_times, 
           time_in_iso = released_t - traced_t) %>% 
    mutate(infectivity=!!y_var) %>% 
    dplyr::select(scenario,stringency,sim,idx,released_test,infectivity, 
                  index_test_delay, -contains("delay"),screening)
  
  
  baseline_summaries <- 
    inner_join(released_times_summaries,
               baseline) %>% 
    rename("baseline_infectivity"   = infectivity,
           #"baseline_released_test" = released_test,
           "baseline_scenario"      = scenario,
           "baseline_stringency"    = stringency) %>%
    dplyr::select(sim,idx, contains("baseline"),
                  -screening,-pathogen,-contains("delay") ,index_test_delay)
  
  
  n_risk_ratios <- released_times_summaries %>% 
    filter(grepl(x = released_test,
                 pattern ="Released after"),
           !grepl(x = released_test,
                  pattern = "\\+")) %>%
    inner_join(baseline_summaries) %>%   
    mutate(ratio=(infectivity)/(baseline_infectivity)) %>% 
    replace_na(list(ratio=1)) %>% 
    nest(data = -c(scenario,baseline_scenario,index_test_delay)) %>%
    mutate(Q = map(.x = data, ~quantile(.x$ratio, probs = probs)),
           M = map_dbl(.x = data, ~mean(.x$ratio))) %>%
    unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    inner_join(stringencies) %>%
    inner_join(input)
  
  ylabA <- sprintf("Rate ratio of infectivity in comparison to\n%s stringency, %i day quarantine, %s scenario",
                   baseline_scenario$stringency,
                   with(baseline_scenario,
                        first_test_delay + screening + 
                          ifelse(is.na(second_test_delay), 0, second_test_delay)),"no testing")
  
  rr_fig_data <-
    plot_data(input, n_risk_ratios, main_scenarios = NULL) 
  
  rr_fig <- rr_fig_data %>%
    make_release_figure(
      x         = .,
      input     = input,
      text_size = text_size,
      xlab      = xlab,
      ylab      = ylabA, 
      log_scale = log_scale,
      hline     = "dashed",
      faceting  = faceting)  
  
  
  file <- paste(names(baseline_scenario),
                baseline_scenario, sep = "-", 
                collapse = " ")
  
  list("png", "pdf") %>%
    map(~ggsave(filename = paste0("results/rr_figs_baseline_",
                                  file,".",.x),
                plot=rr_fig,
                width = 260, 
                height = 80*nrow(distinct(ungroup(n_risk_ratios),
                                          !!lhs(faceting))), units="mm",
                dpi = 320,
                device = ifelse(.x=="pdf",cairo_pdf,
                                "png")))
  
  
  return(rr_fig_data)
}
