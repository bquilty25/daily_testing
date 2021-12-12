

make_plots <- function(
  x, 
  input,
  main_scenarios = NULL,
  log_scale = FALSE,
  #fixed = TRUE,
  text_size = 2.5,
  #trav_vol_manual = NULL,
  xlab = "Days in quarantine\n(including 1 day delay on testing results)",
  sum = FALSE,
  y_var = "days_released_inf",
  faceting = NULL){
  
  #browser()
  
  ylabA = "Number of infectious persons\nreleased per index case"
  
  if (sum){
    ylabB = 
      "Total number of person-days infectiousness\nremaining for released secondary case"
  } else {
    ylabB = 
      "Number of days infectiousness\nremaining per released secondary case"
  }
  
  
  all_grouping_vars <- all.vars(faceting)
  
  # ... should be country, type, pre_board_screening, etc.
  x_summaries <- x %>%
    filter(stage_released == "Infectious") %>%
    make_released_quantiles(., all_grouping_vars) %>%
    inner_join(input)
  
  
  # why are the plo medians coming out backwars in Low?
  figA_data <- plot_data(input, 
                         x_summaries,
                         main_scenarios)
  
  
  # deal with the join and stage_released # need to have scenario guide the labelling
  
  # browser()
  figA <- figA_data %>% 
    #filter(pre_board_screening == "None") %>% 
    # # this filter should be done outside the function
    make_release_figure(
      x_summaries = .,
      input     = input,
      text_size = text_size,
      xlab      = xlab,
      ylab      = ylabA, 
      faceting  = faceting) 
  
  ## person-days
  
  x_days_summaries <- 
    make_released_time_quantiles(x, 
                                 y_var = y_var,
                                 all_grouping_vars,
                                 sum = sum)
  
  
  
  figB_data <- plot_data(input = input, 
                         x_summaries = 
                           x_days_summaries,
                         main_scenarios)
  
  figB <- figB_data %>% 
    #filter(pre_board_screening == "None") %>%  
    make_release_figure(
      x = .,
      input=input,
      xlab = "Days in quarantine\n(including 1 day delay on testing results)",
      text_size = text_size,
      ylab = ylabB,
      faceting = faceting) 
  
  
  fig <- figA + figB + plot_layout(ncol = 1, guide = "collect") +
    plot_annotation(tag_levels = "A",
                    theme = theme(legend.position = "bottom"))
  
  return(fig)
  
}

