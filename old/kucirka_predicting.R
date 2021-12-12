


max_time <- incubation_times %>%
  mutate(tmax = exp_to_onset + onset_to_recov) %>%
  summarise(tmax = max(tmax)) %>% pull(tmax)


max_time <- max_time + max(input$first_test_delay + 
                             replace_na(input$quar_dur, 0))
# predict from day 1 forward
dat_pred <- expand.grid(day = seq(1, ceiling(max_time))) %>%
  dplyr::bind_cols({predict(dat_gam, ., se.fit = T)} %>% as.data.frame) %>%
  dplyr::mutate(L = fit - 1.96*se.fit,
                U = fit + 1.96*se.fit) %>%
  dplyr::mutate_at(.vars = dplyr::vars(fit, L, U),
                   .funs = boot::inv.logit)

# return P(t)
pcr_curve <- dplyr::select(dat_pred, day, pcr_p = fit)
