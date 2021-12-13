source("scripts/utils.r")


traj <- make_trajectories(n_cases = 100,n_sims = 100,variant = "wild")

traj_ <- traj %>% 
  arrange(sim) %>% 
  group_split.(sim) %>% 
  map.(.f=inf_and_test,sampling_freq = c(1)) %>% 
  bind_rows.()

neg_analysis <- traj_ %>% 
  nest(ct, test_t, test_no, test_p, test_label) %>% slice_sample.(n=1000) %>% 
  crossing(n_negatives=1:3) %>% 
  mutate.(earliest_positive = map2.(.f = earliest_pos_neg, .x = data,.y=n_negatives)) %>%
  unnest_wider.(earliest_positive)

neg_analysis %>% 
  ggplot()+
  geom_bar(aes(x=10-(end_iso-start_iso),y=..prop..))+
  scale_x_continuous("Duration of isolation",breaks=scales::breaks_width(1))+
  facet_wrap(~n_negatives,ncol = 1)
