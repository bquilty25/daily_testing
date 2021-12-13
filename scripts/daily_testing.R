source("scripts/utils.r")

variant_char <- tribble(~"variant", ~"prolif_mean", ~"prolif_lb", ~"prolif_ub", ~"clear_mean", ~"clear_lb", ~"clear_ub", ~"min_ct_mean", ~"min_ct_lb", ~"min_ct_ub",
                        "non_voc", 4.2, 3.3, 5.2, 7.3, 6.1, 8.4, 20.1, 18.3, 21.7,
                        "alpha", 3.4, 2.6, 4.5, 6.2, 5.2, 7.4, 21.0, 19.1, 20.9,
                        "delta", 3.0, 2.2, 4.0, 6.2, 5.2, 7.4, 19.8, 18.0, 22.0,
                        "unvacc", 3.5, 3.0, 4.0, 7.5, 6.8, 8.2, 20.7, 19.8, 20.2,
                        "vacc", 3.2, 2.5, 4.0, 5.5, 4.6, 6.5, 20.5, 19.0, 21.0
) %>% as_tidytable()

traj <- variant_char %>% 
  group_split.(variant) %>% 
  map.(~make_trajectories(n_cases = 100,n_sims = 100,asymp_parms=asymp_fraction,seed=seed,variant_info=.x)) %>% 
  bind_rows.()

traj_ <- traj %>% 
  arrange(variant,sim) %>% 
  group_split.(variant,sim) %>% 
  map.(.f=inf_and_test,sampling_freq = c(1)) %>% 
  bind_rows.()

neg_analysis <- traj_ %>% select.(-c(m,infectiousness)) %>% 
  nest.(data=c(ct, test_t, test_no, test_p, test_label)) %>% slice_sample.(n=1000) %>% 
  crossing(n_negatives=1:3) %>% 
  mutate.(earliest_positive = map2.(.f = earliest_pos_neg, .x = data,.y=n_negatives)) %>%
  unnest_wider.(earliest_positive)

neg_analysis %>% 
  ggplot()+
  geom_bar(aes(x=10-(end_iso-start_iso),y=..prop..))+
  scale_x_continuous("Days saved",breaks=scales::breaks_width(1))+
  facet_grid(variant~n_negatives,
             labeller=labeller(n_negatives=function(x){paste0(x, " day(s) of negative tests required for release")}))

#calculate auc between the two dates and compare to overall = infectiousness averted

prop_averted <- neg_analysis %>% 
  unnest.(infectiousness) %>% 
  mutate.(in_iso=between.(t,start_iso,end_iso)) %>% 
  summarise.(iso_inf=sum(culture),.by=c(variant,sim,n_negatives,in_iso)) %>% 
  pivot_wider.(values_from=iso_inf,names_from=in_iso) %>% 
  mutate.(prop_averted=(`TRUE`/(`FALSE`+`TRUE`))) %>% 
  nest_by.(c(variant,n_negatives)) %>% 
  mutate.(q = map(.x = data, ~quantile(.x$prop_averted,probs=c(0.5,0.025,0.975)))) %>% 
  unnest_wider(q) 


prop_averted %>% 
  ggplot()+
  geom_pointrange(aes(x=n_negatives,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
  scale_y_continuous(limits=c(0,1))
