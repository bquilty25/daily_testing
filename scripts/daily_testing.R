source("scripts/utils.r")
pickering %>% pivot_longer.(cols=c(Innova,`SureScreen F`, Encode)) %>% drop_na(value,culture) %>% ggplot(aes(x=factor(culture),fill=factor(value)),position="stack",stat="identity")+geom_bar()+facet_wrap(~name)+labs(x="Culture positive",fill="LFT positive")

pickering %>% pivot_longer.(cols=c(culture,Innova,"SureScreen F", Encode)) %>% ggplot(aes(x=ct,y=value,group=name,colour=name))+geom_point()+geom_smooth(method=glm,method.args=list(family="binomial"))

variant_char <- tribble(~"variant", ~"prolif_mean", ~"prolif_lb", ~"prolif_ub", ~"clear_mean", ~"clear_lb", ~"clear_ub", ~"min_ct_mean", ~"min_ct_lb", ~"min_ct_ub",
                        "non_voc", 4.2, 3.3, 5.2, 7.3, 6.1, 8.4, 20.1, 18.3, 21.7,
                        "alpha", 3.4, 2.6, 4.5, 6.2, 5.2, 7.4, 21.0, 19.1, 20.9,
                        "delta", 3.0, 2.2, 4.0, 6.2, 5.2, 7.4, 19.8, 18.0, 22.0,
                        "unvacc", 3.5, 3.0, 4.0, 7.5, 6.8, 8.2, 20.7, 19.8, 20.2,
                        "vacc", 3.2, 2.5, 4.0, 5.5, 4.6, 6.5, 20.5, 19.0, 21.0,
                        "omicron_est",3.0/2, 2.2/2, 4.0/2, 6.2, 5.2, 7.4, 19.8, 18.0, 22.0,
) %>% 
  as_tidytable()

traj <- variant_char %>% 
  group_split.(variant) %>% 
  map.(~make_trajectories(n_cases = 100,n_sims = 100,asymp_parms=asymp_fraction,seed=seed,variant_info=.x)) %>% 
  bind_rows.()

traj_ <- traj %>% 
  arrange(variant,sim) %>% 
  group_split.(variant,sim) %>% 
  map.(.f=inf_and_test,sampling_freq = c(1)) %>% 
  bind_rows.()

neg_analysis <- traj_ %>%
  select.(-c(m,infectiousness)) %>%
  crossing(n_negatives=1:3,delay=c(0,2,4)) %>% 
  nest.(data=c(ct, test_t, test_no, test_p, test_label)) %>% slice_sample.(n=1000,.by=variant) %>% 
  mutate_rowwise.(iso_interval = earliest_pos_neg(data,n_negatives,delay)) %>% 
  separate.(iso_interval,into = c("start_iso","end_iso")) %>% 
  mutate.(across.(c(start_iso,end_iso),as.integer))

neg_analysis %>% 
  replace_na.(list(end_iso=10)) %>% 
  ggplot()+
  geom_boxplot(aes(y=end_iso-start_iso,x=factor(n_negatives)))+
  labs(x="Day(s) of negative tests required for release",
       y="Duration of isolation")+
  #scale_x_continuous("Duration of isolation",breaks=scales::breaks_width(1))+
  facet_grid(delay~variant,
             labeller=labeller(n_negatives=function(x){paste0(x, " day(s) of negative tests required for release")},
                               delay=function(x){paste0(x, " day(s) from first positive to next test")})
             )

#calculate auc between the two dates and compare to overall = infectiousness averted

prop_averted <- neg_analysis %>% left_join.(traj_ %>% select.(sim,idx,variant,m,infectiousness)) %>% 
  unnest.(infectiousness) %>% 
  mutate.(in_iso=between.(t,start_iso,end_iso),
          iso_dur=end_iso-start_iso) %>% 
  summarise.(iso_inf=sum(culture),.by=c(variant,sim,idx,n_negatives,in_iso,delay,iso_dur)) %>% 
  pivot_wider.(values_from=iso_inf,names_from=in_iso) %>% 
  mutate.(prop_averted=(`TRUE`/(`FALSE`+`TRUE`))) %>% 
  nest_by.(c(variant,n_negatives,delay)) %>% 
  mutate.(q = map.(.x = data, ~quantile(.x$prop_averted,probs=c(0.5,0.025,0.975),na.rm=T))) %>% 
  unnest_wider.(q) 


prop_averted %>% ggplot(aes(x=delay,y=n_negatives,fill=..50.))+
  geom_tile()+
  geom_text(aes(label=paste0(label_percent(accuracy = 0.1)(..50.),"\n(",label_percent(accuracy = 0.1)(..2.5.),",\n", label_percent(accuracy = 0.1)(..97.5.),")")),
            size=3)+
  facet_wrap(~variant)+
  labs(x="Days after positive to resume testing",y="Days of negatives required for release",fill="Infectiousness averted")+
  scale_fill_viridis_c()+
  coord_fixed(2)
  #geom_smooth()
ggsave("infectiousness_averted.png",dpi=600)

prop_averted %>% 
  ggplot()+
  geom_pointrange(aes(x=as.factor(n_negatives),y=..50.,ymin=..2.5.,ymax=..97.5.))+
  #scale_y_continuous(limits=c(0,1))+
  facet_grid(delay~variant)
