source("scripts/utils.r")
pickering %>% pivot_longer.(cols=c(Innova,`SureScreen F`, Encode)) %>% drop_na(value,culture) %>% ggplot(aes(x=factor(culture),fill=factor(value)),position="stack",stat="identity")+geom_bar()+facet_wrap(~name)+labs(x="Culture positive",fill="LFT positive")

pickering %>% pivot_longer.(cols=c(culture,Innova,"SureScreen F", Encode)) %>% ggplot(aes(x=ct,y=value,group=name,colour=name))+geom_point()+geom_smooth(method=glm,method.args=list(family="binomial"))

#Variant characteristics from Kissler et al.
variant_char <- tribble(~"variant", ~"prolif_mean", ~"prolif_lb", ~"prolif_ub", ~"clear_mean", ~"clear_lb", ~"clear_ub", ~"min_ct_mean", ~"min_ct_lb", ~"min_ct_ub",
                        "non_voc", 4.2, 3.3, 5.2, 7.3, 6.1, 8.4, 20.1, 18.3, 21.7,
                        "alpha", 3.4, 2.6, 4.5, 6.2, 5.2, 7.4, 21.0, 19.1, 20.9,
                        "delta", 3.0, 2.2, 4.0, 6.2, 5.2, 7.4, 19.8, 18.0, 22.0,
                        "unvacc", 3.5, 3.0, 4.0, 7.5, 6.8, 8.2, 20.7, 19.8, 20.2,
                        "vacc", 3.2, 2.5, 4.0, 5.5, 4.6, 6.5, 20.5, 19.0, 21.0,
                        "omicron_est",3.0/2, 2.2/2, 4.0/2, 6.2, 5.2, 7.4, 19.8, 18.0, 22.0,
) %>% 
  as_tidytable()

#Make VL trajectories from variant characteristics and asymptomatic fraction
traj <- variant_char %>% 
  filter.(variant=="vacc") %>% 
  group_split.(variant) %>% 
  map.(~make_trajectories(n_cases = 100,n_sims = 100,asymp_parms=asymp_fraction,seed=seed,variant_info=.x)) %>% 
  bind_rows.()

#Calculate daily infectiousness and test positivity, remove never-infectious
traj_ <- traj %>% 
  arrange(variant,sim) %>% 
  group_split.(variant,sim) %>% 
  map.(.f=inf_and_test,sampling_freq = c(1)) %>% 
  bind_rows.() %>% 
  filter.(sum_inf>0)

#create daily testing scenarios 
#delay = how long to wait until resuming testing
#n_negatives = how many consecutive negatives until release
scenarios <- crossing(n_negatives=c(NA,1:3),delay=c(3,5,7)) %>% 
  mutate.(test_to_release=!is.na(n_negatives))

#calculate isolation interval given viral load, test positivity and scenarios
neg_analysis <- traj_ %>%
  select.(-c(m,infectiousness)) %>%
  crossing(scenarios) %>% 
  nest.(data=c(ct, test_t, test_no, test_p, test_label)) %>% #slice_sample.(n=500,.by=c(n_negatives,test_to_release,delay)) %>%
  mutate_rowwise.(iso_interval = earliest_pos_neg(data,test_to_release,n_negatives,delay)) 

qsave(neg_analysis,"neg_analysis.qs")

#cacluate infectious days by days saved
prop_averted <- neg_analysis %>% 
  left_join.(traj_ %>%
               nest(data=c(ct, test_t, test_no, test_p, test_label)) %>% 
               select.(sim,idx,type,variant,m,infectiousness)) %>% 
  mutate.(in_iso=map2.(.f=inf_in_iso,.x=infectiousness,.y=iso_interval)) %>% 
  unnest.(in_iso) %>% 
  filter.(isoFALSE+isoTRUE>0) %>% 
  mutate.(iso_dur=end_iso-start_iso,
          iso_dur=ifelse(iso_dur>10,10,iso_dur),
          days_saved=10-iso_dur) %>% 
  rename.(inf_community=isoFALSE,
          inf_iso=isoTRUE)

qsave(prop_averted,"prop_averted.qs")
  
#create plots
plot_dat <-  prop_averted %>% mutate.(tests_used=ifelse(is.na(n_negatives),0,iso_dur-delay),
                 across.(c(iso_dur,days_saved,inf_community,n_negatives,delay,tests_used),as.factor)) %>% 
  mutate.(n_negatives=fct_explicit_na(n_negatives,"No test"),
          n_negatives=fct_recode(n_negatives,
                                 "3 negatives" = "3",
                                 "2 negatives" = "2",
                                 "1 negative"  = "1"),
          delay=fct_recode(delay,
                           "3 days wait" = "3",
                           "5 days wait" = "5",
                           "7 days wait" = "7"))

plot_a <- plot_dat %>% 
  filter.(variant=="vacc") %>% 
  ggplot(aes(x=days_saved,y=..prop..,group=n_negatives,fill=n_negatives,colour=n_negatives))+
  geom_bar(alpha=0.5,position="dodge")+
  scale_y_continuous("Proportion of infected individuals",position = "left",labels = percent)+
  scale_x_discrete("Days saved vs. 10 day isolation")
  #scale_x_reverse("Days saved vs. 10 day isolation",breaks=breaks_width(-1))

plot_b <- plot_dat %>% 
  filter.(variant=="vacc") %>% 
  ggplot(aes(x=inf_community,y=..prop..,group=n_negatives,fill=n_negatives,colour=n_negatives))+
  geom_bar(alpha=0.5,position="dodge")+
  scale_y_continuous("Proportion of infected individuals",labels = percent)+
  scale_x_discrete("Infectious days in the community")

#tests used post positive test
plot_c <- plot_dat  %>% 
  filter.(variant=="vacc") %>% 
  ggplot(aes(x=tests_used,y=..prop..,group=n_negatives,fill=n_negatives,colour=n_negatives))+
  geom_bar(alpha=0.5,position="dodge")+
  scale_y_continuous("Proportion of infected individuals",labels = percent)+
  scale_x_discrete("Tests used")

plot_a+plot_c+plot_b+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides="collect")&
  theme_minimal()&
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())&
  scale_color_brewer(palette = "Set1")&
  scale_fill_brewer(palette="Set1")&
  ggh4x::facet_nested(delay+n_negatives~.,nest_line = T)

ggsave("main_plot.png",width=210,height=297,units="mm",dpi=600,bg = "white")

#count infectious person hours per sim
prop_averted %>% 
  summarise.(sum_inf_comm=sum(inf_community),
             sum_days_saved=sum(days_saved),
             .by=c(sim,variant,n_negatives,delay,test_to_release)) %>% 
  pivot_longer.(c(sum_inf_comm,sum_days_saved)) %>% 
  group_by(name,variant,n_negatives,delay,test_to_release) %>% 
  summarise_at(vars(value), funs(!!!p_funs)) %>% 
  as_tidytable() %>% 
  ggplot()+
  geom_pointrange(aes(shape=factor(test_to_release),x=interaction(delay,n_negatives),y=`50%`,ymin=`2.5%`,ymax=`97.5%`),position = position_dodge(0.2))+
  facet_wrap(~name,scales = "free")
