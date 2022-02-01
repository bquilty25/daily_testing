source("scripts/utils.r")
pickering %>% pivot_longer.(cols=c(Innova,`SureScreen F`, Encode)) %>% drop_na(value,culture) %>% ggplot(aes(x=factor(culture),fill=factor(value)),position="stack",stat="density")+geom_bar()+facet_wrap(~name)+labs(x="Culture positive",fill="LFT positive")

pickering %>% pivot_longer.(cols=c(culture,Innova,"SureScreen F", Encode)) %>% ggplot(aes(x=vl,y=value,group=name,colour=name))+geom_point()+geom_smooth(method=glm,method.args=list(family="binomial"))

#Variant characteristics from Kissler et al.
variant_char <- tribble(~"variant", ~"prolif_mean", ~"prolif_lb", ~"prolif_ub", ~"clear_mean", ~"clear_lb", ~"clear_ub", ~"max_vl_mean",~"max_vl_lb",~"max_vl_ub",
                        "non_voc", 4.2, 3.3, 5.2, 7.3, 6.1, 8.4, 8.2, 7.7, 11.6,
                        "alpha", 3.4, 2.6, 4.5, 6.2, 5.2, 7.4, 7.9, 8.0, 11.5,
                        "delta", 3.0, 2.2, 4.0, 6.2, 5.2, 7.4, 8.3, 7.7, 11.6,
                        "unvacc", 3.5, 3.0, 4.0, 7.5, 6.8, 8.2, 8.0, 8.2, 11.5,
                        "vacc", 3.2, 2.5, 4.0, 5.5, 4.6, 6.5, 8.1, 7.9, 11.5,
                        "omicron",4.52, 3.61, 5.54, 5.35, 4.78, 6.00, 7.28, 7.02, 7.54,
) %>% 
  as_tidytable() 

#Make VL trajectories from variant characteristics and asymptomatic fraction
traj <- variant_char %>% 
  filter.(variant%in%c("unvacc","vacc","delta","omicron")) %>% 
  group_split.(variant) %>% 
  map.(~make_trajectories(n_sims = 1000,asymp_parms=asymp_fraction,seed=seed,variant_info=.x)) %>% 
  bind_rows.()

#' TODO precalculate test / infection outcomes

#Calculate daily infectiousness and test positivity, remove never-infectious
traj_ <- traj %>% 
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(m = m,start=start,end=end)))  %>% 
  unnest.(infectiousness) %>% 
  mutate.(culture_p        = stats::predict(culture_mod, type = "response", newdata = tidytable(vl=vl)),
          infectious_label = rbernoulli(n=n(),p=culture_p),
          test_p           = stats::predict(innova_mod, type = "response", newdata = tidytable(vl=vl)),
          test_label       = rbernoulli(n=n(),p = test_p)) %>% 
  filter.(any(infectious_label==T),.by=c(sim,variant))

#create daily testing scenarios 
#delay = how long to wait until resuming testing
#n_negatives = how many consecutive negatives until release
scenarios <- crossing(n_negatives=c(NA,1:2),delay=c(3,4,5,6,7)) %>% 
  mutate.(test_to_release=!is.na(n_negatives),scenario_id=row_number.())

#calculate isolation interval given viral load, test positivity and scenarios
neg_analysis_dat <- scenarios %>% 
  crossing(traj_) %>% 
  replace_na.(list(test_label=FALSE))

neg_analysis_dat[,
  start_iso := fifelse(any(test_label), t[which.max(test_label)], Inf),
  by=c("sim","variant","scenario_id")
]

neg_analysis_dat[(t >= start_iso + delay),
   testndelays := c(1, diff(t)),
   by=c("sim","variant","scenario_id")
]

neg_analysis_dat[
  testndelays == 1,
   end_iso := earliest_pos_neg(.SD),
   .SDcols = c("t","test_label","n_negatives","test_to_release","delay","start_iso"),
   by=c("sim","variant","scenario_id")
]


#cacluate infectious days by days saved
prop_averted <- neg_analysis_dat %>% 
  filter.(t>=start_iso) %>% 
  fill.(testndelays,end_iso,.direction="up",.by=c(sim,variant,scenario_id)) %>% 
  mutate.(iso=between(t,start_iso,end_iso-1),.by=c(sim,variant,scenario_id)) %>% 
  summarise.(inf_iso=sum(infectious_label,na.rm = T),.by=c(sim,variant,scenario_id,n_negatives,delay,iso,start_iso,end_iso)) %>% 
  pivot_wider.(names_from = iso,values_from = inf_iso,names_prefix = "iso_") %>% 
  replace_na.(list(iso_FALSE=0,iso_TRUE=0)) %>% 
  mutate.(inf_days=iso_TRUE+iso_FALSE,
          iso_dur=end_iso-start_iso,
          iso_dur_over_10=iso_dur>10,
          iso_dur=ifelse(iso_dur>10,10,iso_dur),
          days_saved=10-iso_dur,
          tests_used=ifelse(is.na(n_negatives),0,iso_dur-delay)) %>% 
  rename.(inf_community=iso_FALSE,
          inf_iso=iso_TRUE) %>% 
  mutate.(vacc_status=ifelse(str_detect(variant,"unvacc"),"Unvaccinated","Vaccinated"),
          omicron=ifelse(str_detect(variant,"omicron"),"Omicron (assump.)","Non-Omicron"))

qsave(prop_averted,"prop_averted.qs")
prop_averted <- qread("prop_averted.qs")

#create plots
plot_dat <-  prop_averted %>% 
  #filter.(variant=="unvacc") %>% 
  mutate.(across.(c(iso_dur,n_negatives,delay,tests_used,variant),as.factor)) #%>% 
  # mutate.(n_negatives=fct_explicit_na(n_negatives,"No test"),
  #         n_negatives=fct_recode(n_negatives,
  #                                "3 negatives" = "3",
  #                                "2 negatives" = "2",
  #                                "1 negative"  = "1"),
  #         delay=fct_recode(delay,
  #                          "3 days wait" = "3",
  #                          "5 days wait" = "5",
  #                          "7 days wait" = "7"))

plot_a <- plot_dat %>% 
  ggplot(aes(x=days_saved,y=..prop..,group=n_negatives,fill=n_negatives,colour=n_negatives))+
  geom_bar(alpha=0.7)+
  scale_y_continuous("Proportion of infected individuals",position = "left")+
  scale_x_continuous("Days saved vs. 10 day isolation",breaks = breaks_width(1))

plot_b <- plot_dat %>% 
  ggplot(aes(x=factor(inf_community),y=..prop..,group=n_negatives,fill=n_negatives,colour=n_negatives))+
  geom_bar(alpha=0.7)+
  labs(x="Infectious days in the community")

#tests used post positive test
plot_c <- plot_dat  %>% 
  ggplot(aes(x=tests_used,y=..prop..,group=n_negatives,fill=n_negatives,colour=n_negatives))+
  geom_bar(alpha=0.7)+
  labs(x="Tests used")

plot_a+plot_c+plot_b+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides="collect")&
  theme_minimal()&
  scale_color_brewer(palette = "Set2")&
  scale_fill_brewer(palette="Set2")&
  ggh4x::facet_nested(delay+n_negatives~variant,nest_line = T)&
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

ggsave("output/main_plot.png",width=210,height=297,units="mm",dpi=600,bg = "white")

#days saved
plot_2a_dat <- prop_averted %>% 
  mutate.(across.(c(iso_dur,n_negatives,delay,tests_used),as.factor)) %>% 
  mutate.(n_negatives=fct_explicit_na(n_negatives,"No test"),
          n_negatives=fct_recode(n_negatives,
                                 "3 days negative" = "3",
                                 "2 days negative" = "2",
                                 "1 day negative"  = "1"),
          n_negatives=fct_relevel(n_negatives,"No test"),
          delay=fct_recode(delay,
                           "3 days wait" = "3",
                           "5 days wait" = "5",
                           "7 days wait" = "7")) %>% 
  pivot_longer.(c(days_saved)) %>% 
  group_by(name,variant,scenario_id,delay,n_negatives,vacc_status,omicron) %>% 
  summarise_each(funs(mean, median, sd, se=sd(.)/sqrt(n()), lower=bayestestR::ci(.,method="ETI")$CI_low,upper=bayestestR::ci(.,method="ETI")$CI_high), value) %>% 
  as_tibble()

(plot_2a <- plot_2a_dat %>%  
    filter(variant%in%c("omicron")) %>% 
    ggplot()+
    coord_flip()+
    geom_pointrange(aes(x=fct_rev(delay),
                        group= fct_rev(n_negatives),
                        colour=n_negatives,
                        y=median,
                        ymin=lower,
                        ymax=upper),
                    position = position_dodge(0.5),
                    fatten = 1)+
    scale_y_continuous(breaks=breaks_width(1),limits=c(0,NA))+
    labs(x="",y="Days saved vs. 10 days isolation\nper individual",
         colour="Number of consecutive days of negative tests required for release")
)

plot_2b_dat <- prop_averted %>% 
  mutate.(across.(c(iso_dur,n_negatives,delay,tests_used),as.factor)) %>% 
  mutate.(n_negatives=fct_explicit_na(n_negatives,"No test"),
          n_negatives=fct_recode(n_negatives,
                                 "3 days negative" = "3",
                                 "2 days negative" = "2",
                                 "1 day negative"  = "1"),
          n_negatives=fct_relevel(n_negatives,"No test"),
          delay=fct_recode(delay,
                           "3 days wait" = "3",
                           "5 days wait" = "5",
                           "7 days wait" = "7")) %>% 
  #filter.(inf_community>0) %>% 
  as_tibble() %>% 
  group_by(variant,n_negatives,delay,vacc_status,omicron) %>% 
  mutate(x=bootf(inf_community)) %>% 
  unnest(x)

(plot_2b <- plot_2b_dat %>% 
    filter(variant%in%c("omicron")) %>% 
  ggplot()+
  coord_flip()+
  geom_pointrange(aes(x=fct_rev(delay),
                      group= fct_rev(n_negatives),
                      colour=n_negatives,
                      y=Mean,
                      ymin=Lower,
                      ymax=Upper
                      ),
                  position = position_dodge(0.5),
                  fatten = 1
                  )+
    scale_y_continuous(labels=scales::percent_format(1L))+
    #scale_y_continuous(trans="pseudo_log",breaks=c(0,10,100,1000,10000,1e5))+
  labs(x="",y="Days infectious in the community\nper 10,000 infected individuals",
       colour="Number of consecutive days of negative tests required for release")
  )

plot_2c_dat <- prop_averted %>% 
  mutate.(across.(c(iso_dur,n_negatives,delay),as.factor)) %>% 
  mutate.(n_negatives=fct_explicit_na(n_negatives,"No test"),
          n_negatives=fct_recode(n_negatives,
                                 "3 days negative" = "3",
                                 "2 days negative" = "2",
                                 "1 day negative"  = "1"),
          n_negatives=fct_relevel(n_negatives,"No test"),
          delay=fct_recode(delay,
                           "3 days wait" = "3",
                           "5 days wait" = "5",
                           "7 days wait" = "7")) %>% 
  summarise.(tests_used=sum(tests_used)*(10000/n()),
             n=n(),
             .by=c(sim,variant,n_negatives,delay,test_to_release,vacc_status,omicron)) %>% 
  pivot_longer.(c(tests_used)) %>% 
  group_by(name,variant,n_negatives,delay,test_to_release,vacc_status,omicron) %>% 
  summarise_each(funs(mean, median, sd, se=sd(.)/sqrt(n()), lower=bayestestR::ci(.,method="ETI")$CI_low,upper=bayestestR::ci(.,method="ETI")$CI_high), value) %>% 
  pivot_wider(values_from = quantile,names_from=prob) %>% 
  as_tidytable()

(plot_2c <- plot_2c_dat %>% 
    filter.(variant%in%c("omicron")) %>% 
    ggplot()+
    coord_flip()+
    geom_pointrange(aes(x=fct_rev(delay),
                        group= fct_rev(n_negatives),
                        colour=n_negatives,
                        y=median,
                        ymin=lower,
                        ymax=upper
    ),
    position = position_dodge(0.5),
    fatten = 1
    )+
    scale_y_continuous()+
    labs(x="",y="Tests used \nper 10,000 infected individuals",
         colour="Number of consecutive days of negative tests required for release")
  
)

plot_2a+plot_2b+plot_2c+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect")&
  #scale_color_brewer(palette="Set2")&
  scale_color_manual(values=c("#22577A",
    "#38A3A5",
    "#57CC99",
    "#80ED99"),guide=guide_legend(title.position = "top"))&
  theme_minimal()&
  theme(panel.grid.major.y = element_blank(),
              axis.ticks = element_line(),
              axis.line.x.bottom = element_line(),
              axis.line.y.left = element_line(),
              legend.position = "bottom")

ggsave("output/main_plot2_vacc.png",width=250,height=100,units="mm",dpi=600,bg = "white")

