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
  filter.(variant%in%c("omicron")) %>% 
  group_split.(variant) %>% 
  map.(~make_trajectories(n_sims = 1000,asymp_parms=asymp_fraction,seed=seed,variant_info=.x)) %>% 
  bind_rows.()

#' TODO precalculate test / infection outcomes

#Calculate daily infectiousness and test positivity, remove never-infectious
traj_ <- traj %>%
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = ceiling(max(end)+5)
  )))  %>%
  unnest.(infectiousness) %>%
  crossing(
    lower_inf_thresh = c(TRUE, FALSE),
    higher_detect_thresh = c(TRUE, FALSE)
  ) %>%
  mutate.(
    culture_p        = stats::predict(
      object = inf_model_choice(lower_inf_thresh),
      type = "response",
      newdata = tidytable(vl = vl)
    ),
    infectious_label = rbernoulli(n = n(),
                                  p = culture_p),
    test_p           = stats::predict(
      object =  test_model_choice(higher_detect_thresh),
      type = "response",
      newdata = tidytable(vl = vl)
    ),
    test_label       = rbernoulli(n = n(),
                                  p = test_p),
    .by = c(lower_inf_thresh, higher_detect_thresh)
  ) %>%
  replace_na.(list(test_label       = FALSE,
                   infectious_label = FALSE))

#create daily testing scenarios 
#delay = how long to wait until resuming testing
#n_negatives = how many consecutive negatives until release
scenarios <- crossing(n_negatives = c(0:2),
                      delay       = c(3, 4, 5, 6, 7)) %>%
             mutate.(scenario_id = row_number.())

#calculate isolation interval given viral load, test positivity and scenarios
neg_analysis_dat <- traj_ %>% crossing.(scenarios)

neg_analysis_dat[,
  start_iso := fifelse(any(test_label), t[which.max(test_label)], sample(c(2,3),size = 1)), #if never test positive, assume day 3 (symp onset)
  by=c("sim","variant","scenario_id", "lower_inf_thresh", "higher_detect_thresh")
]

neg_analysis_dat[(t >= start_iso + delay),
   testndelays := c(1, diff(t)),
   by=c("sim","variant","scenario_id", "lower_inf_thresh", "higher_detect_thresh")
]

neg_analysis_dat[
  testndelays == 1,
   end_iso := earliest_pos_neg(.SD),
   .SDcols = c("t","test_label","n_negatives","delay","start_iso"),
   by=c("sim","variant","scenario_id", "lower_inf_thresh", "higher_detect_thresh")
]


#cacluate infectious days by days saved
prop_averted <- neg_analysis_dat %>%
  #filter.(t >= start_iso) %>%
  fill.(
    testndelays,
    end_iso,
    .direction = "up",
    .by = c(sim, variant, scenario_id, lower_inf_thresh, higher_detect_thresh)
  ) %>%
  mutate.(iso = between(t, start_iso, end_iso - 1),
          .by = c(sim, variant, scenario_id, lower_inf_thresh, higher_detect_thresh)) %>%
  summarise.(
    inf_iso = sum(infectious_label, na.rm = T),
    .by = c(
      sim,
      variant,
      scenario_id,
      lower_inf_thresh, 
      higher_detect_thresh,
      n_negatives,
      delay,
      iso,
      start_iso,
      end_iso
    )
  ) %>%
  pivot_wider.(names_from = iso,
               values_from = inf_iso,
               names_prefix = "iso_") %>%
  #replace_na.(list(iso_FALSE = 0, iso_TRUE = 0)) %>%
  mutate.(
    inf_days = iso_TRUE + iso_FALSE,
    iso_dur = end_iso - start_iso,
    days_saved = 10 - iso_dur,
    tests_used = ifelse(n_negatives==0,0,iso_dur - delay+1)
  ) %>%
  rename.(inf_community = iso_FALSE,
          inf_iso = iso_TRUE) %>%
  mutate.(across.(c(n_negatives, delay), as.factor),
          n_negatives=fct_rev(n_negatives)) #%>% 
  #filter.(higher_detect_thresh==F)

qsave(prop_averted,"prop_averted.qs")
prop_averted <- qread("prop_averted.qs")

#days saved
plot_2a_dat <- prop_averted[,
                            #{x=summarise_func(days_saved);as.list(x)},
                            {x=Hmisc::smean.cl.boot(days_saved,B = 1000);as.list(x)},
                            by=c("variant","n_negatives","delay","lower_inf_thresh", "higher_detect_thresh")]

(plot_2a <- plot_2a_dat %>%  
    ggplot(aes(x=delay,
               group= n_negatives,
               colour=n_negatives,
               fill=n_negatives,
               y=Mean))+
    coord_flip()+
    geom_col(position = position_dodge(0.9))+
    geom_errorbar(aes(x=delay,
                      group= n_negatives,
                      ymin=Lower,
                      ymax=Upper
    ),
    colour="black",
    position = position_dodge(0.9),
    width=.2)+
    scale_y_continuous(breaks=breaks_width(1),limits=c(0,NA))+
    labs(x="",y="Days saved vs. 10 days isolation\nper individual",
         colour="Number of consecutive days of negative tests required for release")
)

plot_2b_dat <- prop_averted[,
                            {x=Hmisc::smean.cl.boot(inf_community>0,B = 1000);as.list(x)},
                            #{x=boot_func(inf_community,stat="prop_greater");as.list(x)},
                            #{x=Hmisc::smean.cl.boot(inf_community>0,B = 1000);as.list(x)},
                            by=c("variant","n_negatives","delay","lower_inf_thresh", "higher_detect_thresh")]

(plot_2b <- plot_2b_dat %>% 
 ggplot(aes(x=delay,
               group= n_negatives,
               colour=n_negatives,
               fill=n_negatives,
               y=Mean*10000))+
    coord_flip()+
    geom_col(position = position_dodge(1))+
    geom_errorbar(aes(x=delay,
                      group= n_negatives,
                      ymin=Lower*10000,
                      ymax=Upper*10000
    ),
    colour="black",
    position = position_dodge(1),
    width=.2)+
    scale_y_continuous(trans = pseudo_log_trans(1, 10),
                       #breaks = c(0,lseq(1,1000,length.out = 4)),
                       labels=scales::percent_format(0.1),
                       #breaks=c(0,0.005,0.05,0.5)
                       )+
  labs(x="",
       y="Proportion infectious after release")
  )


plot_2c_dat <- prop_averted[,
               {x=Hmisc::smean.cl.boot(tests_used,B = 1000);as.list(x)},
               by=c("variant","n_negatives","delay","lower_inf_thresh", "higher_detect_thresh")]

(plot_2c <- plot_2c_dat %>% 
    ggplot(aes(x=delay,
           group= n_negatives,
           colour=n_negatives,
           fill=n_negatives,
           y=Mean*10000))+
    coord_flip()+
    geom_col(position = position_dodge(1))+
    geom_errorbar(aes(x=delay,
                        group= n_negatives,
                        ymin=Lower*10000,
                        ymax=Upper*10000
    ),
    colour="black",
    position = position_dodge(1),
    width=.2
    )+
    scale_y_continuous(trans = pseudo_log_trans(1e7, 10),
                       breaks = c(0,seq(5e3,2.5e4,by=5e3)))+
    labs(x="",y="Tests used \nper 10,000 infected individuals",
         colour="Number of consecutive days of negative tests required for release")
  
)

plot_2a+plot_2b+plot_2c+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect")&
  scale_x_discrete(labels=delay_lab,limits=rev)&
  scale_color_manual(values=c("#22577A",
    "#38A3A5",
    "#57CC99",
    "#80ED99"),
    labels=n_negatives_lab,
    guide=guide_legend(title.position = "top"))&
  scale_fill_manual(values=c("#22577A",
                              "#38A3A5",
                              "#57CC99",
                              "#80ED99"),
                     labels=n_negatives_lab,
                     guide=guide_legend(title.position = "top"))&
  labs(colour="Number of consecutive days of negative tests required for release",
       fill="Number of consecutive days of negative tests required for release")&
  theme_minimal()&
  theme(panel.grid.major.y = element_blank(),
              axis.ticks = element_line(),
              axis.line.x.bottom = element_line(),
              axis.line.y.left = element_line(),
              legend.position = "bottom")&
  facet_rep_grid(lower_inf_thresh ~ higher_detect_thresh,labeller=label_both)

ggsave("output/main_plot2_vacc.png",width=500,height=250,units="mm",dpi=600,bg = "white")

source("scripts/additional_plots.R")


#relative benefit of extra day and test vs. extra day no test
baseline <- prop_averted %>% 
  filter(n_negatives==0) %>% 
  nest(-c(sim,variant,scenario_id,lower_inf_thresh,higher_detect_thresh,n_negatives,delay,inf_community))

prop_averted %>% 
  left_join(baseline,by=c("sim","variant","lower_inf_thresh","higher_detect_thresh")) %>% 
  filter(delay.x==delay.y) %>% 
  mutate(change=inf_community.x/inf_community.y,
         change=ifelse(is.na(change),1,change)) %>% 
  group_by(delay.x,n_negatives.x) %>% 
  mutate(x=map(change,summarise_func))
  ggplot(aes(x=delay.x,y=change,group=n_negatives.x,colour=n_negatives.x))+
  geom_point(position = position_jitterdodge(dodge.width = 1,jitter.height = 0),alpha=0.1)
    
