## curves for paper

set.seed(2020)
traj <- variant_char %>% 
  filter.(variant%in%c("omicron")) %>% 
  group_split.(variant) %>% 
  map.(~make_trajectories(n_sims = 1000,
                          seed=seed,
                          variant_info=.x,
                          browsing = F)) %>% 
  bind_rows.()

#Calculate daily infectiousness and test positivity, remove never-infectious
traj_ <- traj %>% filter.(variant%in%c("omicron")) %>% 
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = ceiling(max(end))
  )))  %>%
  unnest.(infectiousness) %>%
  crossing(
    lower_inf_thresh = c(TRUE, FALSE)
  ) %>%
  rename.("ct"=vl) %>% 
  mutate.(
    vl=convert_Ct_logGEML(ct),
    culture_p        = stats::predict(
      object = inf_model_choice(lower_inf_thresh),
      type = "response",
      newdata = tidytable(vl = vl)
    ),
    infectious_label = rbernoulli(n = n(),
                                  p = culture_p),
    test_p           = stats::predict(
      object =  innova_mod,
      type = "response",
      newdata = tidytable(vl = vl)
    ),
    test_label       = rbernoulli(n = n(),
                                  p = test_p),
    .by = c(lower_inf_thresh)
  ) %>%
  replace_na.(list(test_label       = FALSE,
                   infectious_label = FALSE))

#trajectories_to_plot <- traj_ %>% #filter.(lower_inf_thresh==F) %>%  #filter(variant=="omicron") %>% 
#  nest_by.(sim,variant,lower_inf_thresh) %>% 
#  slice_sample.(20,.by=c(variant,lower_inf_thresh)) %>% 
#  unnest.(data)

curve_plot_a <- traj_ %>% 
  ggplot()+
  geom_line(aes(x=t,y=ct,colour=culture_p,group=sim),alpha=0.5)+
  scale_color_viridis_c(name="Probability of positive culture",
                        option = "inferno",
                        begin=0.2,
                        end=0.8,
                        guide=guide_colorsteps(show.limits = T,
                                               barwidth=20,
                                               title.position="top",
                                               title.hjust=0.5),
                        breaks=seq(0,1,by=0.1),
                        label=percent_format(accuracy = 1L)
                        )+
  ggnewscale::new_scale_colour()+
  #geom_point(aes(x=t,y=ct,colour=infectious_label),alpha=0.5)+
  scale_color_viridis_d(name="Probability of positive culture",
                        option = "inferno",
                        begin=0.2,
                        end=0.8
  )

curve_plot_b <- traj_ %>% 
  ggplot()+
  geom_line(aes(x=t,y=ct,colour=test_p,group=sim),alpha=0.5)+
  scale_color_viridis_c(name="Probability of positive LFT",
                        option = "magma",
                        begin=0.2,
                        end=0.8,
                        guide=guide_colorsteps(show.limits = T,
                                               barwidth=20,
                                               title.position="top",
                                               title.hjust=0.5),
                        breaks=seq(0,1,by=0.1),
                        label=percent_format(accuracy = 1L))+
  ggnewscale::new_scale_colour()+
  #geom_point(aes(x=t,y=ct,colour=test_label),alpha=0.5)+
    scale_color_viridis_d(name="Probability of positive LFT",
                          option = "magma",
                          begin=0.2,
                          end=0.8
                  )


a <- curve_plot_a/curve_plot_b+plot_annotation(tag_levels = "A")&
  scale_x_continuous(name="Days since first detectable by PCR",breaks=breaks_width(5))&
  scale_y_reverse(name = "Ct", 
                  sec.axis = sec_axis(~convert_Ct_logGEML(.), 
                                                 name="log10 RNA copies per ml"),
                  expand = c(0, 0.01))&
  facet_rep_grid(variant~lower_inf_thresh ,labeller=label_both)&
  theme_minimal()&
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.line.y.right = element_line(),
        legend.position = "bottom")


ggsave("output/curve_plot.png",height=297,width=210,units="mm",dpi=600,bg = "white")

b <- traj_ %>% 
  filter(sim%in%c(1:25),lower_inf_thresh==F) %>% 
  ggplot()+
  geom_tile(aes(x=t,y=sim,group=sim,fill=infectious_label),colour="white",size=0.5,alpha=0.75)+
  scale_fill_manual(name="Culture",values=c("#edf8e9","#74c476"),labels=c("No","Yes"))+
  ggnewscale::new_scale_fill()+
  geom_point(aes(x=t,y=sim,group=sim,fill=test_label),shape=21,size=2)+
  #scale_shape_manual(values=c("-","+"))
  scale_fill_manual("LFT",values=c("white","#006d2c"),labels=c("Negative","Positive"))+
  scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks=breaks_width(5),limits=c(NA,20.5))+
  coord_cartesian(expand=F)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_line(),
        axis.text.y = element_blank(),
        axis.title.y= element_blank(),
        axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line(),
        #axis.line.y.left = element_line(),
        #axis.line.y.right = element_line(),
        legend.position = "bottom")

c <- traj_ %>% 
  pivot_longer(c(test_label,infectious_label)) %>% 
  group_by(sim,variant,lower_inf_thresh,name) %>% 
  filter(value) %>% 
  summarise(min_t=min(t),max_t=max(t)) %>% 
  mutate(diff=max_t-min_t)  %>% 
  ggplot()+
  geom_histogram(aes(x=diff,y=stat(density),fill=name),binwidth=1)+
  stat_pointinterval(aes(x=diff))+
  labs(x="Duration of positivity",y="Probability density")+
  scale_x_continuous(breaks=breaks_width(2))+
  scale_y_continuous(limits=c(0,0.25))+
  scale_fill_manual(values=c("#FCA50AFF","#FE9F6DFF"),guide="none")+
  facet_rep_wrap(lower_inf_thresh~name,ncol=1,labeller=labeller(name=c("infectious_label"="Culture",
                                                "test_label"="LFT")))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.line.y.right = element_line(),
        legend.position = "bottom")

#Detection and infectiousness assumptions
pickering %>% 
  pivot_longer.(cols=c(culture,Innova,"SureScreen F", Encode),values_to = "prob") %>% 
  filter(name%in%c("culture","Innova")) %>% 
  ggplot(aes(x=vl,
             y=prob,
             group=name,
             colour=name))+
  geom_smooth(method=glm,method.args=list(family="binomial"),se=F)+
  labs(x="Viral load (log10 RNA copies/ml)",y="Probability")+
  scale_color_manual(values=c("#ff7f00",
                              "#33a02c"),
                     labels=c("Culture","LFT"),
                     guide=guide_legend("Assumption", title.position = "top",title.hjust=0.5))+
  scale_x_continuous(breaks=breaks_width(2.5))+
  theme_minimal()+
  theme(panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "bottom")

ggsave("output/det_inf_assump.png",height=150,width=210,units="mm",dpi=600,bg = "white")

