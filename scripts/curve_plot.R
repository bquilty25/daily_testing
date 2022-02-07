## curves for paper

set.seed(2020)

trajectories_to_plot <- traj_ %>% #filter(variant=="vacc") %>% 
  nest_by.(sim,variant,lower_inf_thresh, higher_detect_thresh) %>% 
  slice_sample.(100,.by=c(variant,lower_inf_thresh, higher_detect_thresh)) %>% 
  unnest.(data)

curve_plot_a <- trajectories_to_plot %>% 
  ggplot()+
  geom_line(aes(x=t,y=vl,colour=culture_p,group=sim),alpha=0.1)+
  geom_point(aes(x=t,y=vl,colour=culture_p),alpha=0.5)+
  scale_x_continuous(name="Days since exposure",breaks=breaks_width(1))+
  scale_color_viridis_b(name="Probability of positive culture (infectiousness)",option = "inferno",begin=0.2,end=0.8,guide=guide_colorsteps(show.limits = T,barwidth=20,title.position="top",title.hjust=0.5),breaks=seq(0,1,by=0.1),label=percent_format(accuracy = 1L))+
  scale_y_continuous(name="log10 RNA copies/ml",breaks=breaks_width(2),limits = c(0,NA),expand = c(0, 0.1))+
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "bottom")

curve_plot_b <- trajectories_to_plot %>% 
  ggplot()+
  geom_line(aes(x=t,y=vl,colour=test_p,group=sim),alpha=0.1)+
  geom_point(aes(x=t,y=vl,colour=test_p),alpha=0.5)+
  scale_x_continuous(name="Days since exposure",breaks=breaks_width(1))+
  scale_color_viridis_b(name="Probability of positive LFT",option = "inferno",begin=0.2,end=0.8,guide=guide_colorsteps(show.limits = T,barwidth=20,title.position="top",title.hjust=0.5),breaks=seq(0,1,by=0.1),label=percent_format(accuracy = 1L))+
  scale_y_continuous(name="log10 RNA copies/ml",breaks=breaks_width(2),limits = c(0,NA),expand = c(0, 0.1))+
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "bottom")

curve_plot_a/curve_plot_b+plot_annotation(tag_levels = "A")&
  facet_rep_grid(variant~lower_inf_thresh + higher_detect_thresh,labeller=label_both)


ggsave("output/curve_plot.png",height=297,width=210,units="mm",dpi=600,bg = "white")


#Detection and infectiousness assumptions
det_plot <- pickering %>% 
  pivot_longer.(cols=c(culture,Innova,"SureScreen F", Encode),values_to = "prob") %>% 
  filter(name=="Innova") %>% 
  select(-name) %>% 
  mutate(vl_2.5=vl+2.5) %>% 
  pivot_longer(c(vl,vl_2.5)) %>% 
  ggplot(aes(x=value,
             y=prob,
             group=name,
             colour=name))+
  geom_smooth(method=glm,method.args=list(family="binomial"),se=F)+
  labs(x="Viral load (log10 RNA copies/ml)",y="Probability of detection")+
  scale_color_manual(values=c("#ff7f00","#fdbf6f"),
                     labels=c("Baseline","Harder to detect \n(Baseline + 2.5 log10 RNA copies/ml)"),
                     guide=guide_legend("Assumption", title.position = "top",title.hjust=0.5))

#Detection and infectiousness assumptions
inf_plot <- pickering %>% 
  pivot_longer.(cols=c(culture,Innova,"SureScreen F", Encode),values_to = "prob") %>% 
  filter(name=="culture") %>% 
  select(-name) %>% 
  mutate(vl_2.5=vl-2.5) %>% 
  pivot_longer(c(vl,vl_2.5)) %>% 
  ggplot(aes(x=value,
             y=prob,
             group=name,
             colour=name))+
  geom_smooth(method=glm,method.args=list(family="binomial"),se=F)+
  labs(x="Viral load (log10 RNA copies/ml)",y="Probability of culture/infectiousness")+
  scale_color_manual(values=c("#33a02c",
                              "#b2df8a"),
                     labels=c("Baseline","Lower infectious dose\n(Baseline - 2.5 log10 RNA copies/ml)"),
                     guide=guide_legend("Assumption", title.position = "top",title.hjust=0.5))
  

det_plot/inf_plot+
  plot_annotation(tag_levels = "A")&
  scale_x_continuous(limits = c(0,12.5))&
  theme_minimal()&
  theme(panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "bottom")

ggsave("output/det_inf_assump.png",height=297*0.75,width=210*0.75,units="mm",dpi=600,bg = "white")

