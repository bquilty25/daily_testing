## curves for paper

set.seed(2020)

trajectories_to_plot <- traj_ %>%
  nest_by.(sim,idx,variant) %>% 
  slice_sample.(100,.by=c(variant)) %>% 
  unnest.(data) %>% 
  mutate.(culture=stats::predict(culture_mod, type = "response", newdata = tidytable(vl=vl)),
          infectious_label = rbernoulli(n=n(),p=culture)) %>% 
  arrange.(sim,idx,variant) %>% 
  mutate.(variant=fct_relevel(variant,c("vacc",
                              "unvacc",
                              "omicron_vacc_est",
                              "omicron_unvacc_est")))

curve_plot_a <- trajectories_to_plot %>% 
  ggplot()+
  geom_line(aes(x=test_t,y=vl,colour=culture,group=interaction(sim,idx)),alpha=0.1)+
  geom_point(aes(x=test_t,y=vl,colour=culture),alpha=0.5)+
  scale_x_continuous(name="Days since exposure",breaks=breaks_width(1))+
  scale_color_viridis_b(name="Probability of positive culture (infectiousness)",option = "magma",begin=0.2,end=0.8,guide=guide_colorsteps(show.limits = T,barwidth=12,title.position="top",title.vjust=1),label=label_percent())+
  scale_y_continuous(name="log10 RNA copies/ml",breaks=breaks_width(2),limits = c(0,NA),expand = c(0, 0))+
  facet_rep_wrap(~variant,
                 labeller=labeller(variant=c("vacc"="Vaccinated",
                                             "unvacc"="Unvaccinated",
                                             "omicron_vacc_est"="Omicron vaccinated (assumed)",
                                             "omicron_unvacc_est"="Omicron unvaccinated (assumed)")))+
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "bottom")

curve_plot_b <- trajectories_to_plot %>% 
  ggplot()+
  geom_line(aes(x=test_t,y=vl,colour=test_p,group=interaction(sim,idx)),alpha=0.1)+
  geom_point(aes(x=test_t,y=vl,colour=test_p),alpha=0.5)+
  scale_x_continuous(name="Days since exposure",breaks=breaks_width(1))+
  scale_color_viridis_b(name="Probability of positive LFT",option = "inferno",begin=0.2,end=0.8,guide=guide_colorsteps(show.limits = T,barwidth=10),label=label_percent())+
  scale_y_continuous(name="log10 RNA copies/ml",breaks=breaks_width(2),limits = c(0,NA),expand = c(0, 0))+
  facet_rep_wrap(~variant,
                 labeller=labeller(variant=c("vacc"="Vaccinated",
                                             "unvacc"="Unvaccinated",
                                             "omicron_vacc_est"="Omicron vaccinated (assumed)",
                                             "omicron_unvacc_est"="Omicron unvaccinated (assumed)")))+
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "bottom")

curve_plot_a/curve_plot_b+plot_annotation(tag_levels = "A")

ggsave("curve_plot.png",width=210,height=297,units="mm",dpi=600,bg = "white")
