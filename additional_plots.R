#create plots
plot_dat <-  prop_averted %>% 
  #filter.(variant=="unvacc") %>% 
  mutate.(across.(c(iso_dur,n_negatives,delay,tests_used,variant),as.factor))

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
  ggh4x::facet_nested(delay+n_negatives~variant+lower_inf_thresh+higher_detect_thresh,nest_line = T,labeller=label_both)&
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

ggsave("output/main_plot.png",width=210*2,height=297*2,units="mm",dpi=600,bg = "white")
