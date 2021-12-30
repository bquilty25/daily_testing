
(plot_2a_supp <- plot_2a_dat %>%  
    #filter.(variant%in%c("vacc")) %>% 
    ggplot()+
    coord_flip()+
    geom_pointrange(aes(x=fct_rev(delay),
                        group= fct_rev(n_negatives),
                        colour=n_negatives,
                        y=mean,
                        ymin=lower,
                        ymax=upper),
                    position = position_dodge(0.5),
                    fatten = 1)+
    scale_y_continuous(breaks=breaks_width(1),limits=c(0,NA))+
    labs(x="",y="Days saved vs. 10 days isolation\nper individual",
         colour="Number of consecutive days of negative tests required for release")
)

(plot_2b_supp <- plot_2b_dat %>% 
    ggplot()+
    coord_flip()+
    geom_pointrange(aes(x=fct_rev(delay),
                        group= fct_rev(n_negatives),
                        colour=n_negatives,
                        y=mean,
                        ymin=lower,
                        ymax=upper
    ),
    position = position_dodge(0.5),
    fatten = 1)+
    scale_y_continuous(trans="pseudo_log",breaks=c(0,10,100,1000,10000))+
    labs(x="",y="Days infectious in the community\nper 10,000 infected individuals",
         colour="Number of consecutive days of negative tests required for release")
  
)


(plot_2c_supp <- plot_2c_dat %>% 
    ggplot()+
    coord_flip()+
    geom_pointrange(aes(x=fct_rev(delay),
                        group= fct_rev(n_negatives),
                        colour=n_negatives,
                        y=mean,
                        ymin=lower,
                        ymax=upper
    ),
    position = position_dodge(0.5),
    fatten = 1)+
    scale_y_continuous()+
    labs(x="",y="Tests used \nper 10,000 infected individuals",
         colour="Number of consecutive days of negative tests required for release")
  
)

plot_2a_supp/plot_2b_supp/plot_2c_supp+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect")&
  #scale_color_brewer(palette="Set2")&
  scale_color_manual(values=c("#22577A",
                              "#38A3A5",
                              "#57CC99",
                              "#80ED99"))&
  facet_rep_grid(vacc_status~omicron,
                 labeller=labeller(variant=c("vacc"="Vaccinated",
                                             "unvacc"="Unvaccinated",
                                             "omicron_vacc_est"="Omicron vaccinated\n(assumed)",
                                             "omicron_unvacc_est"="Omicron unvaccinated\n(assumed)")))&
  guides(colour=guide_legend(title.position="top"))&
  theme_minimal()&
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "bottom")

ggsave("output/plot2_supp.png",width=210,height=297,units="mm",dpi=600,bg = "white")
