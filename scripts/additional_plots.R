plot_dat <- prop_averted 

plot_2a_dat <- plot_dat[,
                        {x=Hmisc::smean.cl.boot(inf_community>0,B = 1000);as.list(x)},
                        #{x=boot_func(inf_community,stat="prop_greater");as.list(x)},
                        #{x=Hmisc::smean.cl.boot(inf_community>0,B = 1000);as.list(x)},
                        by=c("variant","n_negatives","delay","lower_inf_thresh")]

(plot_2a <- plot_2a_dat %>% 
    ggplot(aes(x=delay,
               group= n_negatives,
               fill=n_negatives,
               y=Mean))+
    coord_flip()+
    geom_col(position = position_dodge(0.8),width=0.5)+
    geom_linerange(aes(x=delay,
                       group= n_negatives,
                       ymin=Lower,
                       ymax=Upper
    ),
    colour="black",
    position = position_dodge(0.8),
    width=.2)+
    scale_y_continuous(trans = pseudo_log_trans(1, 10),
                       #breaks = c(0,lseq(1,1000,length.out = 4)),
                       labels=scales::percent_format(trim=F,accuracy = 1L),
                       breaks=pretty_breaks(n=4)
                       #breaks=function(x)c(min(x),c(0,0.001,0.01,0.1),max(x))
    )+
    labs(x="",
         y="Proportion infectious after release")
)

plot_2b_dat <- plot_dat[,
                        {x=summarise_func(days_saved);as.list(x)},
                        #{x=Hmisc::smean.cl.boot(days_saved,B = 1000);as.list(x)},
                        by=c("variant","n_negatives","delay","lower_inf_thresh")]

plot_2b <- plot_dat %>% 
  ggplot()+
  geom_count(aes(x=delay,
                 y=days_saved,
                 group=n_negatives,
                 colour=n_negatives,
                 size=after_stat(prop)),
             position = position_dodge(0.8))+
  geom_pointrange(data=plot_2b_dat,aes(x=delay,
                                       y=Mean,
                                       ymin=Lower,
                                       ymax=Upper,
                                       group=n_negatives),
                  colour="black",
                  fatten = 0.5,
                  position = position_dodge(0.8))+
  coord_flip()+
  scale_y_continuous(breaks=breaks_width(1),limits=c(0,7),expand = expansion(add=0.5))+
  labs(x="",y="Number of days saved versus a 10 day isolation",
       colour="Number of consecutive days of negative tests required for release")

plot_2c_dat <- plot_dat[,
                        {x=summarise_func(tests_used);as.list(x)},
                        #{x=Hmisc::smean.cl.boot(tests_used,B = 1000);as.list(x)},
                        by=c("variant","n_negatives","delay","lower_inf_thresh")]

plot_2c <- plot_dat %>% 
  ggplot()+
  geom_count(aes(x=delay,
                 y=tests_used,
                 group=n_negatives,
                 colour=n_negatives,
                 size=after_stat(prop)),
             position = position_dodge(0.8))+
  geom_pointrange(data=plot_2c_dat,aes(x=delay,
                                       y=Mean,
                                       ymin=Lower,
                                       ymax=Upper,
                                       group=n_negatives),
                  colour="black",
                  fatten = 0.5,
                  position = position_dodge(0.8))+
  coord_flip()+
  scale_y_continuous(breaks=breaks_width(1),limits=c(0,6),expand = expansion(add=0.5))+
  labs(x="",y="Number of tests used",
       colour="Number of consecutive days of negative tests required for release")

plot_2a+plot_2b+plot_2c+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect")&
  scale_x_discrete(labels=delay_lab,limits=rev)&
  scale_color_manual(values=lighten(c("#22577A",
                                      "#38A3A5",
                                      "#57CC99",
                                      "#80ED99"),amount=0.3),
                     labels=n_negatives_lab,
                     guide=guide_legend(title.position = "top",title.hjust = 1,reverse = T))&
  scale_fill_manual(values=lighten(c("#22577A",
                                     "#38A3A5",
                                     "#57CC99",
                                     "#80ED99"),amount=0.3),
                    labels=n_negatives_lab,
                    guide="none",
                    #guide=guide_legend(title.position = "top",title.hjust = 1,reverse = T)
  )&
  scale_size_area("Proportion",
                  max_size = 5,
                  guide=guide_legend(title.position = "top",title.hjust = 0.5),
                  labels=percent_format(accuracy = 1L))&
  labs(colour="Number of consecutive days of negative tests required for release",
       fill="Number of consecutive days of negative tests required for release")&
  theme_minimal()&
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_line(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "bottom",
        legend.title.align=0.5,
        panel.spacing = unit(0.001, "lines"))&
facet_rep_grid(~lower_inf_thresh,labeller=labeller(lower_inf_thresh=inf_thresh_lab),scales="free_x")

ggsave("output/main_plot_sens.png",width=350,height=150,units="mm",dpi=600,bg = "white")

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
