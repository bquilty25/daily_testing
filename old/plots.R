source("packages.R")
source("utils.R")
source("plot_functions.R")
#results_name <- "sum_results"

results<- read_results(results_name)

main <- ribbon_plot(
  results %>% filter(waning == "adhere_100",
                         index_test_delay == 2,
                         delay_scaling == 1) %>% 
    #calculate time until release from exposure for each scenario
    mutate(time_since_exp=ifelse(stringency=="none",
                                 yes=quar_dur,
                                 no=quar_dur + results_delay * delay_scaling)),
  by_type = F,
  custom_facets =  ~ stringency ,
  y_labels = infectivity_labels["infectivity_averted"],
  #ribbon = T,
  colour_var = "assay")
main

ggsave("results/Figure_2_averted.png", width=210, height=150,units="mm",dpi=320)

delay_scaling <- ribbon_plot(
  results %>% 
    filter(waning=="waning_none",
      index_test_delay==2,
      #delay_scaling == 1
    ) %>% 
    mutate(delay_scaling=fct_rev(as.factor(delay_scaling))) %>% View(),
  by_type = F,
  custom_facets =  ~  delay_scaling+stringency,
  y_labels = infectivity_labels[c("infectivity_averted")],
  ribbon = T,
  colour_var = "stringency")
delay_scaling

ggsave("results/delay_scaling.png", width=210, height=100,units="mm",dpi=320)

waning <- ribbon_plot(
  read_results(results_name) %>% 
    filter(#waning=="waning_none",
           index_test_delay==2,
           delay_scaling == 1
    ) #%>% 
    #mutate(waning=fct_rev(as.factor(waning))
          # ),
  by_type = F,
  custom_facets =  ~  waning+stringency,
  y_labels = infectivity_labels[c("infectivity_averted")],
  ribbon = T,
  colour_var = "stringency")
waning

ggsave("results/waning_severe.png", width=210, height=100,units="mm",dpi=320)

post <- ribbon_plot(
  results_dat %>% filter(waning == "waning_none",
                         index_test_delay == 2,
                         delay_scaling == 1
                         ) %>% 
    #calculate time until release from exposure for each scenario
    mutate(time_since_exp=ifelse(stringency=="none",
                                 yes=quar_dur,
                                 no=quar_dur + results_delay * delay_scaling)) %>% 
    mutate(delay_scaling=fct_rev(factor(delay_scaling))),
  by_type = T,
  custom_facets = . ~ stringency ,
  y_labels = infectivity_labels["infectivity_post"],
  ribbon = T,
  colour_var = "stringency")
post <- post + labs(y="Transmission potential after quarantine")


save_plot(post,
          dpi = 320, 
          device = "png",
          prefix = "post",
          base = "plot", 
          width = 210, 
          height = 140)

#adherence sensitivity
adhere_sens_plot <- read_results("adhere_sens") %>% 
  test_labeller() %>% 
  filter(yvar=="infectivity_averted",type=="all") %>% 
  mutate(waning=str_sub(waning,-2)) %>% 
  mutate(waning=ifelse(waning=="00",100,waning),waning=as.factor(as.numeric(waning))) %>% 
  ggplot()+
  geom_hline(aes(yintercept=0.07017903,linetype="Median transmission potential averted assuming 10% adherence to 14-day quarantine in UK (Smith et al. 2020)"))+
  geom_line(aes(group=waning,colour=waning,x=time_since_exp,y=`50%`),show.legend = F)+
  directlabels::geom_dl(aes(x=time_since_exp,y=`50%`,label = paste0(waning,"%"),colour=waning), method = list(dl.combine("last.qp")),size=0.2) +
  facet_grid(~stringency)+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA),
        legend.position = "bottom")+
  labs(x="Time since exposure (days)",y="Median transmission potential averted")+
  scale_x_continuous(limits=c(0,20),breaks=breaks_width(7))+
  scale_y_continuous(limits = c(0,1),breaks=breaks_width(0.25))+
  scale_colour_viridis_d(name="Adherence to\nquarantine (%)",option = "magma",direction = -1,begin = 0.2,end=0.8)+
  scale_linetype_manual(name="",values="dotted")

save_plot(adhere_sens_plot,
          dpi = 320, 
          device = "png",
          prefix = "adhere_sens",
          base = "plot1", 
          width = 210, 
          height = 100)
