source("packages.R")
source("plot_functions.R")

detection_curve <- read.csv("data/pos_curve.csv") %>% 
  mutate(assay="PCR")

detection_prob <- function(upper,lower){
  if (any(is.na(c(upper, lower)))){
    return(list(shape1 = NA,
                shape2 = NA))
  }  
  suppressMessages(rriskDistributions::get.beta.par(
    q = c(lower,  upper),
    p = c(0.025, 0.975), 
    show.output = F, plot = F) %>%
      as.list)
}

detection_parameters <- 
  detection_curve %>% as_tibble() %>% 
  mutate(beta_parms = pmap(.l=list(upper = upper,
                                   lower = lower),
                           .f=detection_prob))

plot_curve <- ggplot(data =detection_curve,
                     aes(x = days_since_infection,
                         y = median)) +
  geom_ribbon(alpha = 0.5,
              aes(ymin = lower,
                  ymax = upper,
                  fill = fct_reorder2(assay,days_since_infection,median)),show.legend =F) +
  geom_line(aes(colour = fct_reorder2(assay,days_since_infection,median)),show.legend = F) +
  xlab("Time since exposure (days)") +
  ylab("Probability of positive PCR test") +
  scale_fill_manual(name="",values=col_pal[2])+
  scale_color_manual(name="",values=col_pal[2])+
  plotting_theme +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format())


infectivity_curve <- data.frame(x=-15:15) %>% 
  mutate(y=dgamma(x=x+infect_shift,shape=infect_shape,rate=infect_rate)) %>% 
  ggplot()+
  geom_line(aes(x=x,y=y,colour="Infectivity profile"),show.legend = F)+
  labs(x="Time since onset (days)",
       y="Infectivity (density)")+
  scale_colour_manual(name="",values = col_pal[3])+
  plotting_theme

plot_curve + infectivity_curve + plot_annotation(tag_levels = "A")

save_plot(device="png",width=300,height=150)
