pacman::p_load(tidyverse)

#### LIVERPOOL INNOVA ANALYSIS ----
innova_liv <- tribble(~left,~right,~pos,~neg,
                      30,35,1,18,
                      25,30,3,9,
                      20,25,18,5,
                      15,20,17,2) %>% 
  pivot_longer(c(pos,neg))  

# random sample from binned intervals 
innova_liv_sim <- innova_liv %>% 
  group_by(left,right,name) %>% 
  uncount(value) %>% 
  ungroup() %>% 
  mutate(id=row_number()) %>% 
  crossing(sim=c(1:100)) %>% 
  mutate(ct=runif(n=n(),min = left,max=right),
         AG=as.numeric(as.factor(name))-1)

# Logisitic prediction models
innova_liv_mod <- glm(AG~ct,data=innova_liv_sim,family="binomial")

innova_liv_mod_plus_5 <- glm(AG~ct,data=innova_liv_sim %>% mutate(ct=ct+2.5),family="binomial")

innova_liv_mod_minus_5 <- glm(AG~ct,data=innova_liv_sim %>% mutate(ct=ct-2.5),family="binomial")

bind_rows(innova_liv_sim %>% mutate(name="Innova LFT",type="Actual (Innova)",value=AG)) %>%
  bind_rows(innova_liv_sim %>% mutate(name="Innova LFT + 2.5 CT adjustment",type="Sensitivity analysis",value=AG,ct=ct+2.5)) %>%
  bind_rows(innova_liv_sim %>% mutate(name="Innova LFT - 2.5 CT adjustment",type="Sensitivity analysis",value=AG,ct=ct-2.5)) %>%
  ggplot(aes(x=ct,y=value,group=name,linetype=type))+
  geom_smooth(method = "glm",se=F,method.args=list(family="binomial"),size=1)+ 
  scale_color_brewer(type="qual",palette = "Set2")+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA),
        legend.position = "bottom")+
  scale_x_continuous(breaks=scales::breaks_width(5))+
  labs(x="Cycle Threshold",
       y="Probability of detection",
       colour="",
       linetype="")

ggsave("results/test_comparison.png",height=200,width=180,dpi=400,units="mm")
