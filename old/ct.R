pacman::p_load(tidyverse,mgcv,brms,tidybayes,modelr,ggeffects)

dt <- read_csv("data/dt.csv")

never_pos <- function(df){
  x_ <- df %>% filter(pcr_result)
  
  if(nrow(x_)==0L){
    res <- "never pos"
  }else{
    res <- "pos"
  }
  
  return(res)
}

 dt <- dt %>% 
  mutate(num_id   = factor(num_id),
         # zeros and NAs mean the same thing here?
         ct       = case_when(is.na(ct) ~ 40, # when missing, it's undetected
                              ct == 0   ~ 40, # when zero it's undetected
                              TRUE      ~ ct),
         inv_ct   = 1/ct,
         prob_det = cut(ct,
                        breaks = c(-Inf, 18, 21.5, 25, 28, 31, 34.5, 40, Inf),
                        labels = c(1, 0.98, 0.96, 0.89, 0.65, 0.32, 0.23, 0.01)),
         prob_det = parse_number(as.character(prob_det))) %>% 
  group_by(num_id) %>% 
  do(add_row(., num_id = unique(.$num_id), x = 0, ct = 40)) %>% 
  mutate(inv_ct=replace_na(inv_ct,0)) %>% 
  nest(-c(num_id)) %>% 
  mutate(res=map(.x=data,.f=never_pos)) %>% 
  filter(res!="never_pos") %>% 
  unnest(data) %>% 
  filter(x>=0)

ggplot(dt,aes(x=x,y=ct,colour=factor(num_id),group=factor(num_id)))+
  geom_point()+
  geom_smooth(method="loess",se=F,size=0.5)+
  #geom_hline(aes(yintercept=35,linetype="PCR threshold"))+
  #geom_hline(aes(yintercept=27,linetype="LFA threshold"))+
  # scale_y_log10("Inferred viral load (RNA copies/mL)",
  #               breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  scale_colour_viridis_d()+
  #plotting_theme+
  theme(legend.position = "right")
  
  m1 <- brm(bf(ct ~ s(x)+s(num_id,bs="re")),
      data = dt,control = list(adapt_delta=0.99))
summary(m1)

conditional_effects(m1, spaghetti = T,nsamples=100,error=FALSE,effects="x")

save_plot(device="png",width=250,height = 150)

nd <- data.frame(x=seq(0,40,by=0.1))


f <-
  fitted(m1, newdata = nd) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(nd) 

