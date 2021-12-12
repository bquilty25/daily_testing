x <- make_incubation_times(1000,pathogen = pathogen,asymp_parms = asymp_fraction) %>% 
  filter(type=="symptomatic") %>% 
  crossing(assay=c("LFA","PCR")) 

x_ <- curves %>%  
    distinct(idx,iter_lfa,days_since_infection,value,peak_timing,.keep_all = T) %>% 
    inner_join(x, by = c("idx","assay")) %>%
    #shift curves by onset 
    mutate(days_since_infection_shift=days_since_infection-peak_timing) %>% 
    nest(data = -c(idx,assay)) %>%
    mutate(model = map(data, ~approxfun(x=.$days_since_infection_shift,
                                        y=.$value))) %>%
    crossing(test_t=seq(-10,10,by=0.1)) %>% 
    mutate(infectious_prob=dgamma(x = test_t+infect_shift,shape = infect_shape,rate=infect_rate)/0.1511372,
           infectious=as.logical(rbinom(n=n(),size=1,prob=infectious_prob))) %>%  
    # mutate(upper_threshold = 30,
    #            test_p      = calc_sensitivity(model           = model,
    #                                           x               = test_t, 
    #                                           upper_threshold = upper_threshold)) %>% 
    # select(-c(upper_threshold,model)) %>% 
     mutate(screen = runif(n(), 0, 1)) %>% 
    mutate(test_p=case_when(infectious&assay=="LFA"~0.768,
                            infectious&assay=="PCR"~0.99,
                            TRUE~0)) %>% 
    mutate(test_label       = detector(pcr = test_p,  u = screen)) %>% 
    mutate(detected_infectious = infectious&test_label)
    

x_ %>% group_by(assay) %>% filter(between(test_t,-5,5)) %>% 
  summarise(mean(detected_infectious))
