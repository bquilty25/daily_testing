trajectories_to_plot <- traj %>%
  #filter(type=="symptomatic") %>% 
  slice_sample.(n=100,.by=c(variant,type)) %>% 
  #ungroup() %>% 
  mutate.(pred = map(.x = m, 
                    ~data.frame(t = seq(0, 25, by=0.25)) %>%
                      mutate(ct = .x(t)))) %>%
  unnest.(pred) %>%
  mutate.(LFT = stats::predict(innova_mod,
                              newdata =
                                data.frame(ct = ct),
                              type = "response"
  ),
  culture=stats::predict(culture_mod,
                         newdata =
                           data.frame(ct = ct),
                         type = "response"
  ),
  PCR=ifelse(ct<35,1,0)) %>% 
  arrange.(variant,type)

trajectories_to_plot %>% 
  #pivot_longer(cols=c(ct,culture,LFT)) %>%  
  ggplot()+
  geom_line(aes(x=t,
                y=ct,
                group=interaction(idx,sim,variant,type)),
            alpha=0.5)+
  scale_fill_brewer()+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  #scale_color_manual(values="black",name="",guide=F)+
  #scale_x_continuous(name="Days since exposure",breaks = breaks_width(7))
  scale_y_reverse(# Features of the first axis
    name = "Cycle\nthreshold")+
  facet_grid(variant~fct_rev(type))+
  limits(y=c(40,NA))
