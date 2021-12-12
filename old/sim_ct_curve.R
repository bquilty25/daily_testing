pacman::p_load(EnvStats)


sim_ct <- data.frame(idx=1:100) %>% 
    crossing(start=0,type=c("symptomatic","asymptomatic")) %>% 
    mutate(end=case_when(type=="symptomatic" ~ rnorm(n=n(),mean=17,sd=2),
                       type=="asymptomatic" ~ rnorm(n=n(),mean=17*0.6,sd=1))) %>% 
    mutate(peak=rlnormTrunc(n=n(),
                            meanlog=1.63,
                            sdlog=0.5,
                            min=0,
                            max=end)) %>% 
    pivot_longer(cols=-c(idx,type),values_to = "x") %>% 
    mutate(y=case_when(name=="start"~40,
                       name=="end"~40,
                     name=="peak"~rnorm(n=n(),mean=25,sd=5))) 

sim_ct%>% 
  ggplot(aes(x=x,y=y,group=idx))+
  geom_point() +
  stat_smooth(geom='line', alpha=0.2, se=FALSE)+
  facet_grid(~type)

save_plot(device="png",prefix = "sim_ct",width = 210,height=150)
  
  
make_newdata <- function(data,model){
 # browser()
  newdata <- data.frame(x=seq(from=min(data$x),to=max(data$x),0.1))
  
  newdata$y_pred <- predict(model,newdata)
  
  return(newdata)
}
  
  #loess
 models <-  sim_ct %>%
    nest(-c(idx,type)) %>%  
    dplyr::mutate(
    # Perform loess calculation on each individual 
    m = purrr::map(data, loess,
                   formula = y ~ x),
    # Retrieve the fitted values from each model
    pred_df=purrr::map2(.x=data,.y=m,make_newdata))

 results <- models %>%
   tidyr::unnest(pred_df)
 
 # Plot with loess line for each group
 ggplot(results, aes(x = x, y = y_pred, group = idx, colour = type)) +
   geom_line(alpha=0.5)+
   facet_grid(~type)
 
 save_plot(device="png",prefix = "sim_ct",width = 210,height=150)
 
 

 