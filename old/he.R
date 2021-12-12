{get.gamma.par(
  q = c(0.8, 3),
  p = c(0.025, 0.975),
  plot = F, show.output = F)} -> delay_parms 

HE <- data.frame(p = seq(1e-5, 1 - 1e-5, length.out = 101)) %>%
  mutate(delay =        
           qgamma(p, 
                  shape = delay_parms[["shape"]], 
                  rate  = delay_parms[["rate"]])) %>%
  add_row(p = 0, delay = 0, .before = 1)
