## distributions for paper

## total


delay_pred <- data.frame(x = seq(0, 10, by = 0.01))

delay_total <- data.frame(sim = 1:10000) %>%
  mutate(Testing = 
           time_to_event(n = n(),
                         mean = P_r[["mean"]],
                         var  = P_r[["var"]])) %>% 
  #sample contact info delay
  mutate(Sourcing = 
           time_to_event(n = n(),
                         mean = P_c[["mean"]],
                         var  = P_c[["var"]])) %>% 
  #sample tracing delay
  mutate(Tracing      =
           time_to_event(n = n(),
                         mean = P_t[["mean"]],
                         var  = P_t[["var"]])) %>%
  mutate(Total = Testing + Sourcing + Tracing)

total_density <- 
  pull(delay_total, Total) %>%
  {density(log(.))} %>%
  {data.frame(x = exp(.$x),
              y = .$y, Event = "Total")}

index_result_delay <- 
  list(regional  = "Table_3",
       mobile    = "Table_4",
       satellite = "Table_5",
       home      = "Table_6") %>%
  map_df(~read_delays(
    sheet = .x, 
    path = tti_file,
    string = "^Number of test results received .* of taking a test$"),
    .id = "source")



#Delay from positive test to getting info on close contacts from index
getting_contact_info <- 
  read_delays("Table_10", path = tti_file,
              string = "^Number of people reached")


#Delay from getting info to tracing contacts
tracing_delay <-
  read_delays("Table_13", path = tti_file,
              string = "^Number of people reached")

delay_histograms <- 
  list(
    `Testing` = index_result_delay,
    `Sourcing`= getting_contact_info,
    `Tracing` = tracing_delay) %>%
  bind_rows(.id = "Event")  %>%
  mutate(Event = fct_inorder(factor(Event))) %>%
  count(Event, t) %>%
  group_by(Event) %>%
  mutate(y = n/sum(n)) %>%
  rename(x = t)

delay_labeller <- function(x){
  fct_recode(x,
             "Return of index case's test" = "Testing",
             "Sourcing of contacts"        = "Sourcing",
             "Tracing of contacts"         = "Tracing",
             "Total"                       = "Total")
}

delay_distributions <-
  list(
    Testing  = P_r,
    Sourcing = P_c,
    Tracing  = P_t) %>%
  map(~mv2gamma(.x$mean, .x$var)) %>%
  map_df(.x = .,
         ~mutate(delay_pred,
                 y = dgamma(x = x,
                            shape = .x$shape,
                            rate  = .x$rate)),
         .id = "Event") %>%
  bind_rows(total_density) %>%
  mutate(Event = fct_inorder(factor(Event)))

delay_plot <- ggplot(data = delay_distributions, aes(x=x, y=y)) + 
  geom_col(data = delay_histograms,
           fill = lshtm_greens[2],
           color = "white", width = 1,
           alpha=0.7) +
  geom_line(color = lshtm_greens[1]) +
  #geom_segment(data = delay_summary,
               #aes(x = `2.5%`, xend = `97.5%`,
               #    y = 0, yend = 0),
               #color = lshtm_greens[1],
               #size = 1) +
  #geom_point(data = delay_summary,
             # aes(x = `50%`, y = 0),
             # pch = 21,
             # fill = "white",
             # color = lshtm_greens[1],
             # size = 1) +
  facet_wrap( ~Event,
              nrow = 2,
              labeller = labeller(Event = delay_labeller)) +
  theme_minimal() + 
  theme(panel.border=element_rect(fill=NA)) +
  xlab("Time to completion (days)") +
  ylab("Density") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(limits = c(0,12),
                     breaks = breaks_width(2))


save_plot(delay_plot, dpi = 320, device = "png",
          prefix = "delays",
          base = "plot", 
          width = 210, 
          height = 140)

ashcroft_plot <- time_to_event_lnorm(n = 1e4, 
                                     inc_parms$mu_inc,
                                     inc_parms$sigma_inc) %>%
  {data.frame(index_onset_t = .)} %>%
  mutate(sec_exposed_t = index_onset_t - infect_shift + 
           rtgamma(n     = n(),
                   a     = infect_shift,
                   b     = Inf, 
                   shape = infect_shape,
                   rate  = infect_rate) 
  ) %>%
  ggplot(data = ., aes(x = sec_exposed_t)) +
  geom_histogram(binwidth = 1, center = 0.5,
                 aes(y = ..density..),
                 fill = lshtm_greens[2],
                 color = "white",
                 alpha = 0.7)  +
  theme_minimal() + 
  theme(panel.border=element_rect(fill=NA)) +
  xlab("Time from exposure (days)") +
  ylab("Density") +
  #coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(limits = c(0,20),
                     breaks = pretty_breaks(4))

save_plot(ashcroft_plot, dpi = 320, device = "png",
          prefix = "ashcroft",
          base = "plot", 
          width = 210, 
          height = 140)


#waning figs
adherence <- tibble(x=0:30) %>% 
  mutate(Decaying = waning_points(x=x,X=c(0,30),Y=c(1,0.158),log=TRUE),
         Constant = waning_points(x, X = 0, Y = 0.75),
         Perfect  = waning_points(x, X=0, Y=1)) %>% 
  pivot_longer(cols=Decaying:Perfect) %>% 
  ggplot()+
  geom_line(aes(x=x,y=value,colour=fct_rev(name)),size=2)+
  theme_minimal() + 
  theme(panel.border=element_rect(fill=NA))+
  labs(x="Time since tracing (days)",
       y="w(t)")+
  scale_y_continuous(limits=c(0,1))+
  scale_color_brewer(name="Adherence",type="qual")

save_plot(adherence,
          dpi = 320, 
          device = "png",
          prefix = "adherence",
          base = "plot", 
          width = 210, 
          height = 140)
  

iso_fill <- c(rep("#00AEC7",5),rep("#FE5000",5),rep("#FFB81C",10))

time_in_iso_plot <- time_to_event_lnorm(n = 1e4, 
                                     inc_parms$mu_inc,
                                     inc_parms$sigma_inc) %>%
  {data.frame(index_onset_t = .)} %>%
  mutate(sec_exposed_t = index_onset_t - infect_shift + 
           rtgamma(n     = n(),
                   a     = infect_shift,
                   b     = Inf, 
                   shape = infect_shape,
                   rate  = infect_rate) 
  ) %>%
  mutate(iso=case_when(sec_exposed_t<=5~"Prior to tracing",
                       sec_exposed_t>10~"After quarantine",
                       TRUE ~ "During quarantine")) %>% 
  ggplot(data = ., aes(x = sec_exposed_t#,fill=iso
                       )) +
  geom_histogram(binwidth = 1, center = 0.5,
                aes(y = ..density..),
                 #position = "identity",
                 fill = iso_fill,
                 color = "white",
                 alpha = 0.7)  +
  theme_minimal() + 
  theme(panel.border=element_rect(fill=NA)) +
  xlab("Time from exposure (days)") +
  ylab("Density") +
  scale_fill_manual(values=c("Prior to tracing", "During quarantine", "Post-quarantine"),show.legend=T)+
  #coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(limits = c(0,20),
                     breaks = pretty_breaks(4))
time_in_iso_plot
save_plot(time_in_iso_plot, dpi = 320, device = "png",
          prefix = "time_in_iso",
          base = "plot", 
          width = 210, 
          height = 140)
