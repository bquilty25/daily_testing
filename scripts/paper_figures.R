## curves for paper

set.seed(145)
trajectories <- make_trajectories(n_cases = 1000)

trajectories_to_plot <- trajectories$models %>%
  #filter(type == "symptomatic") %>%
  mutate(pred = map(.x = m, 
                    ~data.frame(x = seq(0, 25, length.out = 101)) %>%
                      mutate(y = .x(x)))) %>%
  unnest(pred) %>%
  mutate(LFA= boot::inv.logit(
      predict(innova_liv_mod,
              newdata =
                data.frame(ct = y),
              type = "response")
    )) %>% 
  # mutate(LFA = base::cut(y, 
  #                        breaks = c(-Inf,20,25,30,35,Inf),
  #                        labels=c("82.4%","54.5%","8.3%","5.3%","0%"))) %>%
  #mutate(PCR = percent(as.integer(y < 35))) %>%
  gather(key, value, type) %>% 
  mutate(value=fct_relevel(value,"symptomatic"))#%>% 
  #mutate(value=fct_rev(fct_relevel(value,"100%","82.4%","54.5%","8.3%","5.3%","0%")))

mask <- data.frame(min = c(-Inf, 20, 25, 30, 35),
                   max = c(20, 25, 30, 35, Inf),
                   val = c(.824, .545, .083,0.053,0),
                   key = "LFA") %>%
  bind_rows(data.frame(min = c(-Inf, 35),
                       max = c(35, Inf),
                       val = c(1, 0),
                       key = "PCR"))

trajectories_to_plot %>% 
  ggplot(data = ., aes(x = x, y = y)) +
  # annotate(geom = "rect",
  #          ymax = Inf, xmin = -Inf, xmax = Inf, ymin = 30,
  #          color = NA,
  #          fill = "black", alpha = 0.1) +
  geom_line(aes(group = idx,
                color = LFA,
                ),
            alpha = 0.1) +
  theme_bw() +
  facet_wrap(~value,labeller = labeller(value=capitalize),ncol=1)+
  #geom_hline(aes(yintercept=35,linetype="PCR detection threshold"))+
  geom_hline(aes(yintercept=30,linetype="Infectivity threshold"))+
  scale_linetype_manual(name="",values=c("dotted","dashed"))+
  # geom_hline(data = mask,
  #            aes(yintercept = max),
  #            lty = 2) +
   scale_color_viridis_c(option="magma",begin=0.2,end=0.8,
                      name   = "Probability of detection by LFT",guide=guide_colorsteps(ticks=T,barwidth=unit(5,"cm"))) +
  ylab("Ct value") +
  xlab("Time since exposure (days)") +
  plotting_theme +
  theme(panel.spacing = unit(1,"cm")) +
  scale_y_reverse()+coord_cartesian(expand = F)

  #guides(color = guide_legend(override.aes = list(size=2,alpha = 1))) 
  # geom_text(x = 25, y = 29,
  #           hjust = 1,
  #           label = "Infectious") +
  # geom_text(x = 25, y = 31,
  #           hjust = 1, 
  #           label = "Non-infectious")

save_plot(dpi=400,
          device="png",prefix = "trajectories", height = 200,width=180)
