# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("kucirka_fitting.R")
source("parameters.R")

results_name <- "results_list"
assign(results_name,read.fst("results_2020427_all.fst"))

baseline_high <- get(results_name) %>% 
  filter(adherence_iso==0.67,adherence_quar==0.5,delay_scaling==1,quar_dur==14,!tests,is.na(sens_LFA)) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  select(ind_idx, max_overlap,symp_overlap,inf_end,inf_start) %>% 
  group_by(ind_idx) %>% 
  summarise(baseline_prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  mutate(name="baseline_high")

baseline_low <- get(results_name) %>% 
  filter(adherence_iso==0.67,adherence_quar==0.5,delay_scaling==1,quar_dur==0,!tests,is.na(sens_LFA)) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  select(ind_idx, max_overlap,symp_overlap,inf_end,inf_start) %>% 
  group_by(ind_idx) %>% 
  summarise(baseline_prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  mutate(name="baseline_low")

# Figure 2
plot_1 <- 
  get(results_name) %>% 
  filter(
    sens_LFA=="higher"|is.na(sens_LFA),
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    !multiple_tests
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,sampling_freq) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,assay,quar_dur,n_tests,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(label=ifelse(quar_dur==14&stringency=="Post-exposure quarantine only","Baseline",NA)) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_text_repel(aes(label=label,colour=stringency),
                  position=position_dodge(width=0.5),
                  direction    = "y",
                  angle        = 90,
                  segment.size = 0.2,
                  show.legend = F)+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=2,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_log10(limits=c(0.3,2), 
    #breaks = logTicks(n = 10), minor_breaks = logTicks(n = 40)
    )+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed T&T delays")+
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_2 <- get(results_name) %>% 
  filter(
    sens_LFA=="higher"|is.na(sens_LFA),
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    multiple_tests
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,sampling_freq,sens_LFA) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,assay,quar_dur,n_tests,sampling_freq,sens_LFA) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             size=1.5,
             position=position_dodge(width=0.5)) +
  scale_y_log10(limits=c(0.3,2), 
    )+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed T&T delays")+
  scale_colour_manual(name="",values = col_pal[4])+
  plotting_theme

plot_1+plot_2+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 500, 
          device = "png",
          prefix = "Fig2",
          base = "plot", 
          width = 300, 
          height = 150)

#In text
get(results_name) %>% 
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    sens_LFA=="higher"|is.na(sens_LFA)
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,sampling_freq) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,assay,quar_dur,n_tests,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(label=ifelse(quar_dur==14&stringency=="Post-exposure quarantine only","Reference",NA)) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = txtRound,digits=2) %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
  select(assay,stringency,quar_dur,`50%`,ui) %>% 
  htmlTable()

plot_1_delays <- 
  get(results_name) %>%
  filter( sens_LFA=="higher"|is.na(sens_LFA),
  adherence_iso==0.67,
  adherence_quar==0.5,
  !multiple_tests
) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(label=ifelse(delay_scaling==1&quar_dur==14&stringency=="Post-exposure quarantine only","Baseline",NA)) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_text_repel(aes(label=label,colour=stringency),
                  position=position_dodge(width=0.5),
                  direction    = "y",
                  angle        = 90,
                  segment.size = 0.2,
                  show.legend = F)+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             size=2,
             position=position_dodge(width=0.5)) +
  scale_y_log10(limits=c(NA,2.5)
                )+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed T&T delays")+
  facet_nested(nest_line=T,
               ~delay_scaling, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_2_delays <- get(results_name) %>%
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
    sens_LFA=="higher"|is.na(sens_LFA),
    multiple_tests
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             size=1.5,
             position=position_dodge(width=0.5)) +
  scale_y_log10(limits=c(NA,2.5), 
                )+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed T&T delays")+
  scale_colour_manual(name="",values = col_pal[4])+
  facet_nested(nest_line=T,
               ~delay_scaling, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme

plot_1_delays/plot_2_delays+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 300, 
          device = "png",
          prefix = "Fig3",
          base = "plot", 
          width = 300, 
          height = 250)

get(results_name) %>% 
  filter(sens_LFA=="higher"|is.na(sens_LFA),
    adherence_iso==0.67,
    adherence_quar==0.5,
    
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,sampling_freq,delay_scaling) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,assay,quar_dur,n_tests,sampling_freq,delay_scaling) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(label=ifelse(quar_dur==14&stringency=="Post-exposure quarantine only","Reference",NA)) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = txtRound,digits=2) %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
  select(assay,stringency,delay_scaling,quar_dur,`50%`,ui) %>% 
  htmlTable()

# Figure 4
plot_1_adherence <- get(results_name) %>% 
  filter(
    adherence_iso!=0,
    adherence_quar!=0,
    delay_scaling==1,
    sens_LFA=="higher"|is.na(sens_LFA),
    !multiple_tests
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,adherence_iso,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,adherence_iso,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(label=ifelse(adherence_quar==0.5&adherence_iso==0.67&delay_scaling==1&quar_dur==14&stringency=="Post-exposure quarantine only","Baseline",NA)) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_text_repel(aes(label=label,colour=stringency),
                  position=position_dodge(width=0.5),
                  direction    = "y",
                  angle        = 90,
                  segment.size = 0.2,
                  show.legend = F)+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             size=1,
             position=position_dodge(width=0.5)) +
  scale_y_log10(limits=c(0.3,NA))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed T&T delays")+
  facet_nested(nest_line=T,
               adherence_iso~adherence_quar, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence_quar =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine"),
                 adherence_iso =
                   c("1" = "100% adhere\nto isolation",
                     "0.67" =
                       "67% adhere\nto isolation",
                     "0" =
                       "0% adhere\nto isolation")
               )) +
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_2_adherence <-get(results_name) %>% 
  filter(
    adherence_iso!=0,
    adherence_quar==1,
    sens_LFA=="higher"|is.na(sens_LFA),
    delay_scaling==1,
    multiple_tests
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,adherence_iso,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,adherence_iso,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             size=1,
             position=position_dodge(width=0.5)) +
  scale_y_log10(limits=c(0.3,NA))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed T&T delays")+
  scale_colour_manual(name="",values = col_pal[4])+
  facet_nested(nest_line=T, 
               adherence_iso~., labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence_quar =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine"),
                 adherence_iso =
                   c("1" = "100% adhere\nto isolation",
                     "0.67" =
                       "67% adhere\nto isolation",
                     "0" =
                       "0% adhere\nto isolation")
               )) +
  plotting_theme

plot_1_adherence+plot_2_adherence+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 400, 
          device = "png",
          prefix = "Fig4",
          base = "plot", 
          width = 300, 
          height = 150)

get(results_name) %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso!=0,
    adherence_quar==1,sens_LFA=="higher"|is.na(sens_LFA),
    delay_scaling==1,
    multiple_tests
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,adherence_iso,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,adherence_iso,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = txtRound,digits=2) %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
  select(assay,stringency,quar_dur,`50%`,ui) %>% 
  htmlTable()

plot_1_sens <- 
  get(results_name) %>%
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    !multiple_tests,
    assay=="LFA"
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq,sens_LFA) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq,sens_LFA) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(label=ifelse(delay_scaling==1&quar_dur==14&stringency=="Post-exposure quarantine only","Baseline",NA)) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=2,
             position=position_dodge(width=0.5)) +
  scale_y_log10(
                breaks = logTicks(n = 5), minor_breaks = logTicks(n = 40))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed T&T delays")+
  facet_nested(nest_line=T,
               ~sens_LFA, labeller = labeller(
                 type = capitalize,
                 sens_LFA=sens_scaling_labels,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[2])

plot_2_sens <- get(results_name) %>%
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    multiple_tests,
    n_tests!=14
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq,sens_LFA) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq,sens_LFA) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_log10(breaks = logTicks(n = 10), minor_breaks = logTicks(n = 80))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed T&T delays")+
  facet_nested(nest_line=T,
               ~sens_LFA, labeller = labeller(
                 type = capitalize,
                 sens_LFA=sens_scaling_labels,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme+
  scale_colour_manual(name="",values = col_pal[4])

fig_sens <- plot_1_sens+plot_2_sens

p_ranges_y <- c(10^(ggplot_build(fig_sens[[1]])$layout$panel_scales_y[[1]]$range$range),
                10^(ggplot_build(fig_sens[[2]])$layout$panel_scales_y[[1]]$range$range))

fig_sens+plot_annotation(tag_levels = "A")+plot_layout(guides = "collect")&theme(legend.position = "bottom")&ylim(min(p_ranges_y), max(p_ranges_y))

save_plot(dpi = 300, 
          device = "png",
          prefix = "FigS4",
          base = "plot", 
          width = 300, 
          height = 150)

get(results_name) %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
  ) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(ind_idx,stringency,adherence_iso,assay,sens_LFA,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(prop=sum(max_overlap)/sum(inf_end-inf_start)) %>% 
  left_join(baseline_high) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  replace_na(list(prop_ratio=1)) %>% 
  group_by(stringency,adherence_iso,assay,sens_LFA,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = txtRound,digits=2) %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
  select(assay,stringency,quar_dur,`50%`,ui) %>% 
  htmlTable()
