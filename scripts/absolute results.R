# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("kucirka_fitting.R")
source("parameters.R")

results_name <- "results_df"

most_recent_file <- file.info(list.files("results/", full.names = T)) %>% 
  as.data.frame() %>% 
  rownames_to_column()%>% 
  filter(str_detect(rowname,"_all.fst")) %>% 
  slice_max(mtime) %>% 
  pull(rowname)

assign(results_name,read.fst(most_recent_file))

col_pal <- RColorBrewer::brewer.pal(n=4,name = "Dark2")

#Figure 2
plot_a <- get(results_name)%>% 
  filter(adherence_iso==0.67,
         adherence_quar==0.5,
         delay_scaling==1,
        sens_LFA=="higher"|is.na(sens_LFA),
         !multiple_tests
  ) %>%
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Transmission potential averted")+
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_b <-get(results_name) %>% 
  filter(adherence_iso==0.67,
         adherence_quar==0.5,
         delay_scaling==1,
        sens_LFA=="higher"|is.na(sens_LFA),
         multiple_tests) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Transmission potential averted")+
  scale_colour_manual(name="",values = col_pal[4])+
  plotting_theme

plot_a+plot_b+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 400, 
          device = "png",
          prefix = "FigS1",
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
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_quar,adherence_iso,assay,quar_dur,n_tests,delay_scaling,sampling_freq,sens_LFA) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_quar,adherence_iso,assay,quar_dur,n_tests,delay_scaling,sampling_freq,sens_LFA) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.5,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
  select(delay_scaling,assay,adherence_iso,stringency,sens_LFA,quar_dur,`50%`,ui) %>% 
  htmlTable()


### Sensitivity analyses ----

# Delays
plot_a_delays <- get(results_name) %>% 
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
   sens_LFA=="higher"|is.na(sens_LFA),
   !multiple_tests
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Transmission potential averted")+
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

plot_b_delays <- get(results_name) %>% 
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
   sens_LFA=="higher"|is.na(sens_LFA),
    multiple_tests
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,assay,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,assay,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Transmission potential averted")+
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

plot_a_delays+plot_b_delays+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 400, 
          device = "png",
          prefix = "FigS2",
          base = "plot", 
          width = 350, 
          height = 150)

get(results_name) %>% 
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
    sens_LFA=="higher"|is.na(sens_LFA),
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
  select(delay_scaling,assay,stringency,quar_dur,`50%`,ui) %>% 
  htmlTable()

# Adherence
plot_a_adherence <- get(results_name) %>% 
  filter(
    delay_scaling==1,
   sens_LFA=="lower"|is.na(sens_LFA),
    !multiple_tests
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_iso,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_iso,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Transmission potential averted")+
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

plot_b_adherence <- get(results_name) %>% 
  filter(
    delay_scaling==1,
   sens_LFA=="lower"|is.na(sens_LFA),
    multiple_tests
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_iso,assay,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_iso,assay,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Transmission potential averted")+
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

plot_a_adherence+plot_b_adherence+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 400, 
          device = "png",
          prefix = "FigS3",
          base = "plot", 
          width = 300, 
          height = 200)

#In text
get(results_name) %>% 
  filter(
    adherence_quar!=0,
    adherence_iso!=0,
    delay_scaling==1, sens_LFA=="higher"|is.na(sens_LFA)
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_iso,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_iso,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(ui=paste0("(",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
  select(delay_scaling,assay,stringency,quar_dur,`50%`,ui) %>% 
  htmlTable()

### LFA sensitivity analysis
plot_a_lfa <- get(results_name) %>% 
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    !is.na(sens_LFA),
    !multiple_tests
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_iso,adherence_quar,sens_LFA,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_iso,adherence_quar,sens_LFA,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
               ~sens_LFA, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 sens_LFA = sens_scaling_labels,
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
  scale_colour_manual(name="Strategy",values = col_pal[2])

plot_b_lfa <- get(results_name) %>% 
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    assay=="LFA",
    multiple_tests
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_iso,assay,sens_LFA,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_iso,assay,n_tests,sens_LFA,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=stringency),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Transmission potential averted")+
  scale_colour_manual(name="",values = col_pal[4])+
  facet_nested(nest_line=T,
               ~sens_LFA, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 sens_LFA = sens_scaling_labels,
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

plot_a_lfa+plot_b_lfa+plot_annotation(tag_levels = "A")+plot_layout(guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 400, 
          device = "png",
          prefix = "LFA_sens",
          base = "plot", 
          width = 300, 
          height = 150)

get(results_name) %>% 
  filter(
    adherence_quar==0.5,
    adherence_iso==0.67,
    delay_scaling==1, 
    assay=="LFA"
  ) %>%
  mutate(stringency=case_when(multiple_tests&tests~"Daily LFA testing",
                              tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                              tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                              !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,stringency,adherence_iso,adherence_quar,sens_LFA,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(stringency,adherence_iso,adherence_quar,sens_LFA,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
  select(delay_scaling,assay,stringency,sens_LFA,quar_dur,`50%`,ui) %>% 
  htmlTable()

