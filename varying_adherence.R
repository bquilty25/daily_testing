# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("parameters.R")
source("lft_curves.R")


results_name <- "results_df"

most_recent_file <- file.info(list.files("results/", full.names = T)) %>% 
  as.data.frame() %>% 
  rownames_to_column()%>% 
  filter(str_detect(rowname,"_adherence.fst")) %>% 
  slice_max(mtime) %>% 
  pull(rowname)

assign(results_name,read.fst(most_recent_file))

adherence_list <- get(results_name) %>% 
  filter(#adherence_iso==1,
         #adherence_quar==1,
         #assay=="Innova"|is.na(assay),
         quar_dur==10|is.na(quar_dur),
         test_to_tracing==1
         ) %>% 
  select(-c(adherence_quar,adherence_iso)) %>% 
  crossing(adherence_iso=seq(0,1,by=0.1),adherence_quar=seq(0,1,by=0.1)) %>% 
  mutate(adhering_quar=rbinom(n=n(),size = 1,prob = adherence_quar),
         adhering_iso=rbinom(n=n(),size = 1,prob = adherence_iso)) %>% 
  group_by(adherence_quar,adherence_iso,test_to_tracing) %>%
  group_split() 

plan(multisession, workers = 7)
adherence_results <- adherence_list %>% 
  future_map(.x=.,.f=calc_overlap,.progress = T)

adherence_results_df <- bind_rows(adherence_results)
write.fst(adherence_results_df,paste0("results/results_",st,"_adherence.fst"))

assign(results_name,read.fst(most_recent_file))

#x axis adherence_quar, y axis adherence quar, fill = RR vs 14 day quarantine
baseline <- get(results_name) %>%
    filter(quar_dur  == 10, 
           is.na(assay)) %>%
    filter(!(is.infinite(inf_start) | is.infinite(inf_end))) %>% 
    select(ind_idx, max_overlap, inf_end, inf_start,adherence_quar,adherence_iso,test_to_tracing) %>% 
    group_by(ind_idx,adherence_quar,adherence_iso) %>% 
    summarise(all=sum(inf_end-inf_start),
              baseline_prop=sum(max_overlap)/all) 

x_ <- get(results_name) %>%
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  group_by(stringency,quar_dur,n_tests,assay,ind_idx,adherence_quar,adherence_iso,test_to_tracing) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all) %>% 
  left_join(baseline,by=c("ind_idx","adherence_quar","adherence_iso")) %>% 
   mutate(prop_ratio=prop/baseline_prop) %>% 
   replace_na(list(prop_ratio=0)) 
  

x_ %>%  
  filter(assay=="Innova",!is.na(n_tests),n_tests==5) %>% 
  group_by(stringency,quar_dur,n_tests,assay,adherence_quar,adherence_iso,test_to_tracing) %>% 
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile( .$prop_ratio,
                                              probs = probs)),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  dplyr::select(-data) %>% 
  ggplot(aes(x=adherence_quar,y=adherence_iso,fill=`50%`,label=round(`50%`,2)))+
  geom_raster()+
  geom_text(size=2.5)+
  scale_fill_gradient2(midpoint=1,high=muted("red"),low=muted("blue"))+
  scale_x_continuous(breaks=breaks_width(0.1))+
  scale_y_continuous(breaks=breaks_width(0.1))+
    facet_grid(~n_tests,labeller=labeller(n_tests=function(x)paste(x,"days of testing")))+
    coord_equal(expand = F)+
  theme_minimal()+
  theme(#axis.text = element_text(),axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        legend.position = "bottom",
        axis.text.x=element_text(angle=90))+
  labs(x="Proportion adhering to quarantine",y="Proportion that self-isolate\nafter a positive LFT in DCT",fill="Transmission averted by DCT vs.\n10 days of quarantine (median)")

ggsave("results/risk_surface_5.png",height=120,width=210,dpi=400,units="mm")

x_ %>%  
  filter(!is.na(n_tests)) %>% 
  #mutate(assay=fct_relevel(assay,"Innova (-2.5 CT)", "Innova", "Innova (+2.5 CT)")) %>% 
  group_by(stringency,quar_dur,n_tests,assay,adherence_quar,adherence_iso,test_to_tracing) %>% 
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile( .$prop_ratio,
                                              probs = probs)),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  dplyr::select(-data) %>% 
  ggplot(aes(x=adherence_quar,y=adherence_iso,fill=`50%`,label=round(`50%`,2)))+
  geom_raster()+
  #geom_text(size=2.5)+
  scale_fill_gradient2(midpoint=1,high=muted("red"),low=muted("blue"))+
  scale_x_continuous(breaks=breaks_width(0.1))+
  scale_y_continuous(breaks=breaks_width(0.1))+
  facet_grid(assay~n_tests,labeller=labeller(n_tests=function(x)paste(x,"days of testing")))+
  coord_equal(expand = F)+
  theme_minimal()+
  theme(#axis.text = element_text(),axis.ticks = element_line(),
    panel.border = element_rect(fill=NA),
    legend.position = "bottom",
    axis.text.x=element_text(angle=90))+
  labs(x="Proportion adhering to quarantine",y="Proportion that self-isolate\nafter a positive LFT in DCT",fill="Transmission averted by DCT vs.\n10 days of quarantine (median)")

ggsave("results/risk_surface_all_pres.png",height=150,width=200,dpi=300,units="mm")


#increase in adherence 

lm_fit <- function(x){
  FIT <- lm(prop_ratio ~ adherence_quar+adherence_iso, data = x)
  Pred <- predict(FIT, newdata = crossing(adherence_quar =seq(0, 1, by = 0.01),adherence_iso =seq(0, 1, by = 0.01)))
  crossing(adherence_quar =seq(0, 1, by = 0.01),adherence_iso =seq(0, 1, by = 0.01), pred = Pred)
}

DF2 <- x_ %>% 
  filter(!is.na(assay)) %>% 
  group_by(assay,n_tests,ind_idx) %>% 
  nest() %>% 
  ungroup() %>% 
  sample_n(100) %>% 
  unnest() %>% 
  group_by(assay,ind_idx,n_tests) %>% 
  nest() %>% 
  mutate(pred_ratio = map(data,lm_fit)) %>% 
  select(-data) %>% 
  unnest()

#show predicted line and actual points
ggplot() + geom_line(aes(adherence_iso, pred_ratio,group=adherence_quar), data = DF2) + 
  geom_point(aes(year, lifeExp), data = DF) + facet_wrap(~country) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

x_ %>%
  filter(!is.na(n_tests),adherence_quar==0.5) %>%  
  group_by(ind_idx,stringency,quar_dur,assay,adherence_quar,n_tests,test_to_tracing) %>%
  filter(prop_ratio>=1) %>% 
  filter(adherence_iso==min(adherence_iso)) %>% 
  group_by(stringency,quar_dur,n_tests,assay,adherence_quar,test_to_tracing) %>% 
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile(.$adherence_iso,
                                              probs = c(0.5,0.25,0.75)))) %>%
  unnest_wider(Q) %>% 
  mutate(estimate= paste0(round((adherence_quar-`50%`),2)*100,
                                "% (",
                                round((adherence_quar-`75%`),2)*100,
                                "%, ",
                                round((adherence_quar-`25%`),2)*100,
                                "%)")) %>% 
  select(-c(data,`50%`,`25%`,`75%`,test_to_tracing,adherence_quar)) %>% 
  arrange(assay) %>% 
  htmlTable()


x_ %>%
  filter(!is.na(n_tests),assay=="Innova"|is.na(assay)) %>% 
  group_by(stringency,quar_dur,n_tests,assay,adherence_iso,adherence_quar,test_to_tracing) %>% 
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile( .$prop_ratio,
                                              probs = c(0.5,0.025,0.975)))) %>%
  unnest_wider(Q) %>%  
  ggplot(aes(x=adherence_quar-adherence_iso,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_vline(aes(xintercept=0),linetype="dashed")+
  #coord_flip()+
  geom_ribbon(alpha=0.2)+
  geom_line()+
  scale_x_continuous(labels = scales::percent)+
  #scale_y_continuous(labels = scales::percent)+
  facet_grid(test_to_tracing~n_tests,
             labeller=labeller(n_tests=function(x)paste(x,"days of testing"),
                               test_to_tracing=function(x)paste(x,"days since\nindex case's positive test")))+
  labs(x="Difference in adherence",y="Transmission averted by DCT vs 10 days of quarantine")
             