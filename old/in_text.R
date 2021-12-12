# in text numbers

## Reduction in transmission potential
results_name <- "sum_results"
results <- read_results(results_name)
    
results %>% 
    filter(yvar == "infectivity_pre",
           type == "all")

## longer quarantine reduces post-release 

read_results(results_name) %>%
    filter(yvar == "infectivity_pre",
           type == "all",
           quar_dur %in% c(0, 7, 14)) # select on quar_dur rather than time_since_exp

read_results(results_name) %>%
  filter(
yvar == "infectivity_averted",
#delay_scaling==1,
#adherence==0.5,
type == "all",
quar_dur %in% c(0,7,10,14)) %>% 
    select(stringency,assay, delay_scaling, adherence,type, quar_dur, contains("%")) %>%
    mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
    unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
    unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
    mutate(iqr=paste0("(",iqr,")"),
           ui=paste0("(",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
    select(delay_scaling,assay,stringency,adherence,type,quar_dur,`50%`,iqr,ui) %>% 
    htmlTable()

## changing TTI delays

read_results(results_name) %>% 
    filter(yvar == "infectivity_averted",
           type == "all",
           #delay_scaling == 1,
           #stringency=="none",
           waning=="waning_none",
           index_test_delay == 2,
           quar_dur %in% c(0,7,10,14)
           ) %>%
    select(stringency, delay_scaling, quar_dur, contains("%")) %>%
    mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
    unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
    unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
    mutate(iqr=paste0("(",iqr,")"),
           ui=paste0("(",ui,")")) %>% 
    select(delay_scaling,stringency,quar_dur,`50%`,iqr,ui) %>% 
    arrange(delay_scaling) %>% 
    htmlTable()

# waning
read_results(results_name) %>% 
    filter(yvar == "infectivity_averted",
           type == "all",
           delay_scaling == 1,
           #stringency=="one",
           #waning=="waning_none",
           index_test_delay == 2,
           quar_dur %in% c(0,14)
    ) %>%
    select(stringency, waning, quar_dur, contains("%")) %>%
    mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
    unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
    unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
    mutate(iqr=paste0("(",iqr,")"),
           ui=paste0("(",ui,")")) %>% 
    select(stringency,waning,quar_dur,`50%`,iqr,ui) %>% 
    htmlTable()


# Reduced or waning adherence to quarantine

results_sum <- read_results("sum_results")

# what's the difference between waning rates?
filter(results_sum, 
       index_test_delay == 2,
       delay_scaling    == 1,
       yvar             == "infectivity_averted",
       type             == "all") %>%
  select(time_since_exp, `50%`, stringency, waning) %>%
  spread(key = waning, value = `50%` ) %>%
  mutate(diff = waning_constant - waning_canada_total) %>%
  group_by(stringency) %>%
  filter(diff == min(diff) | time_since_exp %in% range(time_since_exp)) %>%
  mutate(diff = round(diff, 2))

# Reducing index cases test delays

filter(results_sum, 
       #index_test_delay %in% c(1,2),
       delay_scaling    == 1,
       time_since_exp   == 14,
       stringency       == "none",
       yvar             == "infectivity_averted",
       waning           == "waning_none",
       type             == "all") %>%
  select(index_test_delay, contains("%")) %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~percent(x = ., accuracy = 1))

filter(results_sum, 
       index_test_delay %in% c(1,2,3),
       delay_scaling    == 1,
       time_since_exp   == 10,
       stringency       == "none",
       yvar             == "infectivity_averted",
       waning           == "waning_none",
       type             == "all") %>%
  select(index_test_delay, contains("%")) %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~percent(x = ., accuracy = 1))

#adherence sensitivity

read_results(results_name) %>% 
  filter(yvar=="infectivity_averted") %>% 
  mutate(waning=str_sub(waning,-2)) %>% 
  mutate(waning=ifelse(waning=="00",100,waning),waning=as.factor(as.numeric(waning))) %>% View()
  
