
# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("kucirka_fitting.R")
source("parameters.R")

run_model <- function(
  n_sims          = 1000,
  n_sec_cases     = 1000, # this shouldn't matter. just needs to be Big Enough
  n_ind_cases     = 10000,
  input,
  trajectories,
  seed            = 145,
  asymp_parms
){
  
  
  conflicted::conflict_prefer("set_names", "purrr")
  conflicted::conflict_prefer("melt", "reshape2")
  map(.x = c("mutate", "select", "filter"), 
      .f = function(x){conflicted::conflict_prefer(name = x, "dplyr")})
  
  #browser()
  
  set.seed(seed)
  
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), input$scenario))
  
  traj <- trajectories$traj %>% 
    select(-y) %>% 
    pivot_wider(names_from = name,values_from=x) %>% 
    select(-c(start,end))
  
  inf <- data.frame(prop_asy = rbeta(n = n_sims,
                                   shape1 = asymp_parms$shape1,
                                   shape2 = asymp_parms$shape2)) 



# Generate index cases' inc times
index_cases <- traj %>% 
  left_join(trajectories$models,by=c("idx","type")) %>% 
  filter(type=="symptomatic",!is.infinite(inf_start),!is.infinite(inf_end)) %>% 
  select(-c(data,m)) %>% 
  sample_n(n_sims) %>% 
  rename(ind_idx = idx) %>% 
  bind_cols(inf) %>% 
  #sample test result delay
  ## sample uniformly between 0 and 1 when 0.5...
  crossing(distinct(input, index_test_delay, delay_scaling,test_to_tracing)) %>%     
  rename("index_onset_t"  = onset_t) %>% 
  mutate(index_testing_t  = index_onset_t + index_test_delay,
         sec_traced_t     = index_onset_t + index_test_delay + test_to_tracing) %>% 
  mutate(sec_exposed_t = runif(n=n(),min = inf_start,max=pmin(inf_end,index_testing_t)))




my_message("Generating secondary cases' trajectories and asymptomatic fraction")
sec_cases <- index_cases %>% 
  nest(data = -c(ind_idx,
                 prop_asy,
                 test_to_tracing,
                 index_onset_t,
                 index_test_delay,
                 index_testing_t,
                 sec_traced_t,
                 sec_exposed_t,
                 delay_scaling)) %>% 
  mutate(prop_asy    = as.list(prop_asy)) %>%
  mutate(sec_cases   = map(.x = prop_asy, 
                           .f  = ~make_sec_cases(as.numeric(.x),
                                                 traj,
                                                 n_sec_cases)
  ))

sec_cases %<>%
  unnest(prop_asy) %>%
  unnest(sec_cases) %>% 
  ungroup()%>% rename_at(.vars = vars(idx,onset_t),
                       .funs = ~paste0("sec_", .)) %>%
  dplyr::select(-data)


my_message("Shifting secondary cases' times relative to exposure to index")
#exposure date relative to index cases exposure
 sec_cases %<>% 
  mutate_at(.vars = vars(sec_onset_t),
            .funs = function(x,y){x + y}, y = .$sec_exposed_t) 


sec_cases <- left_join(input,
                       sec_cases,
                       by = c("index_test_delay", "delay_scaling","test_to_tracing")) %>% 
  mutate(adhering_quar=rbinom(n=n(),size = 1,prob = adherence_quar),
         adhering_iso=rbinom(n=n(),size = 1,prob = adherence_iso)) 


# generate testing times
my_message("Calculating test times")
sec_cases %<>% 
  mutate(test_t = pmap(.l = list(multiple_tests=multiple_tests,
                                 tests=tests,
                                 tracing_t=sec_traced_t,
                                 sampling_freq=sampling_freq,
                                 sec_exposed_t=sec_exposed_t,
                                 quar_dur=quar_dur,
                                 max_tests=14,
                                 n_tests=n_tests),
                      .f = test_times)) %>% 
  unnest(test_t) 


#calc outcomes 
my_message("Calculating outcomes for each secondary case")
sec_cases %<>% calc_outcomes(x  = .) %>% select(-c(m,rx,ry,u.x,u.y))

#find earliest positive test result
sec_cases %<>% 
  nest(test_t, test_q, test_p, test_no, test_label, have_test, ct, detection_range, screen) %>% 
  mutate(earliest_q      =  map(.f = earliest_pos2, 
                                    .x = data)) %>% 
  unnest_wider(earliest_q) %>% 
  rename("earliest_q"=test_q)%>% 
  select(-data)

#shift other timings relative to onset
sec_cases %<>%
  mutate(exposed_q  = sec_exposed_t - sec_exposed_t,
         traced_q   = sec_traced_t - sec_exposed_t, 
         quar_end_q = pmax(traced_q,(exposed_q + quar_dur)),
         onset_q    = sec_onset_t - sec_exposed_t,
         symp_end_q = onset_q + post_symptom_window,
         test_iso_end_q = earliest_q + post_symptom_window
         )
#browser()
# calculate remaining transmission potential averted by positive test
my_message("Calculating remaining transmission potential for each secondary case")
averted <- sec_cases %>% calc_overlap(.)
  
return(averted)

}

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `Daily testing` = 
      crossing(sampling_freq    = 1,
               tests            = TRUE,
               multiple_tests   = TRUE,
               n_tests          = default_testing, 
               assay            = "LFA",
               sens_LFA     = c("higher","lower"),
               quar_dur         = NA), 
    `Post-exposure quarantine only` = 
      crossing(sampling_freq    = NA,
               tests            = FALSE,
               multiple_tests   = FALSE,
               assay            = NA,
               n_tests          = NA,
               quar_dur         = c(0, default_testing[-1])),
    `Post-exposure quarantine with LFA test` = 
      crossing(sampling_freq    = NA,
               tests            = TRUE,
               multiple_tests   = FALSE,
               n_tests          = NA,
               assay            = "LFA",
               sens_LFA     = c("higher","lower"),
               quar_dur         = c(0, default_testing[-1])),
    `Post-exposure quarantine with PCR test` = 
      crossing(sampling_freq    = NA,
               tests            = TRUE,
               multiple_tests   = FALSE,
               n_tests          = NA,
               assay            = "PCR",
               quar_dur         = c(0, default_testing[-1]))
  ) %>% 
    bind_rows(.id = "stringency")) %>% 
  crossing(post_symptom_window = 10,
           index_test_delay    = c(1),  # time to entering quarantine (index cases)
           delay_scaling       = c(1, 0.5, 0),
           adherence_quar      = c(0, 0.5,  1),
           adherence_iso       = c(0, 0.67, 1)) %>% 
  mutate(test_to_tracing       = 3*delay_scaling) %>% 
  filter(#!(delay_scaling!=1&adherence_iso!=0.67&adherence_quar!=0.5)
         #adherence_iso==0.67,adherence_quar==0.5,
        # stringency=="Post-exposure quarantine with LFA test"
         ) %>% 
  mutate(scenario=row_number()) 

input_split <-
  input %>%
  rowwise %>%
  group_split()

trajectories <- make_trajectories(n_cases = 1000)

results_name <- "results_list"


assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f =  ~ run_model(
           input=.x,
           trajectories=trajectories,
           seed = 1000,
           n_sec_cases = 10,
           n_sims = 100,
           asymp_parms = asymp_fraction
         )))

results_df <- get(results_name) %>% 
  bind_rows() %>% 
  as.data.frame() 

st=format(Sys.time(), "%Y%m%d_%H%M%S")
write.fst(get(results_name) %>% 
            bind_rows() %>% 
            select(-c(rx,ry,m)) %>% 
            as.data.frame(),paste0("results/results_",st,"_all.fst"))
