# tables and in text values


n_delay_sim <- 1e5
n_sim <- 100000

delay_table <- list(
  Testing  = P_r,
  Sourcing = P_c,
  Tracing  = P_t) %>%
  map_df(~data.frame(sim = 1:n_sim,
                     value = time_to_event(n    = n_delay_sim,
                                           mean = .x[["mean"]],
                                           var  = .x[["var"]])),
         .id = "Event") %>%
  bind_rows(., group_by(., sim) %>% 
              summarise(value = sum(value)) %>%
              mutate(Event = "Total"))


delay_summary <- delay_table %>%
  group_by(Event) %>%
  #mutate(sim = 1:n_sim) %>%
  nest(data = -Event) %>%
  mutate(summaries = map(data, ~list(M = mean(.x$value),
                                     Q = quantile(.x$value, probs),
                                     S = sd(.x$value)))) %>%
  unnest_wider(summaries) %>%
  unnest_wider(Q) %>%
  select(Mean = M, SD = S, contains("%")) %>%
  mutate(Event = fct_inorder(Event))

write_csv(x = delay_summary, path = paste0("results/",
                                           results_name, "/Table_1.csv"))

results_infectivity_df <-
  results_infectivity %>% 
  map(show_results, reduction = FALSE) %>%
  map(bind_rows) %>%
  map(~mutate(.x, time_in_iso = factor(time_in_iso,
                                       levels = seq(0, 15),
                                       ordered = T))) %>%
  map_df(bind_rows, .id = "Measure")

# transmission potential ahead of tracing
results_infectivity_df %>%
  filter(stringency == "low", 
         index_test_delay == 2,
         screening == FALSE,
         Measure == "pre") %>%
  select(one_of(all.vars(faceting)),
         contains("%")) %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~percent(round(., 2)))

# how much quarantine time can we shave off if we improve testing delays?
c(2,3) %>%
  set_names(., paste(., "days until index case's test")) %>%
  map(~filter(results_infectivity_df, Measure == "averted") %>%
        filter((stringency == "maximum" & screening == FALSE & index_test_delay == .x) |
                 (stringency == "moderate" & index_test_delay == .x - 1)) %>%
        select(-Measure, -delays, -screening) %>%
        mutate_at(.vars = vars(contains("%")),
                  .funs = ~percent(round(., 2))))

# if we can't reduce delays, can we do double testing?
c(1, 2, 3) %>%
  set_names(., paste(., "days until index case's test")) %>%
  map(~filter(results_infectivity_df, Measure == "averted" & index_test_delay == .x) %>%
        filter((stringency == "maximum" & screening == FALSE) |
                 (stringency == "moderate" )) %>%
        select(-Measure, -screening) %>%
        mutate_at(.vars = vars(contains("%")),
                  .funs = ~percent(round(., 2))))

# baseline at Maximum
results_infectivity_df %>%
  filter(index_test_delay == 2,
         stringency == "maximum",
         Measure == "averted") 

results_infectivity_df %>%
  filter(screening == FALSE,
         stringency == "maximum",
         Measure == "averted")  %>%
  select(index_test_delay, time_in_iso, contains("%"))

# Table 2: max and 10 days (with one or two tests)
results_infectivity_df %>%
  filter(#screening == TRUE,
    time_in_iso %in% c(10,14),
    Measure == "averted")  %>%
  select(stringency, index_test_delay,delays, time_in_iso, contains("%"))  %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~round(.,2)) %>%
  select(-time_in_iso)

# baseline with test on day 0
results_infectivity_df %>%
  filter(index_test_delay == 2,
         stringency == "low",
         Measure == "averted") %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~round(.,2)) %>%
  select(-time_in_iso)

# effect of a test in low
results_infectivity$averted %>%
  filter(stringency == "moderate") %>%
  show_results(reduction = TRUE)

# effect of asymptomatics

results_infectivity_type <- 
  results_df %>%
  make_days_plots(input, 
                  faceting = index_test_delay + delay_scaling ~ stringency + type,
                  y_labels = infectivity_labels["infectivity_averted"],
                  base     = paste0(results_name,"_averted_type_"),
                  sum      = F)

# Table 3... Table 2 + type
results_infectivity_type$averted %>%
  split(.$type) %>%
  map(~filter(.x,time_in_iso == 14 | time_in_iso == 1 | 
                (time_in_iso == 10 & stringency == "moderate")) %>%
        show_results(reduction = F)) %>%
  map_df(bind_rows, .id = "type") %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~round(100 * .)) %>%
  arrange(desc(time_in_iso), index_test_delay, desc(type)) %>%
  mutate(`PCR tests` = 
           ifelse(!screening, 
                  "no test",
                  paste("day", sub(pattern = " & NA", 
                                   replacement = "", 
                                   x = delays)))) %>%
  mutate(`PCR tests` = sprintf("%s (%s)", time_in_iso, `PCR tests`)) %>%
  select(-screening, -delays, -time_in_iso) %>%
  group_by_at(.vars = vars(-contains("%"))) %>%
  transmute(Median = sprintf("%.f%%", `50%`),
            `50% CI (IQR)`  = sprintf("(%.f%%, %.f%%)", `25%`, `75%`),
            `95% CI`  = sprintf("(%.f%%, %.f%%)", `2.5%`, `97.5%`)) %>%
  write_csv(., paste0("results/",results_name,"_Table_3.csv"))