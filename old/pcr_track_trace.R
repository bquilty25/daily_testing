
# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("parameters.R")
source("kucirka_fitting.R")

# Filter input scenarios if required
input %<>% filter(
  index_test_delay %in% c(1),
  #delay_scaling    == 1,
  waning=="adhere_100",
  #quar_dur         %in% seq(0,14,by=2),
  #stringency       == "Test upon tracing\nand end of quarantine"
) %>% 
  mutate(scenario=row_number()) 

input_split <-
  input %>% 
  rowwise %>%
  group_split

# Name results and create directories
results_name <- "daily vs single"

if (!dir.exists(here::here("results", results_name))){
  dir.create(here::here("results", results_name))
}

con <- file(here::here("results", results_name, "results.log"))
#sink(con, append=FALSE)
#sink(con, append=TRUE, type="message")

# Run analysis
assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f = ~run_analysis(
           n_sims             = 10,
           n_ind_cases        = 1000,
           n_sec_cases        =  10,
           input              = .x,
           seed               =  145,
           P_r                = P_r,
           P_c                = P_c,
           P_t                = P_t,
           dat_gam            = dat_gam,
           asymp_parms        = asymp_fraction,
           return_full        = T
         )))

#sink() 
#sink(type="message")

# Save results
saveRDS(get(results_name), 
        here::here("results", results_name, "results.rds"))
saveRDS(input,
        here::here("results", results_name, "input.rds"))

# Save data to .csv if required
#results <- read_results(results_name)
#write.csv(results,here::here("results", results_name, "summarised.csv"))

