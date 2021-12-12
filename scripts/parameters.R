# parameters for simulation

# Ashcroft et al. infectivity profile
# https://smw.ch/article/doi/smw.2020.20336
infect_shape = 97.18750 
infect_rate  =  3.71875
infect_shift = 25.62500

# # McAloon et al. incubation period meta-analysis
#https://bmjopen.bmj.com/content/10/8/e039652
inc_parms <- list(mu_inc = 1.63,
                  sigma_inc = 0.5)

gen_shape <- 9.89
gen_rate <- 2.06

pathogen <- list(
  symptomatic = 
    # review paper Byrne et al. (2020) https://doi.org/10.1101/2020.04.25.20079889
    # define T1 as infection to beginning of presymptomatic infectious period
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      inc_parms,
      
      # Li et al https://www.nejm.org/doi/full/10.1056/nejmoa2001316
      # variance calculated by inverting confidence interval
      list(mu_inf    = 9.1,
           sigma_inf = 14.7)),
  
  asymptomatic = 
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      inc_parms,
      # https://doi.org/10.1101/2020.04.25.20079889
      list(
        mu_inf    =  6,
        sigma_inf = 12))) %>%
  map(~data.frame(.x), .id = "type")

# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- rriskDistributions::get.beta.par(
  q = c(0.24,  0.38),
  p = c(0.025, 0.975), 
  show.output = F, plot = F) %>%
  as.list

waning_none <- function(x){
  waning_points(x, X = 0, Y = 1)
}

adhere_10 <- function(x){
  waning_points(x, X = 0, Y = 0.1)
}
adhere_20 <- function(x){
  waning_points(x, X = 0, Y = 0.2)
}
adhere_30 <- function(x){
  waning_points(x, X = 0, Y = 0.3)
}
adhere_40 <- function(x){
  waning_points(x, X = 0, Y = 0.4)
}
adhere_50 <- function(x){
  waning_points(x, X = 0, Y = 0.5)
}
adhere_60 <- function(x){
  waning_points(x, X = 0, Y = 0.6)
}
adhere_70 <- function(x){
  waning_points(x, X = 0, Y = 0.7)
}
adhere_80 <- function(x){
  waning_points(x, X = 0, Y = 0.8)
}
adhere_90 <- function(x){
  waning_points(x, X = 0, Y = 0.9)
}
adhere_100 <- function(x){
  waning_points(x, X = 0, Y = 1)
}

waning_constant <- function(x){
  waning_points(x, X = 0, Y = 0.75)
}

# waning_drop <- function(x){
#   waning_piecewise_linear(x, 0.75, 0.25, 7, 14)
# }

# waning_linear <- function(x){
#   waning_piecewise_linear(x, ymax = 0.75, .16, 0, 8.3)
# }

# waning_canada_community <- function(x){
#   waning_points(x, X = c(0, 30), Y = c(1, 0.541), log = T)
# }

waning_canada_total <- function(x){
  waning_points(x, X = c(0, 30), Y = c(1, 0.158), log = T)
}

smith_uk <- function(x){
  waning_points(x, X = 0, Y = 0.109)
}

default_testing <- c(1, 3, 5, 7, 10, 14)

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `No testing` = 
      data.frame(sampling_freq=NA,
                 tests=F,
                 multiple_tests=NA,
                 assay=NA,
                 quar_dur=10),
    `Daily testing` = 
      data.frame(sampling_freq=1,
                 tests=T,
                 multiple_tests=T,
                 assay="LFA",
                 quar_dur=NA), 
    `Test at end of quarantine only` = 
      data.frame(sampling_freq=NA,
                 multiple_tests=F,
                 assay="LFA",
                 quar_dur = 10)) %>% 
      bind_rows(.id = "stringency")) %>% 
  crossing(post_symptom_window =  10,
           index_test_delay    =  c(1, 2, 3),  # time to entering quarantine
           delay_scaling       =  c(1, 0.5, 0),
           waning              =c("adhere_100"),
           adherence=c(0,0.5,1)) %>% 
  mutate(test_to_tracing=4*delay_scaling) %>% 
  mutate(scenario=row_number()) %>% 
  #calculate time until release from exposure for each scenario
  mutate(time_since_exp=quar_dur)


## Matching PCR and LFA curves -----

standardise <- function(x){
  (x - mean(x))/sd(x)
}

my_dist <- function(x, y, k = 1){
  # k sets by how much the distance should be influenced by time
  X <- unlist(x[["diff_s"]]  - y[["diff_s"]])
  Y <- unlist(x[["value_r"]] - y[["value_r"]])
  D <- sqrt(k*X^2 + Y^2)
  D_min <- D == min(D)
  # to break ties, sample at random
  sample_n(select(filter(y, D_min), iter_pcr = iter), size = 1)
}


if (!file.exists("data/matched_curves.rds")){
  
  PCR_curves <- read_csv("data/posterior_samples_ct_threshold_37.csv")
  
  #find peak timing
  curves_peak <- PCR_curves %>%
    group_by(iter) %>%
    slice_max(value) %>%
    select(diff) %>%
    rename(peak_timing=diff,
           idx=iter) %>%
    ungroup()
  
  
  LFA_curves <- read_csv("data/posterior_samples_ct_threshold_28.csv")
  
  
  
  curves <- list(LFA = LFA_curves,
                 PCR = PCR_curves) %>%
    map(~group_by(.x, iter) %>%
          filter(value == max(value)) %>%
          ungroup) %>%
    #map(~head(.x, 100)) %>%
    map(~mutate(.x,
                value_s = standardise(value),
                diff_s  = standardise(diff),
                value_r = rank(value),
                diff_r  = rank(diff))) %>%
    {bind_cols(.[[1]],
               bind_rows(lapply(X = group_split(rowwise(.[[1]])),
                                FUN = function(x){
                                  my_dist(x = x,
                                          y =.[[2]],
                                          k = 100)
                                })
               ))
    } %>%
    rename(iter_lfa = iter)  %>%
    select(iter_lfa, iter_pcr) %>%
    tibble::rowid_to_column(.) %>%
    nest(data = -c(rowid, iter_lfa)) %>%
    inner_join(LFA_curves, by = c("iter_lfa" = "iter")) %>%
    unnest(data) %>%
    select(-X1) %>%
    rename(LFA = value) %>%
    left_join(PCR_curves, by = c("iter_pcr" = "iter", "diff")) %>%
    select(-X1) %>%
    rename(PCR = value) %>%
    gather(key, value, LFA, PCR) %>%
    rename("assay"=key,
           "idx"=iter_lfa,
           "days_since_infection"=diff) %>%
    left_join(curves_peak)
  
  saveRDS(curves,"data/matched_curves.rds")
  
} else {
  curves <- read_rds("data/matched_curves.rds")
}
