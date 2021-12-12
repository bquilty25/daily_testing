my_message <- function(x, ...){
  message(paste(Sys.time(), x, sep = "    "), ...)
}

# is this not common to many scripts?
main_scenarios <-
  list(`low` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation")),
       `moderate` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation")),
       `high` = 
         crossing(released_test = "Released after second test"),
       `maximum` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation"))
  ) %>%
  bind_rows(.id = "stringency") %>%
  mutate(stage_released = "Infectious",
         stringency = fct_inorder(stringency)) 

probs        <- c(0.025,0.25,0.5,0.75,0.975)

mv2gamma <- function(mean, var){
  list(shape = mean^2/var,
       rate  = mean/var,
       scale = var/mean) 
}

gamma2mv <- function(shape, rate=NULL, scale=NULL){
  if (is.null(rate)){
    rate <- 1/scale
  }
  
  list(mean = shape/rate,
       var  = shape/rate^2)
}


time_to_event <- function(n, mean, var){
  if (var > 0){
    parms <- mv2gamma(mean, var)
    return(rgamma(n, shape = parms$shape, rate = parms$rate))
  } else{
    return(rep(mean, n))
  }
}

time_to_event_lnorm <- function(n, meanlog, sdlog){
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}

gen_screening_draws <- function(x){
  n <- nrow(x)
  
  # generate screening random draws for comparison
  x <- mutate(x, 
              screen_1 = runif(n, 0, 1),  # on arrival
              screen_2 = runif(n, 0, 1))  # follow-up
}

# given infection histories above, what proportion of travellers end up being 
# caught at each step in the screening process?

calc_outcomes <- function(x){
 # browser()
  # generate required times for screening 
  
  # what's the probability of detection at each test time given a value of CT?
  x_ <- x %>%
   inner_join(trajectories$models, by=c("sec_idx"="idx","type")) %>% 
    select(-data) %>% 
    mutate(test_t  = pmax(test_t,sec_traced_t),
           test_q  = test_t-sec_exposed_t,
           ct      = pmap_dbl(.f=calc_sensitivity,list(model=m,x=test_q))) %>% 
    mutate(detection_range = case_when(sens_LFA=="higher" ~ as.numeric(as.character(cut(ct,
                                                                          breaks=c(-Inf,27,30,35,Inf),
                                                                          labels=c("0.95","0.65","0.3","0")))),
                                       sens_LFA=="lower" ~ as.numeric(as.character(cut(ct,
                                                                                       breaks = c(-Inf,20,25,30,35,Inf),
                                                                                       labels=c("0.824","0.545","0.083","0.053","0")))))) %>% 
    mutate(test_p=case_when(assay=="PCR"&ct<35~1,
                            assay=="PCR"&ct>=35~0,
                            assay=="LFA"~detection_range,
                            TRUE~NA_real_)) %>% 
    mutate(screen      = runif(n(), 0, 1)) %>% 
    mutate(test_label  = detector(pcr = test_p,  u = screen)) 
 
  #browser()
  return(x_)
}

detector <- function(pcr, u = NULL, spec = 1){
  
  if (is.null(u)){
    u <- runif(n = length(pcr))
  }
  
  # true positive if the PCR exceeds a random uniform
  # when uninfected, PCR will be 0
  TP <- pcr > u
  
  # false positive if in the top (1-spec) proportion of random draws
  FP <- (pcr == 0)*(runif(n = length(pcr)) > spec)
  
  TP | FP
}


make_delay_label <- function(x,s){
  paste(na.omit(x), s)
}


capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

make_trajectories <- function(n_cases){
  #simulate CT trajectories
  #browser()
  traj <- data.frame(idx=1:n_cases) %>% 
    crossing(start=0,type=c("symptomatic","asymptomatic") %>% 
               factor(x = .,
                      levels = .,
                      ordered = T)) %>% 
    mutate(u = runif(n(),0,1)) %>%
    # duration from: https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30172-5/fulltext
    # scaling of asymptomatics taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    mutate(end=case_when(type == "symptomatic"  ~ qnormTrunc(p = u, mean=17, 
                                                             sd=approx_sd(15.5,18.6), min = 0),
                         type == "asymptomatic" ~ qnormTrunc(p = u, mean=17*0.6, 
                                                             sd=approx_sd(15.5,18.6), min = 0))) %>% 
    # incubation period from https://bmjopen.bmj.com/content/10/8/e039652.info
    mutate(onset_t=qlnormTrunc(p = u,
                               meanlog=1.63,
                               sdlog=0.5,
                               min = 0,
                               max=end)) %>% 
    pivot_longer(cols = -c(idx,type,u),
                 values_to = "x") %>% 
    # peak CT taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    mutate(y=case_when(name=="start"   ~ 40,
                       name=="end"     ~ 40,
                       name=="onset_t" ~ rnorm(n=n(),mean=22.3,sd=4.2))) 
  
  models <- traj %>%
    nest(data = -c(idx,type,u)) %>%  
    dplyr::mutate(
      # Perform loess calculation on each individual 
      m  = purrr::map(data, ~splinefunH(x = .x$x, y = .x$y,
                                       m = c(0,0,0))),
      rx = purrr::map(data, ~range(.x$x)),
      ry = purrr::map(data, ~range(.x$y)),
      
      inf_period=purrr::pmap(.f = infectious_period,
                             .l = list(model = m, 
                                       rx = rx, # range in x needed for when we're outside the range of the sampled end time
                                       ry = ry))) %>% # range in y is for debugging
    unnest_wider(inf_period)
  
  return(list(traj=traj,models=models))
}


## just making sure the proportion of cases are secondary or not
make_sec_cases <- function(prop_asy, trajectories,n_sec_cases){
  
  props <- c("symptomatic"  = (1 - prop_asy),
             "asymptomatic" = prop_asy)
  
  res <- lapply(names(props), 
                function(x){
                  filter(trajectories, type == x) %>%
                    sample_frac(., size = props[[x]])
                })
  
  res <- do.call("rbind",res) %>% 
    sample_n(n_sec_cases)
  
}


make_released_quantiles <- function(x, vars){
  
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  
  dots <- append(dots1, dots2)
  
  x_count <- x %>%
    dplyr::ungroup(.) %>%
    dplyr::select(!!! dots) %>%
    dplyr::group_by_all(.) %>%
    dplyr::count(.) 
  
  x_count %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-n) %>%
    dplyr::ungroup %>%
    as.list %>%
    map(unique) %>%
    expand.grid %>%
    dplyr::left_join(x_count) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    tidyr::nest(data = c(sim, n)) %>%
    dplyr::mutate(
      Q = purrr::map(
        .x = data,
        .f = ~quantile(.x$n, probs = probs)),
      M = purrr::map_dbl(
        .x = data, 
        .f = ~mean(.x$n))) %>%
    tidyr::unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup(.)
}

make_quantiles <- function(x, y_var,vars, sum = TRUE){
  #browser()
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  y_var <- as.name(y_var)
  dots  <- append(dots1, dots2)
  
  if (sum){
    x <- x %>%
      dplyr::select(!!! dots, y_var) %>%
      group_by_at(.vars = vars(-y_var)) %>% 
      summarise(!!y_var := sum(!!y_var, na.rm=T),n=n()) %>% 
      mutate(!!y_var:=!!y_var/n)
  }
  
  x_days <- x %>%
    dplyr::select(!!! dots, !! y_var) 
  
  x_days %>%
    nest(data = c(!!y_var, sim)) %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .x[[y_var]],
                                                probs = probs)),
           M = map_dbl(.x = data, ~mean(.x[[y_var]]))) %>%
    unnest_wider(Q) %>% 
    # mutate(Q = purrr::map(.x = data, ~bayestestR::hdi( .x[[y_var]],
    #                                             ci = c(0.95,0.5))),
    #        Mean = map_dbl(.x = data, ~mean(.x[[y_var]])),
    #        MAP = map_dbl(.x=data, ~as.numeric(bayestestR::point_estimate(.x[[y_var]],centrality="MAP")))) %>% 
    #unnest(Q) %>% 
    #pivot_wider(names_from = CI,values_from=c(CI_low,CI_high)) %>% 
    dplyr::select(-data) 
  
}


delay_to_gamma <- function(x){
  ans <- dplyr::transmute(x, left = t - 0.5, right = t + 0.5) %>%
    dplyr::mutate(right = ifelse(right == max(right), Inf, right)) %>%
    {fitdistcens(censdata = data.frame(.),
                 distr = "gamma", 
                 start = list(shape = 1, rate = 1))} 
  
  gamma2mv(ans$estimate[["shape"]],
           ans$estimate[["rate"]])
}


rtgamma <- function(n = 1, a = 0, b = Inf, shape, rate = 1, scale = 1/rate){
  #browser()
  p_b <- pgamma(q = b, shape = shape, rate = rate)
  p_a <- pgamma(q = a, shape = shape, rate = rate)
  
  u   <- runif(n = n, min = p_a, max = p_b)
  q   <- qgamma(p = u, shape = shape, rate = rate)
  
  return(q)
}


check_unique_values <- function(df, vars){
  # given a data frame and a vector of variables to be used to facet or group, 
  # which ones have length < 1?
  
  l <- lapply(X = vars, 
              FUN =function(x){
                length(unique(df[, x]))
              })
  
  vars[l > 1]
  
}


waning_piecewise_linear <- function(x, ymax, ymin, k, xmax){
  
  if (ymin == ymax){
    Beta = c(0, ymin)
  } else {
    
    Beta <- solve(a = matrix(data = c(xmax, 1,
                                      k,    1),    ncol = 2, byrow = T),
                  b = matrix(data = c(ymin, ymax), ncol = 1))
  }
  
  (x >= 0)*pmin(ymax, pmax(0, Beta[2] + Beta[1]*x))
  
}

waning_points <- function(x, X, Y, log = FALSE){
  
  if (length(X) != length(Y)){
    stop("X and Y must be same length")
  }
  
  if (length(Y) == 1){
    return(rep(Y, length(x)))
  }
  
  if (log){
    Y <- log(Y)
  }
  
  Beta <- solve(a = cbind(X, 1), b = matrix(Y,ncol=1))
  
  Mu <- Beta[2] + Beta[1]*x
  if (log){
    Mu <- exp(Mu)
  }
  (x >= 0)*pmax(0, Mu)
  
}



summarise_simulation <- function(x, faceting, y_labels = NULL){
  
  if(is.null(y_labels)){
    # if none specified, use all.
    y_labels_names <- grep(x=names(x), pattern="^trans_pot_", value = T)
  } else {
    y_labels_names <- names(y_labels)
  }
  
  all_grouping_vars <- all.vars(faceting)
  
  # if (!any(grepl(pattern = "type", x = all_grouping_vars))){
  #   all_grouping_vars <- c(all_grouping_vars, "type")
  # }
  
  x_summaries <-
    as.list(y_labels_names) %>%
    set_names(., .) %>%
    lapply(X = ., 
           FUN = function(y){
             make_quantiles(x,
                            y_var = y, 
                            vars = all_grouping_vars)})
  
  if (any(grepl(pattern = "type", x = all_grouping_vars))){
    
    x_summaries_all <- as.list(y_labels_names) %>%
      set_names(., .) %>%
      lapply(X = ., 
             FUN = function(y){
               make_quantiles(
                 mutate(x,
                        type = "all"),
                 y_var = y, 
                 vars = all_grouping_vars)})
    
    x_summaries <- map2(.x = x_summaries,
                        .y = x_summaries_all,
                        .f = ~bind_rows(.x, .y))
    
  }
  
  bind_rows(x_summaries, .id = "yvar")
  
}

calc_sensitivity <- function(model, x){
  #browser()
  if(!is.na(x)){
  s <- model(x)
  } else {
    s <- NA_real_
  }
  
  return(s)
}



read_results <- function(results_path){
  #browser()
  list(here::here("results", results_path, "results.rds"),
       here::here("results", results_path, "input.rds")) %>%
    map(read_rds) %>%
    map(bind_rows) %>%
    {inner_join(.[[1]], .[[2]])}
}


test_times <- function(multiple_tests,tests,tracing_t,sec_exposed_t,quar_dur,sampling_freq = 1, max_tests = 14, n_tests){
  #browser()
  
  if(multiple_tests){
    test_timings <- data.frame(test_t = seq(from=tracing_t,to=(tracing_t+max_tests-1),by=sampling_freq)) %>% 
    mutate(test_no = paste0("test_", row_number())) %>% 
    mutate(have_test = row_number()<=n_tests) %>% 
    filter(have_test)
  } else {
    test_timings <- data.frame(test_t=sec_exposed_t+quar_dur) %>% 
      mutate(test_no = paste0("test_", row_number())) %>% 
      mutate(have_test = tests)
  }
  
  return(test_timings)
}

earliest_pos <- function(df){
  #browser()
  x_q <- unlist(df[(df$test_label),"test_q"])
  
  if (length(x_q) == 0L){
    return(Inf)
  } else {
    return(min(x_q))
  }
}

earliest_pos2 <- function(df){
  #browser()
  x_q <- df %>% filter(test_label==TRUE) 
  
  if (nrow(x_q) == 0L){
    return(data.frame(test_no="None",test_p=0,test_q=Inf))
  } else {
    return((x_q %>% select(test_no,test_p,test_q) %>% slice_min(test_q)))
  }
}


infectious_period <- function(model, rx, ry){
  #browser()
  newdata <- data.frame(x = seq(from = 0, to = 30, 0.1))
  
  newdata$y_pred <- model(newdata$x)
  
  
  newdata <- mutate(newdata,
                    y_pred = ifelse(x > rx[2], 40, y_pred))
  
  newdata %>%
    mutate(infectious = y_pred <= 30) %>%
    filter(infectious) %>%
    summarise(inf_start = min(x),
              inf_end   = max(x))
  
}

calc_overlap <- function(x){
#browser()
  x <-
    x %>% 
    rowwise() %>% 
    mutate(
      symp_overlap = case_when(
        type == "symptomatic" ~ Overlap(
          x = c(inf_start,
                inf_end),
          y = c(onset_q,
                symp_end_q)
        ) * adhering_iso,
        type == "asymptomatic" ~ 0
      ),
      symp_overlap = ifelse(is.infinite(inf_start) &
                              is.infinite(inf_end), 0, symp_overlap),
      quar_overlap = case_when(
        #If symptomatic before tracing, no quarantine
        !multiple_tests & onset_q < traced_q & type == "symptomatic" ~ 0, 
        #If symptomatic during quarantine, quarantine ends at onset
        !multiple_tests & onset_q > traced_q & onset_q < quar_end_q & type == "symptomatic" ~ 
        Overlap(
          x = c(inf_start,
                inf_end),
          y = c(traced_q,
                onset_q)
        ) * adhering_quar,
        #If symptomatic after quarantine, quarantine is completed in full
        !multiple_tests & onset_q > traced_q & onset_q > quar_end_q & type == "symptomatic" ~ 
          Overlap(
            x = c(inf_start,
                  inf_end),
            y = c(traced_q,
                  quar_end_q)
          ) * adhering_quar,
        #If asymptomatic, quarantine is completed in full
        !multiple_tests & type == "asymptomatic" ~ Overlap(
          x = c(inf_start,
                inf_end),
          y = c(traced_q,
                quar_end_q)
        ) * adhering_quar,
        multiple_tests ~ 0
      ),
      quar_overlap = ifelse(is.infinite(inf_start) &
                              is.infinite(inf_end), 0, quar_overlap),
      test_overlap = case_when(
        #If testing and symptoms occur before the test, no-post test isolation
        tests & onset_q < earliest_q &type =="symptomatic" ~ 0,
        #If testing and symptoms occur after a positive test, end post-positive test isolation and begin symptomatic
        tests & onset_q > earliest_q & onset_q < test_iso_end_q & type == "symptomatic" ~
          Overlap(
        x = c(inf_start,
              inf_end),
        y = c(earliest_q,
              onset_q)
      ) * adhering_iso,
      #If testing and symptoms occur after the end of post-positive test isolation, complete full self-isolation
      tests & onset_q > earliest_q & onset_q > test_iso_end_q & type == "symptomatic" ~
        Overlap(
          x = c(inf_start,
                inf_end),
          y = c(earliest_q,
                test_iso_end_q)
        ) * adhering_iso,
      #If testing and asymptomatic, complete full self-isolation
      tests & type == "asymptomatic" ~
        Overlap(
          x = c(inf_start,
                inf_end),
          y = c(earliest_q,
                test_iso_end_q)
        ) * adhering_iso,
      !tests ~ 0),
      test_overlap = ifelse(is.infinite(inf_start) &
                              is.infinite(inf_end), 0, test_overlap),
      max_overlap = case_when(
        #If quarantine and testing, time in quarantine is sum of quarantine  symptomatic and test self-isolation durations
        tests &
          !multiple_tests ~ symp_overlap + quar_overlap + test_overlap,
        #If daily testing, time in isolation is sum of symptomatic and test self-isolation durations
        tests &
          multiple_tests  ~ symp_overlap + test_overlap,
        #If quarantine only, sum of symptomatic and quarantine overlap
        !tests            ~ symp_overlap + quar_overlap
      )
    )
  
}

approx_sd <- function(x1, x2){
      (x2-x1) / (qnorm(0.95) - qnorm(0.05) )
   }
