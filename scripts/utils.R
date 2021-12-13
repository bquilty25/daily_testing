# Load required packages scripts
pacman::p_load("fitdistrplus","EnvStats","tidyverse","patchwork","here","rriskDistributions","dtplyr","rms","DescTools","MESS","lubridate","lemon","boot","furrr","tidytable","ggtext","fst")

plan(multisession,workers=8)

# # McAloon et al. incubation period meta-analysis
#https://bmjopen.bmj.com/content/10/8/e039652
inc_parms <- list(mu_inc = 1.63,
                  sigma_inc = 0.5)

#### LIVERPOOL INNOVA ANALYSIS ----
innova_liv <- tribble(~left,~right,~pos,~neg,
                      30,35,1,18,
                      25,30,3,9,
                      20,25,18,5,
                      15,20,17,2) %>% 
  pivot_longer(c(pos,neg))  

# random sample from binned intervals 
innova_liv_sim <- innova_liv %>% 
  group_by(left,right,name) %>% 
  uncount(value) %>% 
  ungroup() %>% 
  mutate(id=row_number()) %>% 
  crossing(sim=c(1:100)) %>% 
  mutate(ct=runif(n=n(),min = left,max=right),
         Innova=as.numeric(as.factor(name))-1)

liv_mod <- glm(Innova~ct,data=innova_liv_sim,family="binomial") 

##### KCL ANALYSIS ----
pickering <- readxl::read_xlsx(here("data","pickering_dat.xlsx")) %>% 
  select(-c(`Viral Growth`,...7,...8)) %>% 
  rename("culture"=...6) %>% 
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)ifelse(x=="ND",NA,x)) %>% 
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)case_when(x%in%c(0.5,1,2)~1,
                                         is.na(x)~NA_real_,
                                         TRUE~0)) %>% 
  mutate(id=row_number()) %>% 
  rename(ct=`Ct N1`)

innova_mod <- glm(Innova~ct,data=pickering,family="binomial") 
culture_mod <- glm(culture~ct,data=pickering,family="binomial") 

# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- rriskDistributions::get.beta.par(
  q = c(0.24,  0.38),
  p = c(0.025, 0.975), 
  show.output = F, plot = F) %>%
  as.list

approx_sd <- function(x1, x2){
  (x2-x1) / (qnorm(0.95) - qnorm(0.05) )
}

#bootstrap confidence interval function
boot_ci <- function(x,nrep=100) {

  trueval <- tibble(param=c("mu","size"),
                    mean=c(x$estimate[[2]],
                           x$estimate[[1]])) 
  
  ci <- bootdist(f=x,niter = nrep)$CI %>% 
    as.data.frame() %>% 
    select(-Median) %>% 
    rownames_to_column("param")
  
  left_join(trueval,ci)
}

#https://science.sciencemag.org/content/sci/early/2021/05/24/science.abi5273.full.pdf
#assuming 1 and 3 days are 95% interval:
peak_to_onset <- rriskDistributions::get.norm.par(p=c(0.025,0.975),q=c(1,3),plot = F)

make_trajectories <- function(n_cases=100, n_sims=100, seed=1000,asymp_parms=asymp_fraction,variant=c("delta")){
  
  set.seed(seed)
  #simulate CT trajectories
  #browser()

  inf <- data.frame(sim=1:n_sims) %>% 
    mutate.(prop_asy = rbeta(n = n(),
                            shape1 = asymp_parms$shape1,
                            shape2 = asymp_parms$shape2)) 
  
  inf %<>%  
    mutate.(x = map.(.x = prop_asy,
                   .f = ~make_proportions(prop_asy     = .x,
                                          n_cases      = n_cases))) %>% 
    unnest.(x) 
  
  traj <- inf %>% 
    crossing.(start=0) %>% 
    crossing.(variant=variant) %>% 
    # duration from: https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30172-5/fulltext
    # scaling of asymptomatics taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    # mutate(end=case_when(type == "symptomatic"  ~ qnormTrunc(p = u, mean=17, 
    #                                                          sd=approx_sd(15.5,18.6), min = 0),
    #                      type == "asymptomatic" ~ qnormTrunc(p = u, mean=17*0.6, 
    #                                                          sd=approx_sd(15.5,18.6), min = 0))) %>% 
    # incubation period from https://bmjopen.bmj.com/content/10/8/e039652.info
    mutate.(peak=rnormTrunc(n=n(),mean=3.2,sd=approx_sd(2.4, 4.2),min=0),
           clear=case_when.(type == "symptomatic"  ~ rnormTrunc(n=n(), mean=10.9, 
                                                               sd=approx_sd(7.8, 14.2), min = 0),
                           type == "asymptomatic" ~ rnormTrunc(n=n(), mean=7.8, 
                                                               sd=approx_sd(6.1, 9.7), min = 0)),
           clear=case_when.(variant=="delta"~1.4*clear,
                           #variant=="delta"&name=="end"~0.71*y,
                           TRUE~clear),
           end=peak+clear,
           onset_t=peak+rnormTrunc(n=n(),mean = 2,sd=1.5,min=start-peak,max=end-peak)
           
           #onset_t=qlnormTrunc(p = u,
           #                            meanlog=1.63,
           #                            sdlog=0.5,
           #                            min = 0,
           #                            max=end)
    ) %>% 
    select.(-clear) %>% 
    pivot_longer.(cols = -c(sim,prop_asy,idx,type,variant,onset_t),
                 values_to = "x") %>%
    # peak CT taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    mutate.(y=case_when(name=="start" ~ 40,
                       name=="end"  ~ 40,
                       name=="peak"&variant!="delta" ~ rnorm(n=n(),mean=22.4,sd=approx_sd(20.7, 24)),
                       name=="peak"&variant=="delta" ~ 0.7*rnorm(n=n(),mean=22.4,sd=approx_sd(20.7, 24))
    ))  #multiply Ct by 0.7 for Delta variant (24Ct /34Ct): https://virological.org/t/viral-infection-and-transmission-in-a-large-well-traced-outbreak-caused-by-the-delta-sars-cov-2-variant/724
  
  #browser()
  
  models <- traj %>%
    nest.(data = -c(idx,type,variant,onset_t)) %>%  
    mutate.(
      # Perform approxfun on each set of points
      m  = map.(data, ~approxfun(x=.x$x,y=.x$y))) 
  
  #cannot pivot wider with "m" column - extract and rejoin
  x_model <- models %>% 
    select.(-data)
  
  models <- models %>% 
    select.(-m) %>% 
    unnest.(data,.drop=F) %>%  
    select.(-c(y)) %>% 
    pivot_wider.(names_from=name,values_from = x) %>% 
    left_join.(x_model)

}

inf_curve_func <- function(m,start=0,end=30,trunc_t){
  #browser()
  x <- data.frame(t=seq(start,end,by=0.25)) %>% 
    mutate(u=runif(n=n(),0,1))
  
  #predict CTs for each individual per day
  x$ct <- m(x$t)
  
  #predict culture probability given CTs
  x$culture <-  stats::predict(culture_mod, type = "response", newdata = x)

  sum_inf <- sum(x$culture)
  
  return(list(infectiousness=x,sum_inf=sum_inf))
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

## sample asymp proportions
make_proportions <- function(prop_asy, n_cases){
  
  props <- c("symptomatic"  = (1 - prop_asy),
             "asymptomatic" = prop_asy)
  
  x <- data.frame(type=rbinom(n=n_cases,size=1,prob = prop_asy)) %>% 
    mutate(type=ifelse(type==1,"asymptomatic","symptomatic"),
           idx=row_number())
  
}

propresponsible=function(R0,k,prop){
  qm1=qnbinom(1-prop,k+1,mu=R0*(k+1)/k)
  remq=1-prop-pnbinom(qm1-1,k+1,mu=R0*(k+1)/k)
  remx=remq/dnbinom(qm1,k+1,mu=R0*(k+1)/k)
  q=qm1+1
  1-pnbinom(q-1,k,mu=R0)-dnbinom(q,k,mu=R0)*remx
}

test_times <- function(type,onset_t,sampling_freq=3){
  #browser()
  
  #  if(type=="asymptomatic"){
  #    initial_t <- sample(size=1,x = c(0:29))
  # }else{
  #   initial_t <- sample(size=1,x = c(0:onset_t))
  # }
  
  initial_t <- 0
  
  if(!is.na(sampling_freq)){
  test_timings <- data.frame(test_t = seq(from=initial_t,to=30,by=sampling_freq)) %>% 
    mutate(test_no = paste0("test_", row_number())) 
  } else {
    test_timings <- data.frame(test_t = Inf) %>% 
    mutate(test_no = paste0("test_", row_number())) 
  }

  
  return(test_timings)
}

earliest_pos_neg <- function(df,n_negatives){
  #browser()
  
  x_q <- df %>% filter.(test_label==TRUE)
  
  if (nrow(x_q) == 0L){
    earliest_pos_t <- Inf
  } else {
    earliest_pos_t <- x_q %>% select.(test_no,test_p,test_t) %>% slice_min.(test_t) %>% pull.(test_t)
  }
  
  # earliest_neg_t <- df %>% 
  # filter(test_t>earliest_pos_t,!test_label) %>% 
  # mutate(diff=test_t-lag(test_t)) %>% 
  # mutate(diff=replace_na.(diff,1)) %>% 
  # filter(diff==1) %>% 
  # slice_min(test_t) %>% 
  # slice_max(test_t,n = n_negatives) %>% 
  # pull(test_t)
  
  earliest_neg_t <- df[test_t > earliest_pos_t][test_label==FALSE][, `:=`(diff = test_t - lag(test_t))][, `:=`(diff = replace_na.(diff, 1))][diff == 1, .SD[order(test_t)][frankv(test_t, ties.method = "min", na.last = "keep") <= n_negatives]][, .SD[order(test_t, decreasing = TRUE)][frankv(-test_t, ties.method = "min",  na.last = "keep") <= 1L]][1,test_t]
  
  
  return(c(start_iso=earliest_pos_t,
              end_iso=earliest_neg_t))
}


inf_and_test <- function(traj,sampling_freq=c(NA,3)){
  #browser()
  
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), traj$sim[1]))
  
  traj %>% as.data.frame() %>% 
    mutate(infectiousness = pmap(inf_curve_func, .l = list(m = m,start=start,end=end)))  %>% 
    unnest_wider(infectiousness) %>% 
    ungroup() %>%
    mutate.(norm_sum = (sum_inf - min(sum_inf)) / (max(sum_inf) - min(sum_inf))) %>% 
    #testing
    crossing(sampling_freq = sampling_freq) %>% 
    mutate.(test_times = pmap(
      .f = test_times,
      list(
        sampling_freq = sampling_freq,
        onset_t = onset_t,
        type = type
      )
    )) %>%
    unnest.(test_times,.drop=F) %>%
    mutate.(
      ct = pmap_dbl(.f = calc_sensitivity, list(model = m, x = test_t)),
      test_p = stats::predict(innova_mod, type = "response", newdata = data.frame(ct = ct)),
      test_label = detector(test_p = test_p,  u = runif(n = n(), 0, 1))
    ) 
} 

detector <- function(test_p, u = NULL){
  
  if (is.null(u)){
    u <- runif(n = length(test_p))
  }
  
  # true positive if the test exceeds a random uniform
  # when uninfected, PCR will be 0
  TP <- test_p > u
  
}