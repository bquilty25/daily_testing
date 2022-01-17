# Load required packages 

pacman::p_load("fitdistrplus","EnvStats","tidyverse","patchwork","here","rriskDistributions","lubridate","lemon","tidytable","fst","scales","data.table","qs","ggh4x","ggdist","profvis","bayestestR")

seed <- 1000

##### KCL ANALYSIS ----
pickering <- readxl::read_xlsx(here::here("data","pickering_dat.xlsx")) %>% 
  select(-c(`Viral Growth`,...7,...8)) %>% 
  rename("culture"=...6) %>% 
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)ifelse(x=="ND",NA,x)) %>% 
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)case_when(x%in%c(0.5,1,2)~1,
                                         is.na(x)~NA_real_,
                                         TRUE~0)) %>% 
  mutate(id=row_number()) %>% 
  rename(ct=`Ct N1`) %>% 
  mutate(vl=(-(ct-44.34)/3.134))

innova_mod <- glm(Innova~vl,data=pickering,family="binomial") 
culture_mod <- glm(culture~vl,data=pickering,family="binomial") 

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
# boot_ci <- function(x,nrep=100) {
# 
#   trueval <- tibble(param=c("mu","size"),
#                     mean=c(x$estimate[[2]],
#                            x$estimate[[1]])) 
#   
#   ci <- bootdist(f=x,niter = nrep)$CI %>% 
#     as.data.frame() %>% 
#     select(-Median) %>% 
#     rownames_to_column("param")
#   
#   left_join(trueval,ci)
# }

make_trajectories <- function(
  n_cases = 100, n_sims = 100,
  seed = seed,
  asymp_parms = asymp_fraction,
  variant_info, browsing = FALSE
){
 
  if (browsing) browser()
  set.seed(seed)
  #simulate CT trajectories
  
  #' TODO:
  #' sample(
  #'   factor(c("asymptomatic", "symptomatic")),
  #'   n_sims,
  #'   rbeta(n_sims, shape1 = ..., shape2 = ...)
  
  
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
    crossing.(variant_info) %>% 
    mutate.(prolif=rnormTrunc(n=n(),mean=prolif_mean,sd=approx_sd(prolif_lb,prolif_ub),min=0),
           clear=rnormTrunc(n=n(), mean=clear_mean, sd=approx_sd(clear_lb,clear_ub), min = 0),
           prolif=ifelse(type=="asymptomatic",prolif*0.8,prolif),
           clear=ifelse(type=="asymptomatic",clear*0.8,clear),
           end=prolif+clear,
           onset_t=prolif+rnormTrunc(n=n(),mean = 2,sd=1.5,min=0,max=end)
    ) %>% 
    select.(-c(prolif_mean, prolif_lb, prolif_ub, clear_mean, clear_lb, clear_ub,min_ct_mean,min_ct_lb,min_ct_ub,clear)) %>% 
    pivot_longer.(cols = -c(sim,prop_asy,idx,type,variant,onset_t),
                 values_to = "x") %>%
    mutate.(y=case_when(name=="start" ~ 0,
                        name=="end"   ~ 0,
                        name=="prolif"  ~ rnorm(n=n(),mean=variant_info[,max_vl_mean],sd=approx_sd(variant_info[,max_vl_lb],variant_info[,max_vl_ub])))
            )
            
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
    left_join.(x_model) %>%
    select.(c(idx, sim, variant, type, onset_t, prolif, start, end, m))

}

inf_curve_func <- function(m,start=0,end=30,trunc_t){
  #browser()
  x <- tidytable(t=seq(start,end,by=1)) %>% 
    mutate.(u=runif(n=n(),0,1),
            vl=m(t),
            culture=stats::predict(culture_mod, type = "response", newdata = tidytable(vl=vl)),
            infectious_label = rbernoulli(n=n(),p=culture)) 

  sum_inf <- sum(x$infectious_label)
  
  
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
  
  #' TODO: change to sample
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

#' TODO
#'  - want this outside of loop
#'  - eliminate string computation
test_times <- function(type,onset_t,end,sampling_freq=3){
  #browser()
  
  #  if(type=="asymptomatic"){
  #    initial_t <- sample(size=1,x = c(0:29))
  # }else{
  #   initial_t <- sample(size=1,x = c(0:onset_t))
  # }
  
  initial_t <- 0
  
  if(!is.na(sampling_freq)){
  test_timings <- data.frame(test_t = seq(from=initial_t,to=14,by=sampling_freq)) %>% 
    mutate(test_no = paste0("test_", row_number())) 
  } else {
    test_timings <- data.frame(test_t = Inf) %>% 
    mutate(test_no = paste0("test_", row_number())) 
  }

  
  return(test_timings)
}

earliest_pos_neg <- function(df) {
  #browser()
  
  # #df <- cbind.data.frame(test_t,test_label,n_negatives,test_to_release,delay) %>% as_tidytable()
  # 
  # x_q <- df %>% filter.(test_label==TRUE)
  # 
  # if (nrow(x_q) == 0L){
  #   start_iso <- Inf
  # } else {
  #   start_iso <- x_q  %>% slice_min.(test_t) %>% pull.(test_t)
  # }
  
  if(df$test_to_release[[1]]){
  # find earliest series of negative tests after delay
  end_iso <- df %>% 
#    filter.(test_t>=start_iso+delay,test_label==FALSE) %>%
#    mutate.(diff=test_t-lags.(test_t,n=1)) %>% 
#    replace_na.(list(diff = 1)) %>% 
#    filter.(diff==1) %>% 
    slice_min.(n=df$n_negatives[[1]], order_by=test_t) %>% 
    slice_max.(n=1, test_t) %>% 
    pull.(test_t)
  

} else {
  # release after delay, regardless of test positivity
  end_iso <- df %>% 
    filter.(test_t>=start_iso+delay) %>%
    slice_min.(n=1,order_by=test_t) %>% 
    pull.(test_t)
    
  }
  
  # if(is.infinite(start_iso)){
  #   end_iso <- start_iso
  # } 
  
  return(end_iso)
}


inf_and_test <- function(traj,sampling_freq=c(NA,3),end=end){
#browser()
  
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), traj$sim[1]))
  
  traj %>% as.data.frame() %>% 
    mutate(infectiousness = pmap(inf_curve_func, .l = list(m = m,start=start,end=end)))  %>% 
    unnest_wider(infectiousness) %>% 
    ungroup() %>%
    #testing
    crossing(sampling_freq = sampling_freq) %>% 
    mutate.(test_times = pmap(
      .f = test_times,
      list(
        sampling_freq = sampling_freq,
        onset_t = onset_t,
        type = type,
        end = end
      )
    )) %>%
    unnest.(test_times,.drop=F) %>%
    mutate.(
      vl = pmap_dbl(.f = calc_sensitivity, list(model = m, x = test_t)),
      test_p = stats::predict(innova_mod, type = "response", newdata = data.frame(vl=vl)),
      test_label = detector(test_p = test_p,  u = runif(n = n(), 0, 1))
    ) 
} 

inf_in_iso <- function(df,iso_interval){
  #browser()
  df %>%
    mutate.(iso_interval=iso_interval) %>% 
    separate.(iso_interval,into = c("start_iso","end_iso"),sep = ",",convert=TRUE) %>% 
    mutate.(iso=between.(t,start_iso,end_iso-1)) %>% #1 minus upper bound
    summarise.(inf_iso=sum(infectious_label),.by=c(start_iso,end_iso,iso)) %>%
    pivot_wider.(names_from = iso,values_from = inf_iso,names_prefix = "iso") 
}

detector <- function(test_p, u = NULL){
  
  if (is.null(u)){
    u <- runif(n = length(test_p))
  }
  
  # true positive if the test exceeds a random uniform
  # when uninfected, PCR will be 0
  TP <- test_p > u
  
}

#summarise multiple columns
p <- c(0.025, 0.5, 0.975)

p_names <- map_chr(p, ~paste0(.x*100, "%"))

p_funs <- map.(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)
