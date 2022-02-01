# Load required packages 

pacman::p_load("fitdistrplus","EnvStats","tidyverse","patchwork","here","rriskDistributions","lubridate","lemon","tidytable","fst","scales","data.table","qs","ggh4x","ggdist","profvis","bayestestR","extraDistr")

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
  
  inf <- rbbinom(n = n_sims,
                 size=1,
                 alpha = asymp_parms$shape1,
                 beta = asymp_parms$shape2) %>% 
    as_tidytable() %>% 
    rename.("asymptomatic"=x) %>% 
    mutate.(sim=row_number.(),
            asymptomatic=as.logical(asymptomatic))
  
  traj <- inf %>% 
    crossing.(start=0) %>% 
    crossing.(variant_info) %>% 
    mutate.(prolif=rnormTrunc(n=n(),mean=prolif_mean,sd=approx_sd(prolif_lb,prolif_ub),min=0),
           clear=rnormTrunc(n=n(), mean=clear_mean, sd=approx_sd(clear_lb,clear_ub), min = 0),
           prolif=ifelse(asymptomatic,prolif*0.8,prolif),
           clear=ifelse(asymptomatic,clear*0.8,clear),
           end=prolif+clear,
           onset_t=prolif+rnormTrunc(n=n(),mean = 2,sd=1.5,min=0,max=end)
    ) %>% 
    select.(-c(prolif_mean, prolif_lb, prolif_ub, clear_mean, clear_lb, clear_ub,clear)) %>% 
    pivot_longer.(cols = -c(sim,asymptomatic,variant,onset_t),
                 values_to = "x") %>%
    mutate.(y=case_when(name=="start" ~ 0,
                        name=="end"   ~ 0,
                        name=="prolif"  ~ rnorm(n=n(),mean=variant_info[,max_vl_mean],sd=approx_sd(variant_info[,max_vl_lb],variant_info[,max_vl_ub])))
            )
            
  models <- traj %>%
    nest.(data = -c(asymptomatic,variant,onset_t)) %>%  
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
    select.(c(sim, variant, asymptomatic, onset_t, prolif, start, end, m)) %>% 
    arrange(sim)

}

inf_curve_func <- function(m,start=0,end=30,trunc_t){
  #browser()
  x <- tidytable(t=seq(start,20,by=1)) %>% 
    mutate.(vl=m(t))
  
  return(infectiousness=x)
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

propresponsible=function(R0,k,prop){
  qm1=qnbinom(1-prop,k+1,mu=R0*(k+1)/k)
  remq=1-prop-pnbinom(qm1-1,k+1,mu=R0*(k+1)/k)
  remx=remq/dnbinom(qm1,k+1,mu=R0*(k+1)/k)
  q=qm1+1
  1-pnbinom(q-1,k,mu=R0)-dnbinom(q,k,mu=R0)*remx
}

earliest_pos_neg <- function(df) {
  #browser()
  if(df$test_to_release[[1]]){
  # find earliest series of negative tests after delay
  end_iso <- df[(test_label == FALSE), .SD[order(t)][
    frankv(t, ties.method = "min",na.last = "keep") <= n_negatives]][
    , .SD[order(t, decreasing = TRUE)][
      frankv(-t, ties.method = "min", na.last = "keep") <= 1L]][,t]

} else {
  # release after delay, regardless of test positivity
  end_iso <- df[, .SD[order(t)][
    frankv(t, ties.method = "min", na.last = "keep") <= 1L]][,t]
  
  }
  
  return(end_iso)
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

#summarise multiple columns
p <- c(0.025, 0.5, 0.975)

p_names <- map_chr(p, ~paste0(.x*100, "%"))

p_funs <- map.(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)


bootf <- function(var) {
  #browser()
  rbind(Hmisc::smean.cl.boot(!!var)) %>% 
    as_tidytable()
}
