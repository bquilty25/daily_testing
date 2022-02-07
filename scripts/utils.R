# Load required packages 

pacman::p_load("fitdistrplus","EnvStats","tidyverse","patchwork","here","rriskDistributions","lubridate","lemon","tidytable","fst","scales","data.table","qs","ggh4x","ggdist","profvis","bayestestR","extraDistr","emdbook")

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

innova_mod <- glm(Innova~vl,
                  data=pickering,
                  family="binomial") 

innova_higher_mod <- glm(Innova~vl,
                        data=pickering %>% 
                        mutate(vl=vl+2.5),family="binomial") 

test_model_choice <- function(boolean){
  if(boolean){
    innova_higher_mod
  }else{
    innova_mod
  }
}

culture_mod <- glm(culture~vl,data=pickering,family="binomial") 

culture_lower_mod <- glm(culture~vl,
                        data=pickering %>% 
                          mutate(vl=vl-2.5),family="binomial") 

inf_model_choice <- function(boolean){
  #browser()
  if(boolean){
    culture_lower_mod
  }else{
    culture_mod
  }
}


# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- rriskDistributions::get.beta.par(
  q = c(0.24,  0.38),
  p = c(0.025, 0.975), 
  show.output = F, plot = F) %>%
  as.list

approx_sd <- function(x1, x2){
  (x2-x1) / (qnorm(0.95) - qnorm(0.05) )
}


make_trajectories <- function(
  n_cases = 100, n_sims = 100,
  seed = seed,
  asymp_parms = asymp_fraction,
  variant_info, browsing = FALSE
){
 
  if (browsing) browser()
  
  set.seed(seed)
  #simulate CT trajectories
  
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
  x <- tidytable(t=seq(start,end,by=1)) %>% 
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
  if(df$n_negatives[[1]]>0){
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

boot_func <- function(x,b = 1000,stat="mean",thresh=0){
  #browser()
 if(stat=="mean"){
 x_samp <- matrix(sample(x,size=b*length(x),replace=T),ncol=b,nrow=length(x))
 
 } else if(stat=="prop_greater"){
   
  x_samp <- matrix(sample(x>thresh,size=b*length(x),replace=T),ncol=b,nrow=length(x))
 }else{
   
 }
  
 mean <- colMeans(x_samp)
 Mean <- mean(mean)
 Lower <- quantile(mean,0.025)
 Upper <- quantile(mean,0.975)
 
 return(list(Mean=Mean,
             Lower=Lower,
             Upper=Upper))
}

summarise_func <- function(x){
  #browser()
  Mean <- mean(x)
  
  Median <- median(x)
  
  Lower <- quantile(x,0.025)
  
  Upper <- quantile(x,0.975)
  
  
  return(list(Mean=Mean,
              Median=Median,
              Lower=Lower,
              Upper=Upper))
}

n_negatives_lab <- function(x) {
  
    paste0(x, " days negative")
}  

delay_lab <- function(x) { paste0(x, " days wait")}
