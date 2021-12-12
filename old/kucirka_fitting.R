message("PCR Sensitivity curve")

#dat <- read_csv("https://raw.githubusercontent.com/HopkinsIDD/covidRTPCR/master/data/antibody-test-data.csv")

dat <- readr::read_csv('covidRTPCR/data/antibody-test-data.csv',
                       col_type = cols()) %>%
  dplyr::filter(grepl("PCR", test) & !grepl("\\_4", x = study))  %>%
  dplyr::mutate(n_adj=n+nqp,
                test_pos_adj=test_pos+nqp) %>% 
  ## remove estimates without observations
  dplyr::filter(n_adj > 0,
                ## days needs to be above -5
                day > -5,
                ## only use the nasal swabs from Kujawski, not throat swabs
                !(study == "Kujawski" & test == "RT_PCR_oro")) %>% 
  dplyr::mutate(study_idx=paste(study, test, sep="_") %>% as.factor() %>% as.numeric(),
                pct_pos=test_pos_adj/n_adj) %>%
  dplyr::mutate(day = day + 5) 

dat_long <- dat %>% 
  dplyr::select(study, day, n, test_pos) %>%
  tidyr::uncount(n) %>%
  dplyr::group_by(study, day) %>% 
  dplyr::mutate(id  = dplyr::row_number(),
                pos = 0 + (test_pos >= id)) %>%
  dplyr::select(-id, -test_pos) 

dat_long %<>% 
  dplyr::ungroup() %>% 
  dplyr::bind_rows(., 
                   dplyr::distinct(., study) %>%
                     dplyr::mutate(day = 0, pos = 0)) %>%
  dplyr::mutate(study = factor(study))

# fit a new model to kucirka's data
dat_gam <- mgcv::gam(pos ~ s(day, bs = "ps") ,
                     family = "binomial", data = dat_long)

