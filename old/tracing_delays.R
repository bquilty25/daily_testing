### Test and Trace delays

delay_names <- list("P_r", "P_c", "P_t")

tti_file <- "data/NHS_T_T_data_tables_w22.ods"


read_delays <- function(sheet, path, string){
 # browser()
  x <- readODS::read_ods(path = path,
                         sheet = sheet, 
                         skip = 2)
  names(x)[1] <- "desc"
  
  x %<>%
    .[,names(.) != ""] %>% as_tibble %>%
    filter(grepl(x = desc, 
                 pattern = string)) %>%
    select(desc, matches("[01-31]/[01-12]")) %>%
    select(1, n = ncol(.)) %>% 
    mutate(t=c(0.5,1.5,2.5,3.5)) %>% 
    mutate(n=as.numeric(n)) %>% 
    uncount(n)
  
  x
}

if (!all(map_lgl(delay_names, ~file.exists(paste0(.x, ".json"))))){
  #Delay from taking a test to receiving the result
  # result_delay_regional <- 
  #   readODS::read_ods(tti_file,
  #                      sheet = "Table_3", 
  #                      skip = 2)
  # names(result_delay_regional)[1] <- "desc"
  
  index_result_delay <- 
    list(regional  = "Table_4",
         mobile    = "Table_6",
         satellite = "Table_7",
         home      = "Table_8") %>%
    map_df(~read_delays(
      sheet = .x, 
      path = tti_file,
      string = "^Number of test results received .* of taking a test$"),
      .id = "source")
  
  
  
  #Delay from positive test to getting info on close contacts from index
  getting_contact_info <-
    read_delays("Table_13", path = tti_file,
                string = "^Number of people reached")
  
  
  #Delay from getting info to tracing contacts
  tracing_delay <-
    read_delays("Table_16", path = tti_file,
                string = "^Number of people reached")
  
  
  message("Fitting delay distributions")
  
  # gamma distributions of delays
  P_r <- delay_to_gamma(index_result_delay)
  P_c <- delay_to_gamma(getting_contact_info)
  P_t <- delay_to_gamma(tracing_delay)
  
  map(delay_names, ~jsonlite::write_json(x    = get(.x), 
                                         path = paste0(.x, ".json")))
  
} else {
  
  message("Loading delay distributions")
  
  map(delay_names, ~assign(x = .x[[1]],
                           value = jsonlite::read_json(paste0(.x, ".json"),
                                                       simplifyVector = T),
                           env = .GlobalEnv))
  
  
}




