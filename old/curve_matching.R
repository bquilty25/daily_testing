## matching PCR and LFA curves

PCR_curves <- read_csv("data/posterior_samples_ct_threshold_37.csv")
LFA_curves <- read_csv("data/posterior_samples_ct_threshold_28.csv")

standardise <- function(x){
  (x - mean(x))/sd(x)
}

my_dist <- function(df1, df2, xvar, yvar, k=1){
  #  browser()
  X <- unlist(df1[[xvar]] - df2[[xvar]])
  Y <- unlist(df1[[yvar]] - df2[[yvar]])
  
  D <- sqrt(X^2 + Y^2/k)
  D_min <- D == min(D)
  
  # to break ties, sample at random
  sample_n(select(filter(df2, D_min), iter_pcr = iter), size = 1)
}

curves_list <- list(LFA = LFA_curves,
                    PCR = PCR_curves) %>%
  map(~group_by(.x, iter) %>%
        filter(value == max(value)) %>%
        ungroup) %>%
  #map(~head(.x, 100)) %>%
  map(~mutate(.x, 
              value_s = standardise(value),
              diff_s  = standardise(diff),
              value_r = rank(value),
              diff_r  = rank(diff)))

mapping <- bind_cols(curves_list[[1]],
                     bind_rows(
                       lapply(
                         X = group_split(rowwise(curves_list[[1]])),
                         FUN = function(x){
                           my_dist(df1  = x, 
                                   df2  = curves_list[[2]],
                                   xvar = "diff",
                                   yvar = "value_r",
                                   k = 1e8)
                         })
                     )) %>%
  rename(iter_lfa = iter) 

mapping_indices <- mapping %>%
  select(iter_lfa, iter_pcr)

mapping_plot <- 
  mapping_indices %>%
  sample_n(100) %>%
  tibble::rowid_to_column(.) %>%
  #nest(data = -c(rowid, iter_lfa)) %>%
  inner_join(LFA_curves, by = c("iter_lfa" = "iter")) %>%
  #unnest(data) %>%
  select(-X1) %>%
  rename(value_lfa = value) %>%
  inner_join(PCR_curves, by = c("iter_pcr" = "iter", "diff")) %>%
  select(-X1) %>%
  rename(value_pcr = value) %>%
  gather(key, value, value_lfa, value_pcr) %>%
  ggplot(data = ., aes(x = diff, y = value)) +
  geom_line(aes(color = key)) +
  facet_wrap(~rowid) +
  theme_void() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank())

ggsave(filename = "mapping.pdf", plot = mapping_plot)


# now that we have iter_lfa and iter_pcr we have a map to get the nearest PCR for each LFA curve

