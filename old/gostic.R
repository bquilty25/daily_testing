# gostic figure 2 s1 data for PCR screening detection

gostic <- read_csv("traveller_screening/2020_nCov/Fig2S1_sourceData.csv")

gostic_dat <- gostic %>%
  dplyr::mutate(outcome = stringr::str_extract(outcome, "[^:]*"),
                fever   = forcats::fct_inorder(fever)) %>%
  dplyr::filter(grepl(pattern = "5\\.5", x = meanIncubate)) %>%
  dplyr::mutate(y = yymax - yymin) %>%
  dplyr::group_by(outcome, days.since.exposed, fever) %>%
  dplyr::summarise_at(.vars = dplyr::vars(y),
                      .funs = sum) %>%
  dplyr::mutate(fever = readr::parse_number(as.character(fever))/100) %>%
  dplyr::filter(outcome == "detected") %>%
  dplyr::ungroup(.)

fit_curve <- function(x,y){
  
  b <- matrix(y,ncol = 1)
  A <- cbind(1, matrix(x, ncol=1), matrix(x,ncol=1)^2)
  
  solve(A,b)
  
}

predict_curve <- function(x,y,xnew = x){
  
  X <- fit_curve(x,y)
  
  Anew  <- cbind(1, matrix(xnew, ncol=1), matrix(xnew,ncol=1)^2)
  
  data.frame(x = xnew, y = as.numeric(Anew %*% X))
  
}

gostic_pred <- gostic_dat %>%
  nest(data = -c(days.since.exposed, outcome)) %>%
  mutate(pred = map(.x = data, ~predict_curve(x    = .x$fever, 
                                              y    = .x$y,
                                              xnew = c(0,1)))) %>%
  unnest(pred) %>%
  select(-data) %>%
  rename(fever = x)

gostic_pred %<>%
  dplyr::bind_rows(gostic_dat) %>%
  dplyr::bind_rows(., tidyr::crossing(
    outcome = .$outcome,
    days.since.exposed =
      max(.$days.since.exposed) + 
      seq(0.1, 30, by = 0.1),
    fever = .$fever)) %>%
  dplyr::arrange(outcome, fever, days.since.exposed) %>%
  tidyr::fill(y, .direction = "down") 

gostic_plot <- gostic_pred %>%
  ggplot(data= ., aes(x = days.since.exposed, y=y)) +
  geom_line(aes(group = fever,
                lty = (fever %in% unique(gostic_dat$fever))),
            color = lshtm_greens[2]) +
  theme_bw() +
  xlab("Days since exposure") +
  ylab("Probability of detection by\nsyndromic screening\nand self-reporting") +
  geom_text(data = filter(gostic_pred,
                          days.since.exposed == max(days.since.exposed)),
            aes(label = paste("alpha == ", sprintf("%0.2f", fever), "", sep = "")),
            hjust = 0, parse = T, vjust = 0.5,
            nudge_x = 1,
            size = 3) +
  ylim(c(0,NA)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_linetype_manual(values = c("FALSE" = 1, "TRUE" = 2), guide = FALSE)

list("png", "pdf") %>%
  map(~ggsave(filename = paste0("results/gostic.",.x),
              plot=gostic_plot,
              width = 105, height = 70, units="mm",
              dpi = 320))

gostic_symp_only <- gostic_pred %>% dplyr::filter(fever == 0)
