imported_cases

source("../pcr_entry_screening/packages.R")
source("../pcr_entry_screening/plot_functions.R")

imports <- imported_cases %>%
  ungroup() %>% 
  select(-origin_country_iso_code) %>% 
  rename(iso_code = destination_country_iso_code) %>% 
  #select(-X1) %>% 
  mutate(region = countrycode(iso_code, 
                              origin = "iso3c",
                              destination="un.region.name"),
         subregion = countrycode(iso_code,
                                 origin = "iso3c",
                                 destination = "un.regionsub.name"))

# need to do better region names
# 
# source("find_latest_results.R")
# assign(results_name, read.fst(here::here("results", results_file)))


fill_3 <- list(`Less than 1%` = data.frame(xmin = 1, xmax = 1e6) %>%
                 gather(key, value) %>%
                 mutate(ymin = 0, ymax = 0.01 * value),
               `1%-10%` = data.frame(xmin = 1, xmax = 1e6) %>%
                 gather(key, value) %>%
                 mutate(ymin = 0.01*value, ymax = 0.1 * value),
               `10%-100%` = data.frame(xmin = 1, xmax = 1e6) %>%
                 gather(key, value) %>%
                 mutate(ymin = 0.1*value, ymax = 1 * value),
               `Greater than 100%` = data.frame(xmin = 1, xmax = 1e6) %>%
                 gather(key, value) %>%
                 mutate(ymin = value, ymax = Inf)) %>%
  bind_rows(.id = "band") %>%
  mutate(band = fct_inorder(band))

imports_to_plot <- imports %>% 
  filter(incidence_mid > 0) %>%
  mutate(subregion = fct_recode(.f = subregion,
                                CANZUS = "Northern America",
                                CANZUS = "Australia and New Zealand",
                                `Southern and Central Asia` = "Southern Asia",
                                `Southern and Central Asia` = "Central Asia",
                                `Latin America\nand the Caribbean`="Latin America and the Caribbean")) %>%
  mutate(region = ifelse(subregion == "CANZUS", "Oceania and\nNorthern America", region)) %>%
  mutate(region = factor(region, levels = c("Africa",
                                            "Americas",
                                            "Oceania and\nNorthern America",
                                            "Asia",
                                            "Europe"))) %>% 
  mutate(across(starts_with("imported_cases"),
                .fns = list(synd_screen=function(x)x*(0.16)),
                .names = "{fn}_{col}"))

imports_plot <- ggplot(data = imports_to_plot)+
  geom_ribbon(data = fill_3,
              aes(x = value, ymin = ymin, ymax = ymax, fill = band),
              alpha = 0.5) +
  geom_linerange(aes( x    = incidence_mid,
                      #y    = imported_cases_mid,
                      ymin = imported_cases_low,
                      ymax = imported_cases_high))+
  geom_errorbarh(aes(y     = imported_cases_mid,
                     xmin  = incidence_low,
                     xmax  = incidence_high))+ 
  geom_point(aes(x    = incidence_mid,
                 y    = imported_cases_mid),
             shape=21,
             fill="white")+
  geom_point(aes(x    = incidence_mid,
                 y    = imported_cases_mid),
             shape=21,
             fill="black")+
  geom_label_repel(aes(x     = incidence_mid,
                       y     = imported_cases_mid,
                       label = country),
                   label.size = NA, force = 15,
                   color = "white",
                   fill = rgb(0,0,0,0.4),
                   segment.colour = rgb(0,0,0,0.4),
                   segment.size = 0.5, #segment.alpha = 1,
                   size = 2, 
                   min.segment.length = unit(0.1, "lines"))+
  scale_fill_brewer(palette = "Reds",
                    name = "Imports as a percentage of domestic incidence") +
  coord_equal(ylim = c(NA, 1e4),expand = F) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_x_log10(breaks = 10^seq(0,6,by=2),#trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  #coord_trans(x = "log10", y = "log10") +
  labs(x = "Estimated domestic incidence (infections per day)",
       y = "Estimated number of imported infections arriving in community (infections per day)")+
  plotting_theme+
  theme(panel.spacing = unit(1, "lines")) +
  facet_wrap(~region+subregion,ncol=6)

#imports_plot

save_plot(plot = imports_plot, prefix = "who", base = "imports",
          device = "png", width=300, height=200, dpi=300)

