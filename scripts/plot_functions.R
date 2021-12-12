## functions used for ploting

# colours
covid_pal <- c("#e66101", "#5e3c99", "#0571b0")
covid_pal2 <- set_names(covid_pal, c("All", "Asymptomatic", "Symptomatic"))
lshtm_greens <- rev(c("#00BF6F","#0d5257"))
col_pal <- RColorBrewer::brewer.pal(n=4,name = "Dark2")

# assay_colours <- RColorBrewer::brewer.pal(n = length(assay_labeller),
#                                          name = "Dark2") %>%
#   set_names(names(assay_labeller))

#extrafont::loadfonts()
pdf.options(useDingbats=FALSE)

plotting_theme <- 
  theme_minimal() +
  theme(panel.border    = element_rect(fill=NA),
        legend.position = "bottom",
        legend.box      = "vertical",
        axis.ticks = element_line())

infectivity_labels <-
  c("infectivity_post" =
      "Transmission potential of secondary cases \nafter release",
    "infectivity_averted" = 
      "Transmission potential of secondary cases \naverted as a result of quarantine and testing",
    "infectivity_avertable" = 
      "Remaining transmission potential of secondary cases \n post-tracing",
     "infectivity_quar" = 
       "Transmission potential in community\ndue to imperfect quarantine adherence",
    "infectivity_pre" =
      "Transmission potential of secondary cases \nprior to being traced",
    "infectivity_total" = 
      "Transmission potential of secondary cases \nin community compared to no quarantine or testing",
    "infectivity_mass"  = 
      "Truncated transmission potential\ndistribution mass",
    "infectivity_avertable" = 
      "Transmission potential occurring\nafter quarantine starts",
    "infectivity_post_release_onset" =
      "Transmission potential averted during\npost-quarantine self-isolation"
  )



test_labeller <- function(x){
  
  mutate(x,
         stringency = factor(stringency,
                             levels = c("none",
                                        "one",
                                        "two"),
                             labels = c("No test",
                                        "One test",
                                        "Two tests"),
                             ordered = T))
}

type_labels <- c("asymptomatic" =
                   "Asymptomatic",
                 "symptomatic" =
                   "Pre-symptomatic")

index_test_labeller <- function(x, newline = FALSE){
  paste0("Delay between index case's\nonset and having a test:",
         ifelse(newline, "\n", " "),
         x,
         ifelse(x == 1, " day", " days"))
}

# delay_scaling_labeller <- function(x, newline = FALSE){
#   paste0("Contact tracing delay",
#          ifelse(newline, "\n", " "),
#          "scaling factor: ",
#          x)
# }
delay_scaling_labeller <- function(x){
  dplyr::case_when(
    x==0   ~ "Instant T&T (0 days)",  
    x==0.5 ~ "T&T delays halved (1.5 days)",
    x==1   ~ "Observed T&T delays (3 days)",
    TRUE   ~ "Unknown")
}

sens_scaling_labels <- c(
    "lower"= "Liverpool pilot LFA sensitivity",
    "higher" ="Oxford/PHE evaluation LFA sensitivity")


waning_labeller <- function(x){
  #paste("Adherence to quarantine guidance:\n",
        dplyr::case_when( 
          x == "waning_none"             ~ "Constant 100% adherence",
          x == "waning_constant"         ~ "Constant 75% adherence",
          x == "waning_canada_total"     ~ "Exponential waning from 100%",
          x == "waning_canada_community" ~ "Exponential decay (community only)",
          TRUE ~ "Unknown")
  #)
}

percentage <- function(x, ...){
  
  ans <- scales::percent(x, ...)
  ans[x > 1] <- ""
  
  ans
  
}

pretty_percentage <- function(x){
  ans <- pretty(x)
  ans[ans <= 1]
}

plotting_func <- function(x=results_df,x_var=NULL,row_vars=NULL,col_vars=NULL,group_var=stringency,probs = c(0.025,0.25,0.5,0.75,0.975)){
  #browser()
  sim_dots <- sym("sim")
  x_dots <-  enquo(x_var) 
  row_dots <- enquo(row_vars)
  col_dots  <- enquo(col_vars)
  group_dots <- enquo(group_var)
  dots  <- enquos(sim_dots,x_var,row_vars,col_vars,group_var)
  
  
  col_pal <- RColorBrewer::brewer.pal(n=4,name = "Dark2")
  names(col_pal) <- c("Post-flight quarantine with LFA test","Post-flight quarantine with PCR test","Post-flight quarantine only","Daily LFA testing")
  
  x_ <- x %>%
    filter(!is.na(!!x_dots)) %>%
    group_by(!!!dots) %>% 
    filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
    summarise(all=sum(pmax(inf_end,flight_arr_t)-pmax(inf_start,flight_arr_t)),
              prop=sum(max_overlap)/all)
  
  x_plot <- x_ %>%
    group_by(!!!x_dots,!!!row_dots,!!!col_dots,!!!group_dots) %>% 
    nest() %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .$prop,
                                                probs = probs)),
           Mean = map_dbl(.x=data,
                          ~mean(.$prop,na.rm=T)),
           SD   = map_dbl(.x=data,
                          ~sd(.$prop,na.rm = T))) %>%
    unnest_wider(Q) %>% 
    dplyr::select(-data) %>% 
    ggplot(aes(x = factor(!!x_dots), y = `50%`)) + 
    geom_linerange(aes(ymin = `2.5%`,
                       ymax = `97.5%`,
                       colour=!!group_dots),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
    geom_linerange(aes(ymin = `25%`,
                       ymax = `75%`,
                       colour=!!group_dots),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
    geom_point(aes(y = `50%`,colour=!!group_dots),
               #pch="-",
               size=1.5,
               position=position_dodge(width=0.5)) +
    scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
    facet_nested(nest_line=T,
                 rows = vars(!!row_dots), 
                 cols = vars(!!col_dots),
                 labeller = labeller(
                   type = capitalize,
                   sens_scaling=sens_scaling_labeller,
                   delay_scaling = delay_scaling_labeller,
                   adherence_quar =
                     c("1" = "100% adhere\nto quarantine",
                       "0.5" =
                         "50% adhere\nto quarantine",
                       "0" =
                         "0% adhere\nto quarantine"),
                   adherence_iso =
                     c("1" = "100% adhere\nto isolation",
                       "0.67" =
                         "67% adhere\nto isolation",
                       "0" =
                         "0% adhere\nto isolation")
                 )) +
    plotting_theme+
    scale_colour_manual(breaks = names(col_pal),
                        values= col_pal)
  
  x_plot + labs(x=case_when(quo_text(x_dots)=="quar_dur"~"Quarantine required until n days have passed since exposure",
                            TRUE ~ "Daily LFA tests for n days after tracing"),
                y="Transmission potential averted")
}

make_release_figure <- function(x_summaries,
                                input,
                                xlab = "Days in quarantine",
                                ylab = "",
                                text_size = 2.5,
                                text_angle = 45,
                                h_just = 0,
                                log_scale = FALSE,
                                hline = 0,
                                faceting = NULL,
                                percent = FALSE){
  
  x_summaries %<>% test_labeller # should this be in the facet call?
  
  # how to do presymptomatic
  
  
  facet_vars <- all.vars(faceting)
  
  if ("type" %in% facet_vars){
    x_summaries %<>% mutate(type = factor(type,
                                          levels = c("asymptomatic",
                                                     "symptomatic"),
                                          labels = c("Asymptomatic",
                                                     "Presymptomatic")))
  }
  
  
  figure <-  
    ggplot(data=x_summaries, aes(x = time_since_exp, 
                                 y = `50%`, 
                                 color = stringency)) +
    geom_hline(aes(yintercept=1), linetype=hline)+
    geom_linerange(aes(ymin  = `2.5%`, 
                       ymax  = `97.5%`,
                       group = stringency),
                   position  = position_dodge2(width = 1),
                   alpha     = 0.3,
                   size      = 3) +
    geom_linerange(aes(ymin  = `25%`,
                       ymax  = `75%`,
                       group = stringency),
                   position  = position_dodge2(width = 1),
                   alpha     = 0.5,
                   size      = 3) +
    geom_point(pch           = "-",
               size          = 12,
               position      = position_dodge2(width = 1),
               aes(y         = `50%`,
                   group     = stringency)) +
    scale_x_continuous(breaks = breaks_width(2))+
    scale_color_manual(name = "Number of negative tests required for release",
                       values = covid_pal)+
    theme_minimal()+
    theme(axis.ticks = element_line(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(fill=NA),
          legend.position = "bottom",
          strip.placement = "outside") +
    ggplot2::labs(x = xlab,
                  y = ylab) +
    xlab("Days in quarantine\n(including 1 day delay on testing results)")
  
  figure <- figure + 
    facet_nested(nest_line = TRUE,
                 facets = faceting,
                 labeller = labeller(index_test_delay = index_test_labeller,
                                     delay_scaling    = delay_scaling_labeller,
                                     waning           = waning_labeller))
  
  
  return(figure)
  
}

plot_data <- function(input, 
                      x_summaries,
                      main_scenarios = NULL){
  
  dat <- x_summaries  %>%
    inner_join(input) %>% # should really carry these through when summarising
    mutate(time_since_exp_ = 
             ifelse(is.na(time_since_exp),
                    0,
                    time_since_exp),
           time_in_iso = 
             first_test_delay + 
             time_since_exp_+
             screening)
  
  
  if (!is.null(main_scenarios)){
    main_scenarios %<>% dplyr::select(-one_of("released_test")) %>% distinct
    dat <- left_join(dat, main_scenarios)
  }
  
  dat %>%  
    #tidyr::nest(data = -c(first_test_delay, time_since_exp)) %>%
    dplyr::mutate(delays = paste(first_test_delay, "&",
                                 first_test_delay + time_since_exp)) %>%
    #tidyr::unnest(data) %>%
    dplyr::mutate(time_in_iso = factor(time_in_iso, 
                                       levels = sort(unique(.$time_in_iso)),
                                       ordered = T)) %>%
    dplyr::filter(M!=0) %>%  # if the mean is zero, this group is empty
    return
}

make_arrivals_table <- function(x, table_vars = c("country")){
  x %>%
    mutate(sim = factor(sim, levels = 1:n_arrival_sims)) %>%
    group_by_at(.vars = vars(sim, one_of(table_vars))) %>%
    count(.drop = F) %>%
    group_by_at(.vars = vars(-n, -sim)) %>%
    nest() %>%
    mutate(Q = map(data, ~quantile(.x$n, probs = probs))) %>%
    unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    group_by_at(.vars = vars(one_of(table_vars))) %>%
    transmute(value = sprintf("%0.0f (95%%: %0.0f, %0.0f)", `50%`, `2.5%`, `97.5%`)) %>%
    mutate_at(.vars = vars(-c(country, value)), 
              .funs = str_to_title) %>%
    spread(country, value)
}


save_plot <- function(plot   = ggplot2::last_plot(),
                      prefix = stringi::stri_rand_strings(1, length = 8),
                      base   = NULL, # to help identify which figure in the paper
                      device = NULL, # additional devices, e.g. "png" 
                      width  = 210, 
                      height = 210,
                      dpi    = 300,
                      units  = "mm"){
  file <- paste0("results/",
                 paste(prefix, base, sep = "_"),
                 ".pdf")
  
  pdf(file = file,
      width  = measurements::conv_unit(x = width,  units, "inch"),
      height = measurements::conv_unit(x = height, units, "inch"),
      useDingbats = FALSE)
  
  print(plot)
  
  dev.off()
  
  if (length(device) > 0L){
    #img <- pdftools::pdf_render_page(file, dpi = dpi)
    
    img <- magick::image_read_pdf(file, density = dpi)
    
    purrr::map(.x = device,
               .f = 
                 ~magick::image_convert(image = img, format = .x) %>%
                 magick::image_write(., path = sub(pattern = "pdf",
                                                   replacement = .x,
                                                   x = file)))
  }
}



make_days_plots <- 
  function(x, 
           main_scenarios   = NULL,
           plot             = TRUE,
           log_scale        = FALSE,
           text_size        = 2.5,
           xlab             = "Days since exposure\n(including 1 day delay on testing results)",
           sum              = FALSE,
           y_labels         = NULL, # pass in y_vars as a named list
           faceting         = NULL,
           dir              = stringi::stri_rand_strings(1, 8),
           base             = stringi::stri_rand_strings(1, 8)){
    
    if (!dir.exists(paste0("results/",dir))){
      dir.create(paste0("results/",dir))
    }
    
    all_grouping_vars <- all.vars(faceting)
    
    if (sum){
      y_labels <- sub("^Average", "Total", y_labels)
    } 
    
    x_days_summaries <-
      as.list(names(y_labels)) %>%
      lapply(X = ., FUN = function(y){
        map_df(.x = x,
               ~make_released_time_quantiles(.x,
                                             y_var = y, sum = sum,
                                             vars = all_grouping_vars))})
    ## end summaries
    
    fig_data <- x_days_summaries %>% 
      map(~plot_data(input = input, # should we still pass it in? recover?
                     x_summaries = 
                       .x,
                     main_scenarios))
    if (plot){
      
      figs <- map2(
        .x = fig_data,
        .y = y_labels,
        .f = ~make_release_figure(
          x         = .x,
          #input     = input,
          xlab      = xlab,
          text_size = text_size,
          ylab      = .y,
          faceting  = faceting,
          percent   = TRUE) )
      
      
      if (length(figs) > 1){
        legend <- cowplot::get_legend(figs[[1]])
        figs   <- lapply(X = figs, function(x){x + theme(legend.position = "none")})
        
        fig    <- cowplot::plot_grid(plotlist = figs, nrow = 1,
                                     labels = LETTERS[1:length(figs)])
        
        fig    <- cowplot::plot_grid(fig, legend, ncol = 1, rel_heights = c(1, .2))
        
      } else {
        fig <- figs[[1]]
      }
      
      
      list("png", "pdf") %>%
        map(~ggsave(filename = paste0("results/",dir,"/days_plots_",base,".",.x),
                    plot=fig,
                    width  = 60*nrow(distinct(fig_data[[1]][,get.vars(rhs(faceting))]))*
                      length(fig_data), 
                    height = 120*nrow(distinct(fig_data[[1]][,get.vars(lhs(faceting))])), 
                    dpi = 600,
                    units = "mm",
                    device = ifelse(.x == "pdf",
                                    cairo_pdf,
                                    "png")))
    }
    
    return(
      set_names(fig_data, sub(pattern = "infectivity_", 
                              replacement = "", x = names(y_labels)))
      
    )
    
  }


summarise_results <- function(x, reduction = TRUE){
  if (!is.logical(reduction)){
    stop("reduction must be logical")
  }
  x <- mutate_at(x, 
                 .vars = vars(contains("%")),
                 .funs = function(x){
                   reduction*(1 - x) +
                     (1 - reduction)*x})
  if (reduction){
    percentages <- grep(x = names(x), pattern = "%")
    x[,rev(percentages)] <- x[,percentages]  
  }
  
  x
  
}

show_results <- function(x, reduction = TRUE){
  dplyr::select(x, delays,
                one_of(all.vars(faceting)),
                screening, time_in_iso,
                contains("%")) %>%
    group_by(stringency, index_test_delay) %>%
    group_split %>%
    map(summarise_results, reduction = reduction) 
}



ribbon_plot <-
  function(x, 
           y_labels   = NULL, 
           colour_var = "stringency",
           by_type    = FALSE,
           custom_facets = NULL,
           bars  =TRUE
  ){

    
    if (is.null(custom_facets)){
      f_lhs <- c("waning", "index_test_delay", "delay_scaling")
      f_rhs <- c("yvar","stringency")
      
    } else {
      f_lhs <- all.vars(lhs(custom_facets))
      f_rhs <- all.vars(rhs(custom_facets))
    }
    
    if (!by_type){
      x <- filter(x, type == "all")
    } else {
      f_lhs <- c(f_lhs, "type") 
    }
    
    if (!any(f_rhs == "yvar") & length(y_labels) > 0){
      f_rhs <- f_rhs <- c("yvar", f_rhs) 
    }
    
    if (!all(f_lhs == ".")){
      f_lhs <- grep(pattern = ".", x = f_lhs, fixed = T, invert = T, value = T)
    }
    
    
    if (!is.null(y_labels)){
      x <- filter(x, yvar %in% names(y_labels))
      x <- mutate(x, yvar = factor(yvar, levels = names(y_labels), ordered = T))
    }
    
    # here we want to drop anything that's only got one value
    
    drop_unused_terms <- function(y,x){
      
      map_chr(.x = y, .f = function(v){
        if (length(unique(x[[v]])) > 1L | v == "."){
          v
        } else {NA_character_}
      }) %>% na.omit %>% c %>% unique
    }
    
    f_lhs <- drop_unused_terms(f_lhs, x)
    f_rhs <- drop_unused_terms(f_rhs, x)
    
    
    faceting_new <-
      as.formula(
        paste(
          paste(f_lhs,
                collapse = " + "),
          paste(f_rhs,
                collapse = " + "),
          sep = " ~ "
        )
      )
     
    
    x %<>% 
      test_labeller() %>% 
      mutate(type       = capitalize(type))
    
    xlims <- range(x$time_since_exp)
    
    colour_var_sym <- sym(colour_var)
    
    the_plot <- 
      ggplot(data = x, aes(x = time_since_exp,
                           y = M,
                           #color = !!colour_var,
                           fill  = !!colour_var_sym)) +
      facet_nested(nest_line = TRUE,
                   drop      = TRUE,
                   facets    = faceting_new,
                   labeller  = labeller(index_test_delay = index_test_labeller,
                                        delay_scaling    = 
                                          function(x){delay_scaling_labeller(x)},
                                        waning           = waning_labeller,
                                        #stringency       = test_labeller,
                                        type             = capitalize,
                                        yvar             = infectivity_labels,
                                        .multi_line      = TRUE)) +
      theme_minimal() +
      theme(legend.position = "bottom",
            panel.border = element_rect(fill=NA),
            axis.ticks = element_line()) + 
      xlab(expression("Quarantine required until"~italic("t")~"days have passed since exposure")) +
      ylab("Median transmission potential averted") 
    
    if(bars == TRUE){
     the_plot <-  the_plot + 
       geom_linerange(aes(ymin = `2.5%`,
                      ymax = `97.5%`,
                      colour=!!colour_var_sym),size=1,alpha=0.3) +
      geom_linerange(aes(ymin = `25%`,
                      ymax = `75%`,
                     colour=!!colour_var_sym),size=1,alpha=0.5)
    } else {
       the_plot <- the_plot +
         geom_line(aes(y=`2.5%`,colour=!!colour_var_sym),linetype="dashed")+
         geom_line(aes(y=`97.5%`,colour=!!colour_var_sym),linetype="dashed")
     }
    
    the_plot <- the_plot +
      geom_point(aes(y = `50%`,
                    color = !!colour_var_sym),
                 pch="-",
                 size=5) +
      scale_x_continuous(minor_breaks = seq(xlims[1], xlims[2], by = 1),
                         breaks       = seq(xlims[1], xlims[2], by = 7))+
      scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1))
    
    
    if (colour_var == "stringency"){
      the_plot <- the_plot +
        scale_color_manual(name = "Number of tests required before release",
                           values = covid_pal) +
        scale_fill_manual(name = "Number of tests required before release",
                          values = covid_pal)
    } else{ 
      the_plot <- the_plot + 
        scale_color_manual(name = "Type of infection",
                           values = covid_pal2) +
        scale_fill_manual(name = "Type of infection",
                          values = covid_pal2)
    } 
    
    the_plot
  } 

#' log scale
#'
#' Creates a function which returns ticks for a given data range. It uses some
#' code from scales::log_breaks, but in contrast to that function it not only
#' the exponentials of the base b, but log minor ticks (f*b^i, where f and i are 
#' integers), too.
#'
#' @param n Approximate number of ticks to produce
#' @param base Logarithm base
#'
#' @return
#'
#' A function which expects one parameter:
#'
#' * **x**: (numeric vector) The data for which to create a set of ticks.
#'
#' @export
logTicks <- function(n = 5, base = 10){
  # Divisors of the logarithm base. E.g. for base 10: 1, 2, 5, 10.
  divisors <- which((base / seq_len(base)) %% 1 == 0)
  mkTcks <- function(min, max, base, divisor){
    f <- seq(divisor, base, by = divisor)
    return(unique(c(base^min, as.vector(outer(f, base^(min:max), `*`)))))
  }
  
  function(x) {
    rng <- range(x, na.rm = TRUE)
    lrng <- log(rng, base = base)
    min <- floor(lrng[1])
    max <- ceiling(lrng[2])
    
    tck <- function(divisor){
      t <- mkTcks(min, max, base, divisor)
      t[t >= rng[1] & t <= rng[2]]
    }
    # For all possible divisors, produce a set of ticks and count how many ticks
    # result
    tcks <- lapply(divisors, function(d) tck(d))
    l <- vapply(tcks, length, numeric(1))
    
    # Take the set of ticks which is nearest to the desired number of ticks
    i <- which.min(abs(n - l))
    if(l[i] < 2){
      # The data range is too small to show more than 1 logarithm tick, fall
      # back to linear interpolation
      ticks <- pretty(x, n = n, min.n = 2)
    }else{
      ticks <- tcks[[i]]
    }
    return(ticks)
  }
}
