Regression plots
================
Lucas Nell
23 Oct 2017

-   [Loading model data](#loading-model-data)
-   [Function to calculate confidence intervals](#function-to-calculate-confidence-intervals)
-   [Individual plots for models by species only](#individual-plots-for-models-by-species-only)
    -   [Function to create base plots](#function-to-create-base-plots)
    -   [Creating plot objects](#creating-plot-objects)
-   [Individual plots for models by species and intestinal segment](#individual-plots-for-models-by-species-and-intestinal-segment)
    -   [Objects to create base plots](#objects-to-create-base-plots)
    -   [Creating plot objects](#creating-plot-objects-1)
-   [Individual plots for clearance and absorption](#individual-plots-for-clearance-and-absorption)
-   [Creating and saving final plots](#creating-and-saving-final-plots)
    -   [Function to combine plots](#function-to-combine-plots)
-   [Session info](#session-info)

``` r
# Packages needed for this script
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(phylolm)
    library(ape)
    library(ggplot2)
    library(grid)
    # library(gridExtra)
})
# The `R/get_data.R` file provides functions to retrieve morphometric, clearance, and
# absorption data.
# See `tidy_csvs.md` for more info.
source('R/get_data.R')

# setting `ggplot2` theme
theme_set(theme_classic() %+replace% 
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 10),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))
```

Loading model data
==================

``` r
models <- list(absorp = read_rds('output/models_absorp.rds'),
               pos = read_rds('output/models_pos.rds'),
               spp = read_rds('output/models_spp.rds'))

data <- list(absorp = get_df('absorp') %>% as_tibble,
             pos = lapply(c('prox','med', 'dist'), 
                          function(p) {get_df(.df = 'pos', .pos = p) %>% 
                                  mutate(pos = p)}) %>% 
                 bind_rows %>% as_tibble %>% 
                 select(pos, everything()) %>% 
                 gather('measure', 'value', -pos, -taxon, -diet, -species) %>% 
                 mutate(pos = factor(pos, levels = c('prox','med', 'dist'), 
                                     labels = c('Proximal', 'Medial', 'Distal'))),
             spp = get_df('spp') %>% as_tibble,
             clear = get_df('clear') %>% as_tibble)
```

Function to calculate confidence intervals
==========================================

``` r
# Creates data frame containing 95% CI based on bootstrapping for one model
mod_ci <- function(.model){
    
    stopifnot(is(.model, 'phylolm'))
    
    y_measure <- {paste(.model$formula) %>% discard(~ grepl('~', .x))}[1]
    pos_name <- paste(.model$call) %>% 
        keep(~ grepl('_df', .x)) %>% 
        gsub(pattern = '_df', replacement = '')
    if (pos_name == 'spp') pos_name <- NA
    
    if ('log_mass' %in% colnames(.model$X)) {
        new_data <- rbind(c(1.0, 1.0), c(1.0, 0.0), 
                          rep(mean(.model$X[,'log_mass']), 2))
        new_data_df <- data.frame(taxon = c(1.0, 0.0), 
                                  log_mass = rep(mean(.model$X[,'log_mass']), 2))
    } else {
        new_data <- rbind(c(1.0, 1.0), c(1.0, 0.0))
        new_data_df <- data.frame(taxon = c(1.0, 0.0))
    }
    
    # Column names coinciding with phylogenetic parameters from models lambda, BM, and
    # OUfixed:
    phylo_cols <- which(colnames(.model$bootstrap) %in% c('lambda', 'sigma2', 'alpha'))
    
    ci_matrix <- .model$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(x[-phylo_cols] %*% new_data)
              }) %>% 
        apply(1,
              function(x) as.numeric(quantile(x, probs = c(0.025, 0.975)))) %>% 
        t
    
    out_df <- cbind(predict(.model, new_data_df), ci_matrix) %>% 
        as_data_frame %>% 
        rename_(.dots = setNames(colnames(.), c('estimate', 'low', 'high'))) %>% 
        mutate(measure = y_measure,
               taxon = factor(c(1.0, 0.0), levels = c(0,1), 
                              labels = c('Rodent', 'Bat'))) %>% 
        select(taxon, measure, everything()) %>% 
        mutate(pos = pos_name) %>% 
        select(taxon, pos, measure, everything())

    return(out_df)
}
```

Individual plots for models by species only
===========================================

Function to create base plots
-----------------------------

``` r
taxon_only_no_mass <- function(.model, y_axis_title, title = NULL) {
    y_name <- {paste(.model$formula) %>% discard(~ grepl('~', .x))}[1]
    .p <- mod_ci(.model) %>%
        ggplot(aes(taxon, estimate)) +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, size = 0.5) +
        geom_segment(aes(yend = estimate, 
                         x = as.numeric(taxon) - 0.1, 
                         xend = as.numeric(taxon) + 0.1)) +
        geom_point(data = data_frame(y_name = as.numeric(.model$y), 
                                     taxon = factor(as.integer(.model$X[,'taxonBat']), 
                                                    levels = c(0,1), 
                                                    labels = c('Rodent', 'Bat'))),
                   aes(y = y_name, shape = taxon),
                   position = position_jitter(width = 0.2, height = 0),
                   color = 'black', size = 2, fill = 'gray60') +
        scale_shape_manual(values = c(1, 21)) +
        theme(legend.position = 'none', axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(color = 'black', size = 10, margin = margin(6)),
              axis.title.y = element_text(margin = margin(t = 0, r = -8, 
                                                          b = 0, l = 0))) +
        ylab(y_axis_title)
    if (!is.null(title)) {
        min_x <- min(ggplot_build(.p)$layout$panel_ranges[[1]]$x.range) + 
            0.01 * diff(ggplot_build(.p)$layout$panel_ranges[[1]]$x.range)
        max_y <- max(ggplot_build(.p)$layout$panel_ranges[[1]]$y.range) +
            0.01 * diff(ggplot_build(.p)$layout$panel_ranges[[1]]$y.range)
        .p <- .p + 
            geom_text(data = NULL, label = title, 
                      x = min_x * 1.1, y = max_y * 0.99, hjust = 0, vjust = 1, 
                      size = 14 * (25.4/72), fontface = 'plain')

    }
    return(.p)
}
```

Creating plot objects
---------------------

``` r
# Figure 1A
fig1a <- taxon_only_no_mass(models$spp$int_length_mass, 
                            expression(atop("Intestinal length / body" ~ mass^{0.4},
                                            "(" * cm / g^{0.4} * ")")),
                            'A')
# Figure 1B
fig1b <- taxon_only_no_mass(models$spp$nsa_mass,
                            expression(atop("NSA / body" ~ mass^{0.75},
                                            "(" * cm^2 / g^{0.75} * ")")),
                            'B')
# Figure 4
fig4 <- taxon_only_no_mass(models$spp$vill_area_mass, 
                           expression(atop("Villous surface area / body" ~ mass^{0.75},
                                           "(" * cm^2 / g^{0.75} * ")")))
# Figure 6
fig6 <- taxon_only_no_mass(models$spp$log_total_enterocytes,
                   expression("Total enterocytes (" %*% 10^9 * ")")) +
    theme(axis.title.y = element_text(margin = margin(0, 5.5, 0, 0))) +
    scale_y_continuous(breaks = log(c(5e8, 1e9, 1.5e9)), labels = seq(0.5, 1.5, 0.5),
                       limits = log(c(1, 1.55e9))) +
    coord_trans(y = 'exp')
# Mention that bars represent model predictions at mean log(body mass) among all species
```

Individual plots for models by species and intestinal segment
=============================================================

Objects to create base plots
----------------------------

``` r
# Making data frame of confidence intervals
# (Nesting by parameter, not position, bc the former is how they'll be plotted.)
pos_ci <- lapply(names(models$pos$prox), 
                 function(n) {
                     bind_rows(list(mod_ci(models$pos$prox[[n]]),
                                    mod_ci(models$pos$med[[n]]),
                                    mod_ci(models$pos$dist[[n]])))
                 }) %>% 
    bind_rows %>% 
    mutate(pos = factor(pos, levels = c('prox', 'med', 'dist'), 
                        labels = c('Proximal', 'Medial', 'Distal')))

# Table of y-axis names for each parameter
plot_names <- read_csv('og,new
log_intestinal_diameter,"Intestinal ~ diameter ~ \'(cm)\'"
villus_height,"Villus ~ height ~ \'(mm)\'"
villus_width,"Villus ~ width ~ \'(mm)\'"
crypt_width,"Crypt ~ width ~ \'(mm)\'"
sef,"Surface ~ enlargement ~ factor ~ \'(SEF)\'"
enterocyte_diameter,"Enterocyte ~ diameter ~ \'(\u03BCm)\'"
log_enterocyte_density,"Enterocyte ~ density ~ \'(\' %*% 10^6 * \')\'"
')

spp_pos_plot <- function(.measure) {
    .p <- pos_ci %>%
        filter(measure == .measure) %>% 
        mutate(pos = as.numeric(pos)) %>% 
        group_by(taxon) %>% 
        mutate(pos = pos + ifelse(taxon == 'Bat', 0.2, -0.2)) %>% 
        ungroup %>% 
        ggplot(aes(pos, group = taxon)) + 
        geom_point(data = data$pos %>% 
                       filter(measure == .measure) %>% 
                       mutate(pos = as.numeric(pos)) %>% 
                       group_by(taxon) %>% 
                       mutate(pos = pos + ifelse(taxon == 'Bat', 0.2, -0.2)) %>% 
                       ungroup, 
                   aes(y = value, shape = taxon),
                   color = 'black', size = 2, fill = 'gray60',
                   position = position_jitter(0.075, 0)) +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, size = 0.5) +
        geom_segment(aes(y = estimate, yend = estimate,
                         x = pos - 0.1, xend = pos + 0.1), 
                     size = 0.5) +
        theme(legend.position = 'none', legend.margin = margin(0,0,0,0),
              axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), legend.title = element_blank()) +
        scale_shape_manual(values = c(1, 21)) +
        ylab(eval(parse(text = plot_names[plot_names$og == .measure,]$new))) +
        scale_x_continuous(breaks = 1:3, 
                           labels = c('Proximal', 'Medial', 'Distal'))
    
    return(.p)
}


add_title <- function(.p, .title) {
    if (.p$labels$x == 'pos') {
        suppressMessages({
            .p <- .p + 
                scale_x_continuous(breaks = 1:3, limits = c(0.4, 3.43),
                                   labels = c('Proximal', 'Medial', 'Distal'))
        })
    }
    x_range <- ggplot_build(.p)$layout$panel_ranges[[1]]$x.range
    y_range <- ggplot_build(.p)$layout$panel_ranges[[1]]$y.range
    if (!is.null(.p$coordinates$trans$x)) {
        x_range <- .p$coordinates$trans$x$inverse(x_range)
    }
    if (!is.null(.p$coordinates$trans$y)) {
        y_range <- .p$coordinates$trans$y$inverse(y_range)
    }
    min_x <- min(x_range) + 0.02 * diff(x_range)
    max_y <- max(y_range) - 0.02 * diff(y_range)
    .p <- .p + 
        geom_text(data = NULL, label = .title, 
                  x = min_x, y = max_y, hjust = 0, vjust = 1, 
                  size = 14 * (25.4/72), fontface = 'plain')
    return(.p)
}
```

Creating plot objects
---------------------

``` r
# Plots for each parameter.
# I'm avoiding facets bc they make `ggplotGrob`s annoying to combine
pos_plots <- lapply(unique(data$pos$measure), spp_pos_plot)
names(pos_plots) <- unique(data$pos$measure)

# Figure 1c
fig1c <- {pos_plots$log_intestinal_diameter +
        theme(legend.position = 'bottom', 
              axis.text.x = element_text(color = 'black', size = 10, 
                                     margin = margin(6))) +
        scale_y_continuous(breaks = log(seq(0.6, 1.2, 0.2)), 
                           labels = seq(0.6, 1.2, 0.2)) +
        coord_trans(y = 'exp')} %>% 
    add_title('C')

# Figure 2a
fig2a <- {pos_plots$villus_height +
    theme(legend.position = 'top')} %>% 
    add_title('A')

# Figure 2b
fig2b <- pos_plots$villus_width %>% 
    add_title('B')

# Figure 2c
fig2c <- {pos_plots$crypt_width +
        theme(axis.text.x = element_text(color = 'black', size = 10, 
                                     margin = margin(6)))} %>% 
    add_title("C")

# figure 3
fig3 <- pos_plots$sef +
    theme(legend.position = 'top', 
          axis.text.x = element_text(color = 'black', size = 10, 
                                     margin = margin(6)))

# Figure 5a
fig5a <- {pos_plots$enterocyte_diameter +
    scale_y_continuous(breaks = seq(2e-3, 10e-3, 2e-3), labels = seq(2, 10, 2))} %>% 
    add_title('A')

# Figure 5b
fig5b <- {pos_plots$log_enterocyte_density +
    theme(axis.text.x = element_text(color = 'black', size = 10, margin = margin(6)),
          legend.position = 'bottom')} %>% 
    add_title('B')
```

Individual plots for clearance and absorption
=============================================

``` r
fig7a <- data$clear %>%
    mutate(sef = exp(log_sef), clear = exp(log_clear),
           spp_type = interaction(diet, taxon) %>% 
               recode_factor(`Carb.Bat` =       "Bat, carb", 
                             `Protein.Bat` =    "Bat, protein",
                             `Carb.Rodent` =    "Rodent, carb", 
                             `Protein.Rodent` = "Rodent, protein")) %>% 
    ggplot(aes(sef, clear, shape = spp_type)) +
    geom_point(color = 'black', size = 2, fill = 'gray60') +
    # geom_point(color = 'black', size = 2) +
    guides(shape = guide_legend('Taxon, diet:', nrow = 4)) +
    theme(legend.position = c(0.01, 0.875),
          legend.title = element_text(size = 10, face = 'italic'),
          legend.title.align = 0,
          legend.justification = c(0, 1),
          legend.margin = margin(0,0,0,0),
          legend.key.size = unit(0.75, 'lines')) +
    scale_shape_manual(values = c(21, 25, 1, 6)) +
    scale_x_continuous("Surface enlagement factor (SEF)", 
                       trans = 'log', breaks = c(8, 12, 18)) +
    scale_y_continuous(expression("L-arabinose clearance (\u03BC" * l ~ min^{-1} ~ 
                                      cm^{-2} * ")"),
                       trans = 'log', breaks = c(1, 3, 9))
fig7a <- fig7a %>% 
    add_title('A')

fig7b <- taxon_only_no_mass(models$absorp, 
                   expression(atop(
                       "Fractional absorption /",
                       "total intestinal surface (cm"^2 ~ g^{0.75} * ")")),
                   "B")
```

Creating and saving final plots
===============================

Function to combine plots
-------------------------

``` r
# Combining multiple figures using the naming scheme `figX` where `X` is the 
# figure number I'm interested in plotting
# It also works for single figures.
one_fig <- function(fig_num, which_size = 'first') {
    grob_list <- c(lapply(ls(envir = .GlobalEnv)[grepl(paste0('^fig', fig_num), 
                                                       ls(envir = .GlobalEnv))], 
                        function(n) ggplotGrob(eval(parse(text = n)))),
                   size = which_size)
    # grid.newpage()
    grid.draw(do.call(rbind, grob_list))
}
# Employs the above function, plus saves the output
save_fig <- function(fig_num, which_size = 'first', .seed = NULL, ...) {
    file_name <- sprintf('figs/fig%02d.pdf', fig_num)
    if (!is.null(.seed)) set.seed(.seed)
    # pdf(file_name, ...)
    quartz(type = 'pdf', file = file_name, family = 'Helvetica', ...)
    one_fig(fig_num, which_size)
    invisible(dev.off())
}
```

``` r
save_fig(1, width = 3.875, height = 3.125 * 3, .seed = 1)
save_fig(2, 'last', width = 3.875, height = 3.125 * 3, .seed = 2)
save_fig(3, width = 3.875, height = 3.125, .seed = 3)
save_fig(4, width = 3.875, height = 3.125, .seed = 4)
save_fig(5, width = 3.875, height = 3.125 * 2, .seed = 5)
save_fig(6, width = 3.875, height = 3.125, .seed = 6)
save_fig(7, 'last', width = 3.875, height = 3.125 * 2, .seed = 7)
```

Session info
============

This outlines the package versions I used for this script.

    ## Session info -------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.4.2 (2017-09-28)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2017-10-23

    ## Packages -----------------------------------------------------------------

    ##  package    * version date       source        
    ##  ape        * 4.1     2017-02-14 CRAN (R 3.4.0)
    ##  assertthat   0.2.0   2017-04-11 CRAN (R 3.4.0)
    ##  backports    1.1.1   2017-09-25 CRAN (R 3.4.2)
    ##  base       * 3.4.2   2017-10-04 local         
    ##  bindr        0.1     2016-11-13 CRAN (R 3.4.0)
    ##  bindrcpp   * 0.2     2017-06-17 CRAN (R 3.4.0)
    ##  colorspace   1.3-2   2016-12-14 CRAN (R 3.4.0)
    ##  compiler     3.4.2   2017-10-04 local         
    ##  datasets   * 3.4.2   2017-10-04 local         
    ##  devtools     1.13.3  2017-08-02 CRAN (R 3.4.1)
    ##  digest       0.6.12  2017-01-27 CRAN (R 3.4.0)
    ##  dplyr      * 0.7.4   2017-09-28 CRAN (R 3.4.2)
    ##  evaluate     0.10.1  2017-06-24 CRAN (R 3.4.1)
    ##  ggplot2    * 2.2.1   2016-12-30 CRAN (R 3.4.0)
    ##  glue         1.1.1   2017-06-21 CRAN (R 3.4.0)
    ##  graphics   * 3.4.2   2017-10-04 local         
    ##  grDevices  * 3.4.2   2017-10-04 local         
    ##  grid       * 3.4.2   2017-10-04 local         
    ##  gtable       0.2.0   2016-02-26 CRAN (R 3.4.0)
    ##  hms          0.3     2016-11-22 CRAN (R 3.4.0)
    ##  htmltools    0.3.6   2017-04-28 cran (@0.3.6) 
    ##  knitr        1.17    2017-08-10 CRAN (R 3.4.1)
    ##  labeling     0.3     2014-08-23 CRAN (R 3.4.0)
    ##  lattice      0.20-35 2017-03-25 CRAN (R 3.4.2)
    ##  lazyeval     0.2.0   2016-06-12 CRAN (R 3.4.0)
    ##  magrittr     1.5     2014-11-22 CRAN (R 3.4.0)
    ##  memoise      1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods    * 3.4.2   2017-10-04 local         
    ##  munsell      0.4.3   2016-02-13 CRAN (R 3.4.0)
    ##  nlme         3.1-131 2017-02-06 CRAN (R 3.4.2)
    ##  parallel     3.4.2   2017-10-04 local         
    ##  phylolm    * 2.5     2016-10-17 CRAN (R 3.4.0)
    ##  pkgconfig    2.0.1   2017-03-21 CRAN (R 3.4.0)
    ##  plyr         1.8.4   2016-06-08 CRAN (R 3.4.0)
    ##  purrr      * 0.2.4   2017-10-18 CRAN (R 3.4.2)
    ##  R6           2.2.2   2017-06-17 CRAN (R 3.4.0)
    ##  Rcpp         0.12.13 2017-09-28 CRAN (R 3.4.2)
    ##  readr      * 1.1.1   2017-05-16 CRAN (R 3.4.0)
    ##  rlang        0.1.2   2017-08-09 CRAN (R 3.4.1)
    ##  rmarkdown    1.6     2017-06-15 CRAN (R 3.4.0)
    ##  rprojroot    1.2     2017-01-16 cran (@1.2)   
    ##  scales       0.5.0   2017-08-24 CRAN (R 3.4.1)
    ##  stats      * 3.4.2   2017-10-04 local         
    ##  stringi      1.1.5   2017-04-07 CRAN (R 3.4.0)
    ##  stringr      1.2.0   2017-02-18 CRAN (R 3.4.0)
    ##  tibble       1.3.4   2017-08-22 CRAN (R 3.4.1)
    ##  tidyr      * 0.7.2   2017-10-16 CRAN (R 3.4.2)
    ##  tidyselect   0.2.2   2017-10-10 CRAN (R 3.4.2)
    ##  tools        3.4.2   2017-10-04 local         
    ##  utils      * 3.4.2   2017-10-04 local         
    ##  withr        2.0.0   2017-07-28 CRAN (R 3.4.1)
    ##  yaml         2.1.14  2016-11-12 cran (@2.1.14)
