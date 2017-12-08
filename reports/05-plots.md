Regression plots
================
Lucas Nell
08 Dec 2017

-   [Loading model data](#loading-model-data)
-   [Standardize plot titles](#standardize-plot-titles)
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
    library(gridExtra)
})
# Functions `get_df`, `get_tr`, `filter_tr`, `cp_mat`
source('R/get_data.R')
# Functions `pval`, `ci`, `summ_df`, `jack_phylolm`, and `predict_ci`
source('R/model_summaries.R')
# Custom version of ape::corphylo
suppressMessages(devtools::load_all('corphyloCpp'))

# setting default `ggplot2` theme
theme_set(theme_classic() %+replace% 
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 11, face = 'italic'),
                    legend.title = element_blank(),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))
```

> Below, the `predict_ci` function creates a data frame of predictions and 95% CI from a `phylolm` object. See its documentation in [`R/model_summaries.R`](R/model_summaries.R) for more information.

Loading model data
==================

`phylolm` objects saved from [`04-phylo_regressions`](04-phylo_regressions.md):

``` r
models <- list(absorp = read_rds('output/models_absorp.rds'),
               pos = read_rds('output/models_pos.rds'),
               spp = read_rds('output/models_spp.rds'))
# I'm not plotting this one here
models$pos$prox$crypt_width_pagel <- NULL
```

Data frames used for each model fit:

> See [`../R/get_data.R`](../R/get_data.R) for what `get_df` does.

``` r
data <- list(absorp = get_df('absorp') %>% as_tibble,
             pos = lapply(c('prox','med', 'dist'), 
                          function(p) {get_df(.df = 'pos', .pos = p) %>% 
                                  mutate(pos = p)}) %>% 
                 bind_rows %>% 
                 as_tibble %>% 
                 select(pos, everything()) %>% 
                 gather('measure', 'value', -pos, -clade, -diet, -species, -log_mass) %>% 
                 mutate(pos = factor(pos, levels = c('prox','med', 'dist'), 
                                     labels = c('Proximal', 'Medial', 'Distal'))),
             spp = get_df('spp') %>% as_tibble,
             clear = get_df('clear') %>% as_tibble)
```

Standardize plot titles
=======================

The below function adds titles to the top left corner of plots.

-   `.p`: `ggplot` object to add the title to.
-   `.title`: String to add as the plot title.
-   `.mult`: List specifying multipliers for offset for x- and y-axes. Offsets are the amount of space from the upper-left corner of the plot the upper-left corner of the title text.

``` r
add_title <- function(.p, .title, .mult = list(x = 1, y = 1), .data = NULL) {
    
    if (missing(.title)) stop("Why is .title missing?")
    if (is.null(.title)) return(.p)
    
    x_range <- ggplot_build(.p)$layout$panel_ranges[[1]]$x.range
    y_range <- ggplot_build(.p)$layout$panel_ranges[[1]]$y.range
    if (!is.null(.p$coordinates$trans$x)) {
        x_range <- .p$coordinates$trans$x$inverse(x_range)
    }
    if (!is.null(.p$coordinates$trans$y)) {
        y_range <- .p$coordinates$trans$y$inverse(y_range)
    }
    min_x <- min(x_range) + .mult$x * 0.02 * diff(x_range)
    max_y <- max(y_range) - .mult$y * 0.02 * diff(y_range)
    .p <- .p +
        geom_text(data = .data, 
                  label = .title,
                  x = min_x, y = max_y,
                  color = 'black',
                  hjust = 0, vjust = 1, 
                  size = 14 * (25.4/72), # <-- 25.4/72 is to convert from mm to pt
                  fontface = 'plain')
    return(.p)
}
```

Individual plots for models by species only
===========================================

Function to create base plots
-----------------------------

This creates the base plots for those models that have only clade on the x-axis (i.e., those organized by species only—not by intestinal segment). All these models include log(mass) as a covariate so are plotted with log(mass) on the x-axis.

Function to create each plot depending on the `phylolm` model (`.model`), y-axis title (`y_axis_title`), plot title (`plot_title`), y-axis break points (`y_breaks`), y-axis labels (`y_labels`), and y-axis limits (`y_limits`). Only the first two are required.

``` r
clade_only_plot <- function(.model, 
                            y_axis_title, 
                            plot_title = NULL, 
                            y_breaks = ggplot2::waiver(), 
                            y_labels = ggplot2::waiver(), 
                            y_limits = NULL) {
    
    if (!'log_mass' %in% colnames(.model$X)) {
        stop("Model must include log(mass) as a covariate")
    }
    y_name <- {paste(.model$formula) %>% discard(~ grepl('~', .x))}[1]
    # Base plot
    .p <- predict_ci(.model) %>%
        # Exponentiate variables for plotting so we can more transparently
        # provide axis labels of the non-transformed numbers.
        # The same is done below when plotting the raw data.
        mutate_at(vars(log_mass, estimate, low, high), exp) %>% 
        ggplot(aes(log_mass, estimate, color = clade)) +
        # 95% CI envelopes
        geom_ribbon(aes(ymin = low, ymax = high, group = clade), 
                    fill = 'gray80', color = NA, alpha = 0.5) +
        # Raw data points
        geom_point(data = data_frame(estimate = as.numeric(.model$y), 
                                     log_mass = as.numeric(.model$X[,'log_mass']),
                                     clade = factor(as.integer(.model$X[,'cladeBat']), 
                                                    levels = c(0,1), 
                                                    labels = c('Rodent', 'Bat'))) %>%
                       mutate_at(vars(log_mass, estimate), exp),
                   size = 2) +
        # Regression fit
        geom_line() +
        theme(legend.position = 'none', legend.title = element_blank()) +
        # All these plots include log(mass) on the x-axis so I'm including this here.
        scale_x_continuous('Body mass (g)', trans = 'log', breaks = 10 * 4^(0:2),
                           limits = exp(c(2.043877, 5.199205)), expand = c(0,0)) +
        scale_y_continuous(y_axis_title, trans = 'log', breaks = y_breaks, 
                           labels = y_labels, limits = y_limits) +
        scale_color_manual(values = c('#9ecae1', '#de2d26'))
    
    # Add title if one is provided
    .p <- add_title(.p, plot_title)
    
    return(.p)
}
```

Creating plot objects
---------------------

For all the plots below...

1.  Both axes are on the log scale.
2.  Envelopes represent 95% CI for model predictions via parametric bootstrapping.

``` r
# Figure 1A
fig1a <- clade_only_plot(models$spp$log_intestinal_length,
                         "Intestinal length (cm)", 
                         y_breaks = 15 * 2^(0:2), 
                         plot_title = 'A') +
    theme(legend.position = c(0.05, 1), legend.justification = c(0, 1),
          axis.title.x = element_blank(), axis.text.x = element_blank())
# Figure 1B
fig1b <- clade_only_plot(models$spp$log_nsa,
                         expression("NSA (" * cm^2 * ")"), 
                         y_breaks = 5 * 2^(0:3), 
                         plot_title = 'B')
# Figure 4
fig4 <- clade_only_plot(models$spp$log_vill_surface_area,
                        expression("Villous surface area (" * cm^2 * ")"),
                        y_breaks = 50 * 3^(0:2)) +
    theme(legend.position = c(0.05, 1), legend.justification = c(0, 1))

# Figure 6
fig6 <- clade_only_plot(models$spp$log_total_enterocytes,
                        expression("Total enterocytes" %*% 10^{-9}),
                        # CHANGING UNITS HERE (from enterocytes to 
                        # 1e9 enterocytes:
                        y_breaks = 200e6 * 2^(0:3), y_labels = 0.2 * 2^(0:3)) +
    theme(legend.position = c(0.05, 1), legend.justification = c(0, 1))
```

Individual plots for models by species and intestinal segment
=============================================================

Objects to create base plots
----------------------------

Making data frame of confidence intervals. (Nesting by parameter, not position, bc the former is how they'll be plotted.)

``` r
pos_ci <- lapply(names(models$pos$prox), 
                 function(n) {
                     bind_rows(list(predict_ci(models$pos$prox[[n]]),
                                    predict_ci(models$pos$med[[n]]),
                                    predict_ci(models$pos$dist[[n]])))
                 }) %>% 
    bind_rows %>% 
    mutate(pos = factor(pos, levels = c('prox', 'med', 'dist'), 
                        labels = c('Proximal', 'Medial', 'Distal')))
```

Table of y-axis names for each parameter:

``` r
plot_names <- rbind(c("log_intestinal_diameter", "Intestinal ~ diameter ~ '(cm)'"),
                    c("log_villus_height", "Villus ~ height ~ '(mm)'"),
                    c("villus_width", "Villus ~ width ~ '(mm)'"),
                    c("crypt_width", "Crypt ~ width ~ '(mm)'"),
                    c("log_sef", "Surface ~ enlargement ~ factor ~ '(SEF)'"),
                    c("enterocyte_diameter", "Enterocyte ~ diameter ~ '(\u03BCm)'"),
                    c("log_enterocyte_density", 
                      "Enterocyte ~ density  %*% 10^{-6} ~ '(' * cm^{-2} * ')'")) %>% 
    as_tibble %>% 
    rename(og = V1, new = V2)
```

Function to create each plot depending on the input measurement name (`.measure`), custom y-axis break points (`y_breaks`) and labels (`y_labels`), and plot title (`plot_title`).

``` r
clade_pos_plot <- function(.measure, 
                           y_breaks = ggplot2::waiver(), 
                           y_labels = ggplot2::waiver(), 
                           plot_title = NULL) {
    # Getting model-prediction data frame for all three intestinal segments
    predict_df <- lapply(c('prox', 'med', 'dist'), 
                         function(seg_) {
                             df_ <- predict_ci(models$pos[[seg_]][[.measure]])
                             # Adding log_mass to models that don't include it
                             if (nrow(df_) == 2) {
                                 lm_ <- data$pos %>% 
                                     filter(grepl(seg_, pos, ignore.case = TRUE),
                                            measure == .measure)
                                 lm_ <- sapply(c('Bat', 'Rodent'),
                                               function(c_) {
                                                   range(lm_$log_mass[lm_$clade == c_])
                                               })
                                 df_ <- bind_rows(df_, df_) %>% 
                                     mutate(log_mass = c(t(lm_)),
                                            # signif = whether log(mass) was included 
                                            # in model (bc its P < 0.05)
                                            signif = 'no')
                             } else {
                                 df_ <- df_ %>% 
                                     mutate(signif = 'yes')
                             }
                             return(df_)
                         }) %>% 
        bind_rows %>% 
        mutate(signif = factor(signif, levels = c('no', 'yes')),
               pos = factor(pos, levels = c('prox', 'med', 'dist'), 
                            labels = c('Proximal', 'Medial', 'Distal')))
    
    # Is the y-axis log-transformed?
    y_logged <- grepl('log', .measure)

    # Exponentiating log(mass) so it's more transparent to create
    # axis labels on the non-transformed scale.
    # Doing this also for y variable if y_logged == TRUE
    if (y_logged) {
        predict_df <- predict_df %>% 
            mutate_at(vars(log_mass, estimate, low, high), exp)
        raw_data <- data$pos %>% 
            filter(measure == .measure) %>% 
            mutate_at(vars(log_mass, value), exp)
    } else {
        predict_df <- predict_df %>% 
            mutate_at(vars(log_mass), exp)
        raw_data <- data$pos %>% 
            filter(measure == .measure) %>% 
            mutate_at(vars(log_mass), exp)
    }
    
    # Provide y-axis breaks if none provided
    if (is(y_breaks, "waiver") & y_logged) {
        y_breaks <- predict_df %>% 
            summarize(min = min(low), max = max(high)) %>% 
            {seq(.$min, .$max, length.out = 4)} %>% 
            signif(digits = 3)
    }

    # y-axis transformation, depending on whether it's logged
    y_trans <- ifelse(y_logged, "log", "identity")
    
    # Title from `plot_names` data frame
    y_axis_title <- eval(parse(text = plot_names[plot_names$og == .measure,]$new))
    
    .p <- predict_df %>%
        ggplot(aes(log_mass, estimate, color = clade)) +
        geom_ribbon(aes(ymin = low, ymax = high, group = clade), 
                    fill = 'gray80', color = NA, alpha = 0.5) + 
        geom_line(aes(linetype = signif)) +
        geom_point(data = raw_data, aes(y = value), size = 2) +
        facet_wrap(~ pos) +
        scale_linetype_manual(values = c(2, 1), guide = FALSE, drop = FALSE) +
        scale_x_continuous('Body mass (g)', trans = 'log', breaks = 10 * 4^(0:2),
                           limits = exp(c(2.043877, 5.199205)), expand = c(0,0)) +
        scale_y_continuous(y_axis_title, trans = y_trans, breaks = y_breaks, 
                           labels = y_labels) +
        scale_color_manual(values = c('#9ecae1', '#de2d26')) +
        theme(legend.position = 'none', 
              axis.line = element_line(colour = "black", size = 0.5))
    
    # Add title if it's provided.
    # .mult is set to 3 for the x-axis to compensate for 3 facets
    # .data is provided to make sure it only shows up in the first facet
    .p <- add_title(.p, .title = plot_title, .mult = list(x = 3, y = 1),
                    .data = data_frame(pos = sort(unique(predict_df$pos))[1]))
    
    return (.p)
}
```

Creating plot objects
---------------------

X-axis is on log scale for all plots.

Y-axis is on log scale for all plots *except* the following:

-   fig2b
-   fig2c
-   fig5a

``` r
# Figure 1c
fig1c <- clade_pos_plot('log_intestinal_diameter', y_breaks = 0.3 * 2^(0:2),
                        plot_title = 'C')

# Figure 2a
fig2a <- clade_pos_plot('log_villus_height', y_breaks = 0.2 * 2 ^(0:2),
                        plot_title = 'A') +
        theme(legend.position = 'top', legend.margin = margin(0,0,0,0),
              axis.text.x = element_blank(), axis.title.x = element_blank())

# Figure 2b
fig2b <- clade_pos_plot('villus_width', y_breaks = seq(0.04, 0.12, 0.04),
                        plot_title = 'B') +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          strip.background = element_blank(), strip.text = element_blank())

# Figure 2c
fig2c <- clade_pos_plot('crypt_width', y_breaks = c(0.02, 0.04),
                        plot_title = 'C') + 
    theme(strip.background = element_blank(), strip.text = element_blank())

# figure 3
fig3 <- clade_pos_plot('log_sef', y_breaks = 5 * 2^(0:2)) +
    theme(legend.position = 'top')

# Figure 5a
fig5a <- clade_pos_plot('enterocyte_diameter',
               # CHANGING UNITS HERE (from mm to µm):
               y_breaks = seq(6e-3, 10e-3, 2e-3), y_labels = seq(6, 10, 2),
               plot_title = 'A') +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())

# Figure 5b
fig5b <- clade_pos_plot('log_enterocyte_density',
               # CHANGING UNITS HERE (from enterocytes / cm^2 to 
               # 1e6 enterocytes / cm^2):
               y_breaks = 8e6 * 3^(0:2), y_labels = 8 * 3^(0:2),
               plot_title = 'B') +
    theme(strip.background = element_blank(), strip.text = element_blank())
```

Individual plots for clearance and absorption
=============================================

All axes are on the log scale for all plots.

``` r
clear_label <- bquote(.("L-arabinose clearance (\u03BC") * 
                          l ~ min^{-1} ~ cm^{-2} * ")")
absorp_label <- bquote(.("Absorption") %*% 10^{3} ~ .("/ ") ~
                                .("intest. area (") * cm^{-2} * .(")"))

fig7a <- data$clear %>%
    mutate(sef = exp(log_sef), clear = exp(log_clear)) %>% 
    ggplot(aes(sef, clear, shape = diet, color = clade)) +
    geom_point(size = 3, fill = 'gray60') +
    theme(legend.position = 'top',
          legend.title = element_text(size = 10, face = 'bold.italic'),
          legend.margin = margin(0,0,4,0),
          legend.justification = c(1, 0),
          legend.key.size = unit(0.75, 'lines')) +
    guides(shape = 'none', color = guide_legend('Clade:')) +
    scale_shape_manual(values = c(16, 17)) +
    scale_color_manual(values = c('#9ecae1', '#de2d26')) +
    scale_x_continuous("Surface enlagement factor (SEF)", 
                       trans = 'log', breaks = 8 * 1.5^(0:2)) +
    scale_y_continuous(clear_label,
                       trans = 'log', breaks = 1 * 3^(0:2))
fig7a <- fig7a %>%
    add_title('A')


fig7b <- data$clear %>%
    filter(!is.na(log_enterocyte_density)) %>% 
    mutate_at(vars(log_enterocyte_density, log_clear), exp) %>% 
    ggplot(aes(log_enterocyte_density, log_clear, shape = diet, color = clade)) +
    geom_point(size = 3, fill = 'gray60') +
    theme(legend.position = 'top',
          legend.title = element_text(size = 10, face = 'bold.italic'),
          legend.margin = margin(0,0,4,0),
          legend.justification = c(-0.05, 0),
          legend.key.size = unit(0.75, 'lines'),
          axis.title.y = element_blank(), axis.text.y = element_blank()) +
    guides(shape = guide_legend('Diet:'), color = 'none') +
    scale_shape_manual(values = c(16, 17)) +
    scale_color_manual(values = c('#9ecae1', '#de2d26')) +
    scale_x_continuous(eval(parse(
        text = plot_names[plot_names$og == 'log_enterocyte_density',]$new)), 
                       trans = 'log', breaks = 1e7 * 3^(0:2), labels = 1e1 * 3^(0:2)) +
    scale_y_continuous(clear_label,
                       trans = 'log', breaks = 1 * 3^(0:2))
fig7b <- fig7b %>% 
    add_title('B')


fig7c <- clade_only_plot(models$absorp,
                         y_axis_title = absorp_label,
                         plot_title = "C", 
                         # CHANGING UNITS HERE:
                         y_breaks = 0.001 * 4^(0:2), y_labels = 1 * 4^(0:2),
                         y_limits = c(7e-4, 0.02))



fig7d <- data$absorp %>%
    mutate(diet = sapply(species, function(s) data$spp$diet[data$spp$species == s]),
           diet = factor(ifelse(diet == "Protein", paste(diet), "Carb"),
                         levels = c('Carb', 'Protein'))) %>% 
    mutate_at(vars(log_total_enterocytes, log_absorp), exp) %>% 
    ggplot(aes(log_total_enterocytes, log_absorp, shape = diet, color = clade)) +
    geom_point(size = 3, fill = 'gray60') +
    theme(legend.position = 'none', 
          axis.title.y = element_blank(), axis.text.y = element_blank()) +
    scale_shape_manual(values = c(16, 17)) +
    scale_color_manual(values = c('#9ecae1', '#de2d26')) +
    scale_x_continuous(expression("Total enterocytes" %*% 10^{-9}),
                       trans = 'log', breaks = 200e6 * 2^(0:3), labels = 0.2 * 2^(0:3)) +
    scale_y_continuous(absorp_label, trans = 'log', limits = c(7e-4, 0.0172),
                       # CHANGING UNITS HERE:
                       breaks = 0.001 * 4^(0:2), labels = 1 * 4^(0:2))
fig7d <- fig7d %>% 
    add_title('D')
```

Creating and saving final plots
===============================

Function to combine plots
-------------------------

``` r
# Combining multiple figures using the naming scheme `figX` where `X` is the 
# figure number I'm interested in plotting
# It also works for single figures.
one_fig <- function(fig_num, fig1_heights = c(7.9, 9, 10)) {
    grob_list <- lapply(ls(envir = .GlobalEnv)[grepl(paste0('^fig', fig_num), 
                                                     ls(envir = .GlobalEnv))], 
                        function(n) {
                            ggplotGrob(eval(parse(text = n)))
                        })
    if (fig_num == 1) {
        grid.draw(arrangeGrob(grobs = grob_list, ncol = 1, nrow = length(grob_list), 
                              heights = fig1_heights))
    } else if (fig_num %in% c(2, 5)) {
        grob_list <- lapply(grob_list, 
               function(gg) {
                   colnames(gg) <- paste0(seq_len(ncol(gg)))
                   return(gg)
               })
        grid.draw(do.call(gtable_rbind, grob_list))
    } else if (fig_num == 7) {
        grid.draw(rbind(do.call(cbind, c(grob_list[1:2], size = 'last')),
                        do.call(cbind, c(grob_list[3:4], size = 'last')),
                        size = 'last'))
    } else {
        grid.draw(grob_list[[1]])
    }
    invisible()
}

# Employs the above function, plus saves the output
save_fig <- function(fig_num, fig1_heights = c(7.9, 9, 10), .seed = NULL, ...) {
    file_name <- sprintf('figs/fig%02d.pdf', fig_num)
    if (!is.null(.seed)) set.seed(.seed)
    quartz(type = 'pdf', file = file_name, family = 'Helvetica', ...)
    one_fig(fig_num, fig1_heights)
    invisible(dev.off())
}
```

``` r
save_fig(1, width = 3.875, height = 3.125 * 3, .seed = 1)
save_fig(2, width = 3.875, height = 3.125 * 3, .seed = 2)
save_fig(3, width = 3.875, height = 3.125, .seed = 3)
save_fig(4, width = 3.875, height = 3.125, .seed = 4)
save_fig(5, width = 3.875, height = 3.125 * 2, .seed = 5)
save_fig(6, width = 3.875, height = 3.125, .seed = 6)
save_fig(7, width = 3.125 * 2, height = 3.125 * 2, .seed = 7)
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
    ##  date     2017-12-08

    ## Packages -----------------------------------------------------------------

    ##  package     * version date       source        
    ##  ape         * 5.0     2017-10-30 CRAN (R 3.4.2)
    ##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.0)
    ##  backports     1.1.1   2017-09-25 CRAN (R 3.4.2)
    ##  base        * 3.4.2   2017-10-04 local         
    ##  bindr         0.1     2016-11-13 CRAN (R 3.4.0)
    ##  bindrcpp    * 0.2     2017-06-17 CRAN (R 3.4.0)
    ##  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.0)
    ##  commonmark    1.4     2017-09-01 CRAN (R 3.4.1)
    ##  compiler      3.4.2   2017-10-04 local         
    ##  corphyloCpp * 1.0     <NA>       local         
    ##  datasets    * 3.4.2   2017-10-04 local         
    ##  devtools      1.13.3  2017-08-02 CRAN (R 3.4.1)
    ##  digest        0.6.12  2017-01-27 CRAN (R 3.4.0)
    ##  dplyr       * 0.7.4   2017-09-28 CRAN (R 3.4.2)
    ##  evaluate      0.10.1  2017-06-24 CRAN (R 3.4.1)
    ##  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.4.0)
    ##  glue          1.2.0   2017-10-29 CRAN (R 3.4.2)
    ##  graphics    * 3.4.2   2017-10-04 local         
    ##  grDevices   * 3.4.2   2017-10-04 local         
    ##  grid        * 3.4.2   2017-10-04 local         
    ##  gridExtra   * 2.3     2017-09-09 CRAN (R 3.4.1)
    ##  gtable        0.2.0   2016-02-26 CRAN (R 3.4.0)
    ##  hms           0.3     2016-11-22 CRAN (R 3.4.0)
    ##  htmltools     0.3.6   2017-04-28 cran (@0.3.6) 
    ##  knitr         1.17    2017-08-10 CRAN (R 3.4.1)
    ##  lattice       0.20-35 2017-03-25 CRAN (R 3.4.2)
    ##  lazyeval      0.2.1   2017-10-29 CRAN (R 3.4.2)
    ##  magrittr      1.5     2014-11-22 CRAN (R 3.4.0)
    ##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods     * 3.4.2   2017-10-04 local         
    ##  munsell       0.4.3   2016-02-13 CRAN (R 3.4.0)
    ##  nlme          3.1-131 2017-02-06 CRAN (R 3.4.2)
    ##  parallel      3.4.2   2017-10-04 local         
    ##  phylolm     * 2.5     2016-10-17 CRAN (R 3.4.0)
    ##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.0)
    ##  plyr          1.8.4   2016-06-08 CRAN (R 3.4.0)
    ##  purrr       * 0.2.4   2017-10-18 CRAN (R 3.4.2)
    ##  R6            2.2.2   2017-06-17 CRAN (R 3.4.0)
    ##  Rcpp          0.12.13 2017-09-28 CRAN (R 3.4.2)
    ##  readr       * 1.1.1   2017-05-16 CRAN (R 3.4.0)
    ##  rlang         0.1.4   2017-11-05 CRAN (R 3.4.2)
    ##  rmarkdown     1.6     2017-06-15 CRAN (R 3.4.0)
    ##  roxygen2      6.0.1   2017-02-06 CRAN (R 3.4.0)
    ##  rprojroot     1.2     2017-01-16 cran (@1.2)   
    ##  scales        0.5.0   2017-08-24 CRAN (R 3.4.1)
    ##  stats       * 3.4.2   2017-10-04 local         
    ##  stringi       1.1.5   2017-04-07 CRAN (R 3.4.0)
    ##  stringr       1.2.0   2017-02-18 CRAN (R 3.4.0)
    ##  tibble        1.3.4   2017-08-22 CRAN (R 3.4.1)
    ##  tidyr       * 0.7.2   2017-10-16 CRAN (R 3.4.2)
    ##  tidyselect    0.2.3   2017-11-06 CRAN (R 3.4.2)
    ##  tools         3.4.2   2017-10-04 local         
    ##  utils       * 3.4.2   2017-10-04 local         
    ##  withr         2.1.0   2017-11-01 CRAN (R 3.4.2)
    ##  xml2          1.1.1   2017-01-24 CRAN (R 3.4.0)
    ##  yaml          2.1.14  2016-11-12 cran (@2.1.14)
