Visualizing phylogenetic tree
================
Lucas Nell
11 Sep 2017

Visualizing tree
================

Here is how I created the phylogenetic tree.

``` r
# Commented areas are for when I was basing tip points on `sp_df` data frame values
x_end <- 32
gg_tr <- ggtree(tr)
gg_tr$data$x <- gg_tr$data$x - max(gg_tr$data$x)
phylo <- gg_tr + #%<+% {sp_df %>% select(species, everything())} +
    theme_tree2(axis.title.x = element_text(size = 14),
                legend.position = c(0.25, 0.75), legend.box = 'horizontal',
                legend.background = element_rect(color = NA, fill = NA)) +
    geom_tiplab(aes(x = x + 3), size = 2, fontface = 'bold.italic') +
    geom_tippoint(aes(x = x)) +  # , size = log_mass, color = diet)) +
    geom_text(data = data_frame(x = rep(x_end + 1, 2), y = c(14, 5), 
                                label = c('Rodents', 'Bats')), 
              aes(label = label), angle = -90, vjust = 0, size = 6) +
    geom_segment(data = data_frame(x = rep(x_end, 2), xend = rep(x_end, 2), 
                                   y = c(1, 10), yend = c(9, 18)), 
              aes(xend = xend, yend = yend)) +
    scale_x_continuous('Time (mya)', limits = c(-100, x_end + 2),
                       breaks = seq(-100, 0, 25), labels = seq(100, 0, -25))# +
    # scale_color_manual('Diet', values = c('#e66101','#fdb863','#b2abd2','#5e3c99'),
    #                      guide = guide_legend(direction = "vertical",
    #                                             title.position = 'top', 
    #                                             title.hjust = 0.5)) +
    # scale_size_continuous('log(mass)', 
    #                      guide = guide_legend(direction = "horizontal",
    #                                           title.position = 'top', 
    #                                           title.hjust = 0.5))
```

Saving this figure.

``` r
ggsave('figs/phylo.pdf', phylo, width = 6, height = 4)
```

Session info
============

This outlines the package versions I used for this script.

    ## Session info -------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.4.1 (2017-06-30)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2017-09-11

    ## Packages -----------------------------------------------------------------

    ##  package    * version date       source        
    ##  ape        * 4.1     2017-02-14 CRAN (R 3.4.0)
    ##  assertthat   0.2.0   2017-04-11 CRAN (R 3.4.0)
    ##  backports    1.1.0   2017-05-22 CRAN (R 3.4.0)
    ##  base       * 3.4.1   2017-07-07 local         
    ##  bindr        0.1     2016-11-13 CRAN (R 3.4.0)
    ##  bindrcpp   * 0.2     2017-06-17 CRAN (R 3.4.0)
    ##  colorspace   1.3-2   2016-12-14 CRAN (R 3.4.0)
    ##  compiler     3.4.1   2017-07-07 local         
    ##  datasets   * 3.4.1   2017-07-07 local         
    ##  devtools     1.13.3  2017-08-02 CRAN (R 3.4.1)
    ##  digest       0.6.12  2017-01-27 CRAN (R 3.4.0)
    ##  dplyr      * 0.7.3   2017-09-09 CRAN (R 3.4.1)
    ##  evaluate     0.10.1  2017-06-24 CRAN (R 3.4.1)
    ##  ggplot2    * 2.2.1   2016-12-30 CRAN (R 3.4.0)
    ##  ggtree     * 1.8.1   2017-04-30 Bioconductor  
    ##  glue         1.1.1   2017-06-21 CRAN (R 3.4.0)
    ##  graphics   * 3.4.1   2017-07-07 local         
    ##  grDevices  * 3.4.1   2017-07-07 local         
    ##  grid         3.4.1   2017-07-07 local         
    ##  gtable       0.2.0   2016-02-26 CRAN (R 3.4.0)
    ##  hms          0.3     2016-11-22 CRAN (R 3.4.0)
    ##  htmltools    0.3.6   2017-04-28 cran (@0.3.6) 
    ##  jsonlite     1.5     2017-06-01 CRAN (R 3.4.0)
    ##  knitr        1.17    2017-08-10 CRAN (R 3.4.1)
    ##  lattice      0.20-35 2017-03-25 CRAN (R 3.4.1)
    ##  lazyeval     0.2.0   2016-06-12 CRAN (R 3.4.0)
    ##  magrittr     1.5     2014-11-22 CRAN (R 3.4.0)
    ##  memoise      1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods    * 3.4.1   2017-07-07 local         
    ##  munsell      0.4.3   2016-02-13 CRAN (R 3.4.0)
    ##  nlme         3.1-131 2017-02-06 CRAN (R 3.4.1)
    ##  parallel     3.4.1   2017-07-07 local         
    ##  pkgconfig    2.0.1   2017-03-21 CRAN (R 3.4.0)
    ##  plyr         1.8.4   2016-06-08 CRAN (R 3.4.0)
    ##  purrr        0.2.3   2017-08-02 CRAN (R 3.4.1)
    ##  R6           2.2.2   2017-06-17 CRAN (R 3.4.0)
    ##  Rcpp         0.12.12 2017-07-15 CRAN (R 3.4.1)
    ##  readr        1.1.1   2017-05-16 CRAN (R 3.4.0)
    ##  rlang        0.1.2   2017-08-09 CRAN (R 3.4.1)
    ##  rmarkdown    1.6     2017-06-15 CRAN (R 3.4.0)
    ##  rprojroot    1.2     2017-01-16 cran (@1.2)   
    ##  rvcheck      0.0.9   2017-07-10 CRAN (R 3.4.1)
    ##  scales       0.5.0   2017-08-24 CRAN (R 3.4.1)
    ##  stats      * 3.4.1   2017-07-07 local         
    ##  stringi      1.1.5   2017-04-07 CRAN (R 3.4.0)
    ##  stringr      1.2.0   2017-02-18 CRAN (R 3.4.0)
    ##  tibble       1.3.4   2017-08-22 CRAN (R 3.4.1)
    ##  tidyr        0.7.1   2017-09-01 CRAN (R 3.4.1)
    ##  tools        3.4.1   2017-07-07 local         
    ##  treeio     * 1.0.2   2017-05-01 Bioconductor  
    ##  utils      * 3.4.1   2017-07-07 local         
    ##  withr        2.0.0   2017-07-28 CRAN (R 3.4.1)
    ##  yaml         2.1.14  2016-11-12 cran (@2.1.14)
