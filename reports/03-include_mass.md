Test whether to include log\_mass in phylogenetic linear regression
================
Lucas Nell
11 Sep 2018

This file determines whether to use `log_mass` in `phylolm` regressions
and `cor_phylo` correlations. These analyses take about half an hour to
run.

All p-values below are for whether the coefficient for log(mass) are not
equal to zero.

# Necessary data:

``` r
# Morphometrics by species
spp_df <- get_df('spp')
tr <- get_tr('spp')
# Morphometrics by species and segment
seg_types <- c('prox', 'med', 'dist')
pos_df_list <- lapply(seg_types, get_df, .df = 'pos')
names(pos_df_list) <- seg_types
# Absorption by species
absorp_df <- get_df('absorp')
absorp_se_df <- get_df('absorp', .stat = 'se')  # <-- contains standard errors
absorp_tr <- get_tr('absorp')
# Clearance by species
clear_df <- get_df('clear')
clear_se_df <- get_df('clear', .stat = 'se')  # <-- contains standard errors
clear_tr <- get_tr('clear')
```

# `SEF` on `Diet`

`phylolm` call and output:

``` r
set.seed(29851644)
diet_fit <- phylolm(log_sef ~ diet + log_mass, data = spp_df, phy = tr,
                    model = 'lambda', boot = 2000)
```

P-value:

    ## P = 0.501

# `Absorption` on `Clade`

`phylolm` call and output:

``` r
set.seed(1092141389)
absorp_fit <- suppressWarnings(  # gives warning about lambda being very low
    phylolm(log_absorp ~ clade + log_mass, data = absorp_df, phy = absorp_tr, 
            model = 'lambda', boot = 2000)
)
```

P-value:

    ## P = 0.000

# `Morphometrics` on `Clade`

`phylolm` calls:

``` r
spp_ys <- c("log_intestinal_length", "log_nsa", "log_vill_surface_area",
            "log_total_enterocytes")
set.seed(357885189)
spp_fits <- lapply(
    spp_ys,
    function(y) {
        f <- paste(y, '~ clade + log_mass')
        suppressWarnings(
            do.call("phylolm", list(as.formula(f), data = as.name("spp_df"),
                                    phy = as.name("tr"), model = 'lambda',
                                    boot = 2000))
        )
    })
names(spp_fits) <- spp_ys
```

P-values:

    ##                           P
    ## log_intestinal_length 0.028
    ## log_nsa               0.000
    ## log_vill_surface_area 0.000
    ## log_total_enterocytes 0.000

# `Morphometrics` on `Clade`, separately by segment

`phylolm`
calls:

``` r
pos_ys <- c('log_intestinal_diameter', 'log_villus_height', 'villus_width', 
            'crypt_width', 'log_sef', 'enterocyte_diameter', 'log_enterocyte_density')
set.seed(632929430)
pos_fits <- lapply(
    seg_types,
    function(pos) {
        # Assigning to obj named <pos>_df so that the call identifies the position
        assign(paste0(pos, '_df'), pos_df_list[[pos]])
        lapply(
            pos_ys,
            function(y) {
                f <- paste(y, ' ~ clade + log_mass')
                arg_list <- list(
                    formula = as.formula(f),
                    data = as.name(paste0(pos, "_df")),
                    phy = as.name("tr"), model = 'lambda')
                # Some models don't find the peak likelihood unless specifying a
                # starting value of 0.1.
                LL_nostart <- suppressWarnings(do.call("phylolm", arg_list))$logLik
                arg_list_wstart <- c(arg_list, starting.value = quote(list(lambda = 0.1)))
                LL_wstart <- suppressWarnings(do.call("phylolm", arg_list_wstart))$logLik
                if (LL_wstart > LL_nostart) arg_list <- arg_list_wstart
                # Adding number of bootstrap replicates
                arg_list <- c(arg_list, boot = 2000)
                # Now call phylolm
                suppressWarnings(do.call("phylolm", arg_list))
            })
    })
names(pos_fits) <- seg_types
for (i in 1:length(pos_fits)) names(pos_fits[[i]]) <- pos_ys; rm(i)
```

Because there are so many models, I’m writing whether to include mass to
a csv file:

``` r
include_df <- lapply(
    seg_types, 
    function(p) {
        data_frame(pos = p, y = names(pos_fits[[p]]),
                   include = sapply(pos_fits[[p]], pval, 
                                    params = 'log_mass') < 0.05)
    }) %>% 
    bind_rows
write_csv(include_df, 'output/include_mass_pos.csv')
```

P-values:

    ## proximal:

    ##                             P
    ## log_intestinal_diameter 0.000
    ## log_villus_height       0.002
    ## villus_width            0.009
    ## crypt_width             0.881
    ## log_sef                 0.022
    ## enterocyte_diameter     0.823
    ## log_enterocyte_density  0.203

    ## medial:

    ##                             P
    ## log_intestinal_diameter 0.006
    ## log_villus_height       0.039
    ## villus_width            0.062
    ## crypt_width             0.379
    ## log_sef                 0.322
    ## enterocyte_diameter     0.621
    ## log_enterocyte_density  0.208

    ## distal:

    ##                             P
    ## log_intestinal_diameter 0.000
    ## log_villus_height       0.032
    ## villus_width            0.193
    ## crypt_width             0.500
    ## log_sef                 0.329
    ## enterocyte_diameter     0.782
    ## log_enterocyte_density  0.462

# `Clearance` and `SEF`

In `cor_phylo`, you can include covariate(s) that might be having
confounding effects.

Below, I’m trying out whether either clearance or SEF needs body mass as
a covariate by including it separately for clearance, then for SEF.

``` r
# Making data frame with all info, including measurement errors:
clear_df_bm <- clear_df
clear_df_bm$log_clear_se <- clear_se_df$log_clear
clear_df_bm$log_sef_se <- clear_se_df$log_sef
# For this comparison, I have to remove one row that doesn't have body mass
clear_df_bm <- clear_df_bm[!is.na(clear_df$log_mass),]

# Phylogenetic tree, removing species with no body mass data
clear_tr_bm <- ape::drop.tip(
    clear_tr,
    tip = clear_tr$tip.label[!clear_tr$tip.label %in% clear_df_bm$species])

# Options that are the same for both calls
opts <- list(traits = quote(list(log_clear, log_sef)),
             meas_errors = quote(list(log_clear_se, log_sef_se)),
             species = quote(species), data = quote(clear_df_bm),
             phy = quote(clear_tr_bm), boot = 2000,
             method = "nelder-mead-r",
             constrain_d = TRUE,
             max_iter = 1e6, keep_boots = "fail")
# Now separating into different lists:
opts <- list(clear = c(opts, covariates = quote(list(log_clear = log_mass))),
             sef = c(opts, covariates = quote(list(log_sef = log_mass))))

# cor_phylo run with bootstrapping (takes ~1.5 min)
set.seed(1844365955)
clear_sef <- map(c("clear", "sef"), ~ do.call(cor_phylo, opts[[.x]])) %>%
    setNames(c("clear", "sef"))
```

Some bootstrap replicates did not converge, so I’m refitting them using
a higher threshold for the reciprocal condition number of two matrices
inside the likelihood function. This makes the optimization process more
strongly “bounce away” from badly conditioned matrices. From trial and
error, two sets of refits (using `rcond_threshold` values of `5e-4` and
`5e-3`) seem to make all the replicates converge and provide sensible
results.

``` r
cp_boot_refits <- list(
    clear = list(
        one = refit_boots(clear_sef$clear, rcond_threshold = 5e-4),
        two = NA
    ),
    sef = list(
        one = refit_boots(clear_sef$sef, rcond_threshold = 5e-4),
        two = NA
    ))
cp_boot_refits$clear$two <- refit_boots(clear_sef$clear,
                                inds = which(map_lgl(cp_boot_refits$clear$one,
                                                     ~ .x$convcode != 0)),
                                rcond_threshold = 5e-3)
cp_boot_refits$sef$two <- refit_boots(clear_sef$sef,
                              inds = which(map_lgl(cp_boot_refits$sef$one,
                                                   ~ .x$convcode != 0)),
                              rcond_threshold = 5e-3)
```

P-values for including body mass for clearance and SEF, respectively:

    ## P = 0.998

    ## P = 0.676

# Session info

This outlines the package versions I used for this
    script.

    ## Session info -------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.5.1 (2018-07-02)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2018-09-11

    ## Packages -----------------------------------------------------------------

    ##  package      * version date       source         
    ##  ape          * 5.1     2018-04-04 CRAN (R 3.5.0) 
    ##  assertthat     0.2.0   2017-04-11 CRAN (R 3.5.0) 
    ##  backports      1.1.2   2017-12-13 CRAN (R 3.5.0) 
    ##  base         * 3.5.1   2018-07-05 local          
    ##  bindr          0.1.1   2018-03-13 CRAN (R 3.5.0) 
    ##  bindrcpp     * 0.2.2   2018-03-29 CRAN (R 3.5.0) 
    ##  codetools      0.2-15  2016-10-05 CRAN (R 3.5.1) 
    ##  compiler       3.5.1   2018-07-05 local          
    ##  crayon         1.3.4   2017-09-16 CRAN (R 3.5.0) 
    ##  datasets     * 3.5.1   2018-07-05 local          
    ##  devtools       1.13.6  2018-06-27 CRAN (R 3.5.0) 
    ##  digest         0.6.16  2018-08-22 CRAN (R 3.5.0) 
    ##  dplyr        * 0.7.6   2018-06-29 CRAN (R 3.5.1) 
    ##  evaluate       0.11    2018-07-17 CRAN (R 3.5.0) 
    ##  future         1.9.0   2018-07-23 CRAN (R 3.5.0) 
    ##  future.apply   1.0.1   2018-08-26 CRAN (R 3.5.0) 
    ##  globals        0.12.2  2018-08-25 CRAN (R 3.5.0) 
    ##  glue           1.3.0   2018-07-17 CRAN (R 3.5.0) 
    ##  graphics     * 3.5.1   2018-07-05 local          
    ##  grDevices    * 3.5.1   2018-07-05 local          
    ##  grid           3.5.1   2018-07-05 local          
    ##  hms            0.4.2   2018-03-10 CRAN (R 3.5.0) 
    ##  htmltools      0.3.6   2017-04-28 CRAN (R 3.5.0) 
    ##  knitr          1.20    2018-02-20 CRAN (R 3.5.0) 
    ##  lattice        0.20-35 2017-03-25 CRAN (R 3.5.1) 
    ##  listenv        0.7.0   2018-01-21 CRAN (R 3.5.0) 
    ##  magrittr       1.5     2014-11-22 CRAN (R 3.5.0) 
    ##  Matrix         1.2-14  2018-04-13 CRAN (R 3.5.1) 
    ##  memoise        1.1.0   2017-04-21 CRAN (R 3.5.0) 
    ##  methods      * 3.5.1   2018-07-05 local          
    ##  nlme           3.1-137 2018-04-07 CRAN (R 3.5.1) 
    ##  nloptr         1.0.4   2017-08-22 CRAN (R 3.5.0) 
    ##  parallel       3.5.1   2018-07-05 local          
    ##  phylolm      * 2.6     2018-05-31 CRAN (R 3.5.0) 
    ##  phyr         * 0.1.5   2018-09-11 local          
    ##  pillar         1.3.0   2018-07-14 CRAN (R 3.5.0) 
    ##  pkgconfig      2.0.2   2018-08-16 CRAN (R 3.5.0) 
    ##  purrr        * 0.2.5   2018-05-29 CRAN (R 3.5.0) 
    ##  R6             2.2.2   2017-06-17 CRAN (R 3.5.0) 
    ##  Rcpp           0.12.18 2018-07-23 cran (@0.12.18)
    ##  readr        * 1.1.1   2017-05-16 CRAN (R 3.5.0) 
    ##  rlang          0.2.2   2018-08-16 CRAN (R 3.5.0) 
    ##  rmarkdown      1.10    2018-06-11 CRAN (R 3.5.0) 
    ##  rprojroot      1.3-2   2018-01-03 CRAN (R 3.5.0) 
    ##  stats        * 3.5.1   2018-07-05 local          
    ##  stringi        1.2.4   2018-07-20 CRAN (R 3.5.0) 
    ##  stringr        1.3.1   2018-05-10 CRAN (R 3.5.0) 
    ##  tibble         1.4.2   2018-01-22 CRAN (R 3.5.0) 
    ##  tidyr        * 0.8.1   2018-05-18 CRAN (R 3.5.0) 
    ##  tidyselect     0.2.4   2018-02-26 CRAN (R 3.5.0) 
    ##  tools          3.5.1   2018-07-05 local          
    ##  utils        * 3.5.1   2018-07-05 local          
    ##  withr          2.1.2   2018-03-15 CRAN (R 3.5.0) 
    ##  yaml           2.2.0   2018-07-25 CRAN (R 3.5.0)
