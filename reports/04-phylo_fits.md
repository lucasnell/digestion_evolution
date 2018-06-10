Phylogenetic linear regressions and correlations
================
Lucas Nell
10 Jun 2018

This file conducts the linear regressions using `phylolm::phylolm` and
computes correlations using `phyr::cor_phylo`.

# Retrieve data

For more information on the functions `get_df` and `get_tr` below (plus
`filter_tr` and `cp_mat` used later), see
[`R/get_data.R`](R/get_data.R).

Note that absorption and clearance data need standard errors as well as
means.

``` r
# Tree for all morphometric and diet analyses
tr <- get_tr('spp')
# Morphometrics by species
spp_df <- get_df('spp')
# Morphometrics by species and segment
seg_types <- c('prox', 'med', 'dist')
pos_df_list <- lapply(seg_types, get_df, .df = 'pos')
names(pos_df_list) <- seg_types
# Absorption by species
absorp_df <- get_df('absorp')
absorp_se_df <- get_df('absorp', .stat = 'se')
absorp_tr <- get_tr('absorp')
# Clearance by species
clear_df <- get_df('clear')
clear_se_df <- get_df('clear', .stat = 'se')
clear_tr <- get_tr('clear')
```

# `phylolm`

The following sections are regressions using `phylolm::phylolm`.

## `SEF` on `Diet`

`phylolm` call:

``` r
set.seed(581120)
diet_fit <- phylolm(log_sef ~ diet, data = spp_df, phy = tr,
                    model = 'lambda', boot = 2000)
```

Saving output:

``` r
readr::write_rds(diet_fit, 'output/models_diet.rds')
```

## `Absorption` on `Clade`

“Absorption” here means `Fractional absorption / (total intestinal
surface)`, where `total intestinal surface = NSA * SEF`

`phylolm` call:

``` r
set.seed(454094511)
absorp_fit <- suppressWarnings(  # gives warning about lambda being very low
    phylolm(log_absorp ~ clade + log_mass, data = absorp_df, phy = absorp_tr, 
            model = 'lambda', boot = 2000)
)
```

Saving output:

``` r
readr::write_rds(absorp_fit, 'output/models_absorp.rds')
```

## `Morphometrics` on `Clade`

List of `Morphometrics`:

  - Intestinal length
  - NSA
  - Villus surface area
  - Total number of enterocytes (log-transformed)
      - Calculated as such: `log(NSA * enterocyte_density)`

> log(body mass) as covariate for all

These are the column names for the above parameters:

``` r
spp_ys <- c("log_intestinal_length", "log_nsa", "log_vill_surface_area",
            "log_total_enterocytes")
```

`phylolm` calls:

``` r
# (takes ~7.5 min)
set.seed(88754829)
spp_fits <- lapply(
    spp_ys,
    function(y_) {
        f <- paste(y_, '~ clade + log_mass')
        suppressWarnings(
            do.call("phylolm", list(as.formula(f), data = as.name("spp_df"),
                                    phy = as.name("tr"), model = 'lambda',
                                    boot = 2000))
        )
    })
names(spp_fits) <- spp_ys
```

Saving output:

``` r
readr::write_rds(spp_fits, 'output/models_spp.rds')
```

## `Morphometrics` on `Clade`, separately by segment

(Segment = proximal, medial, or distal)

List of `Y`s:

  - Intestinal diameter (log-transformed, body mass as covariate)
  - Villus height (body mass as covariate)
  - Villus width
  - Crypt width
  - Surface enlargement factor (SEF)
  - Enterocyte diameter
  - Enterocytes per cm^2 NSA (log-transformed)

Below are the column names for these parameters and all the segment
types.

``` r
pos_ys <- c('log_intestinal_diameter', 'log_villus_height', 'villus_width', 
            'crypt_width', 'log_sef', 'enterocyte_diameter', 'log_enterocyte_density')
seg_types <- c('prox', 'med', 'dist')
```

Below is a data frame including whether or not to include `log_mass` as
a covariate. This determination was based on whether `log_mass` had a
significant effect when it was included in the model, where p-values
were based on parametric bootstrapping (see
[`docs/03-include_mass`](docs/03-include_mass.md)).

``` r
include_mass <- read_csv('output/include_mass_pos.csv', col_types = 'ccl')
include_mass
```

    ## # A tibble: 21 x 3
    ##    pos   y                       include
    ##    <chr> <chr>                   <lgl>  
    ##  1 prox  log_intestinal_diameter TRUE   
    ##  2 prox  log_villus_height       TRUE   
    ##  3 prox  villus_width            TRUE   
    ##  4 prox  crypt_width             FALSE  
    ##  5 prox  log_sef                 TRUE   
    ##  6 prox  enterocyte_diameter     FALSE  
    ##  7 prox  log_enterocyte_density  FALSE  
    ##  8 med   log_intestinal_diameter TRUE   
    ##  9 med   log_villus_height       TRUE   
    ## 10 med   villus_width            FALSE  
    ## # ... with 11 more rows

`phylolm` call:

``` r
# (takes ~32 min)
set.seed(25413535)
pos_fits <- lapply(
    seg_types,
    function(pos_) {
        # Assigning to obj named <pos_>_df so that the call identifies the position
        assign(paste0(pos_, '_df'), get_df('pos', .pos = pos_))
        lapply(
            pos_ys,
            function(y_) {
                f <- paste(y_, ' ~ clade')
                # Whether to include log_mass covariate
                imc <- {include_mass %>% filter(pos == pos_, y == y_)}$include
                if (imc[1]) f <- paste(f, '+ log_mass')
                arg_list <- list(
                    as.formula(f),
                    data = as.name(paste0(pos_, "_df")),
                    phy = as.name("tr"), model = 'lambda')
                # Some models don't find the peak likelihood unless specifying a
                # starting value of 0.1.
                LL_nostart <- suppressWarnings(do.call("phylolm", arg_list))$logLik
                arg_list_wstart <- c(arg_list, starting.value = quote(list(lambda = 0.1)))
                LL_wstart <- suppressWarnings(do.call("phylolm", arg_list_wstart))$logLik
                if (LL_wstart > LL_nostart) arg_list <- arg_list_wstart
                arg_list <- c(arg_list, boot = 2000)
                # Now create the final phylolm object
                suppressWarnings(do.call("phylolm", arg_list))
            })
    })
names(pos_fits) <- seg_types
for (i in 1:length(pos_fits)) names(pos_fits[[i]]) <- pos_ys; rm(i)
```

The model for `crypt_width ~ clade` in the proximal segment has a higher
log likelihood at a very high phylogenetic signal (`logLik = 68.14` at
`lambda = 0.993`) than at very low signal (`logLik = 66.88` at `lambda
= 1e-7`). However, this model is sensitive to starting values, which
suggests multiple peaks in the likelihood profile. Moreover, the models
for the other segments show very low phylogenetic signal (`1e-7` for
both), and this model re-run with the Ornstein-Uhlenbeck model for
phylogenetic error (“OU”; `OUfixedRoot` in `phylolm`) has a higher log
likelihood and shows a much lower phylogenetic signal (`logLik = 68.54`
and `alpha = 0.0156`). Thus the model likely had convergence issues
using Pagel’s lambda, so I’m replacing the original model with one using
the OU error model below. I’m saving the original one to report it, too.

``` r
pos_fits$prox$crypt_width_pagel <- pos_fits$prox$crypt_width
prox_df <- get_df('pos', .pos = 'prox')
set.seed(1340481016)
pos_fits$prox$crypt_width <- update(pos_fits$prox$crypt_width, 
                                       model = "OUfixedRoot")
```

Saving output:

``` r
readr::write_rds(pos_fits, 'output/models_pos.rds')
```

# `cor_phylo`

From the original manuscript:

> … we used reduced major axis regression (model II regression)… because
> both variables \[X and Y\] were subject to error

Instead of an RMA regression, I’ll be using `phyr::cor_phylo`, which is
similar to `ape::corphylo` but faster and can conduct parametric
bootstrapping. P-values are calculated using bootstrap replicates.

I used the same P-values to determine that I do not need to use
`log_mass` as an independent variable for any of these fits (see
[`docs/03-include_mass`](docs/03-include_mass.md) for more info).

All variables under this section are log-transformed.

## `Clearance` and `SEF`

Clearance = “paracellular probe L-arabinose clearance”

``` r
# Making data frame with both means and SEs:
clear_sef_df <- bind_cols(clear_df %>% select(species, log_sef, log_clear),
                          clear_se_df %>% select(log_sef, log_clear) %>% 
                              rename(log_sef_se = log_sef, log_clear_se = log_clear))

set.seed(1844365955)
clear_sef <- cor_phylo(traits = list(log_sef, log_clear),
                       meas_errors = list(log_sef_se, log_clear_se),
                       species = species, phy = clear_tr, data = clear_sef_df,
                       method = "nelder-mead-nlopt",
                       boot = 2000, max_iter = 1e6)
```

Saving `cor_phylo` object:

``` r
readr::write_rds(clear_sef, 'output/models_corphylo.rds')
```

# Assembling all output into one table

I ran `summ_df` on all models (both `phylolm` and `cor_phylo`) above.
This function summarizes both of these object classes into a single data
frame. See [`R/model_summaries.R`](R/model_summaries.R) for more info.

``` r
mod_summaries <- bind_rows(
    list(
        summ_df(diet_fit),
        summ_df(absorp_fit),
        bind_rows(lapply(spp_fits, summ_df)),
        bind_rows(
            lapply(names(pos_fits), function(p) {
                bind_rows(lapply(pos_fits[[p]], summ_df, .pos = p))
            })),
        summ_df(clear_sef)
        ))
```

I lastly write this summary to a csv file.

``` r
write_csv(mod_summaries, 'output/models_summaries.csv')
```

# Session info

This outlines the package versions I used for these
    analyses.

    ## Session info -------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.4.4 (2018-03-15)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2018-06-10

    ## Packages -----------------------------------------------------------------

    ##  package      * version date       source        
    ##  ape          * 5.1     2018-04-04 CRAN (R 3.4.4)
    ##  assertthat     0.2.0   2017-04-11 CRAN (R 3.4.0)
    ##  backports      1.1.2   2017-12-13 CRAN (R 3.4.3)
    ##  base         * 3.4.4   2018-03-15 local         
    ##  bindr          0.1.1   2018-03-13 CRAN (R 3.4.4)
    ##  bindrcpp     * 0.2.2   2018-03-29 CRAN (R 3.4.4)
    ##  cli            1.0.0   2017-11-05 cran (@1.0.0) 
    ##  codetools      0.2-15  2016-10-05 CRAN (R 3.4.4)
    ##  compiler       3.4.4   2018-03-15 local         
    ##  crayon         1.3.4   2017-09-16 CRAN (R 3.4.1)
    ##  datasets     * 3.4.4   2018-03-15 local         
    ##  devtools       1.13.5  2018-02-18 CRAN (R 3.4.3)
    ##  digest         0.6.15  2018-01-28 CRAN (R 3.4.3)
    ##  dplyr        * 0.7.5   2018-05-19 CRAN (R 3.4.4)
    ##  evaluate       0.10.1  2017-06-24 CRAN (R 3.4.1)
    ##  future         1.8.1   2018-05-03 CRAN (R 3.4.4)
    ##  future.apply   0.2.0   2018-05-01 CRAN (R 3.4.4)
    ##  globals        0.11.0  2018-01-10 CRAN (R 3.4.3)
    ##  glue           1.2.0   2017-10-29 CRAN (R 3.4.2)
    ##  graphics     * 3.4.4   2018-03-15 local         
    ##  grDevices    * 3.4.4   2018-03-15 local         
    ##  grid           3.4.4   2018-03-15 local         
    ##  hms            0.4.2   2018-03-10 CRAN (R 3.4.4)
    ##  htmltools      0.3.6   2017-04-28 cran (@0.3.6) 
    ##  knitr          1.20    2018-02-20 CRAN (R 3.4.3)
    ##  lattice        0.20-35 2017-03-25 CRAN (R 3.4.4)
    ##  listenv        0.7.0   2018-01-21 CRAN (R 3.4.3)
    ##  magrittr       1.5     2014-11-22 CRAN (R 3.4.0)
    ##  Matrix         1.2-14  2018-04-09 CRAN (R 3.4.4)
    ##  memoise        1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods      * 3.4.4   2018-03-15 local         
    ##  nlme           3.1-137 2018-04-07 CRAN (R 3.4.4)
    ##  nloptr         1.0.4   2014-08-04 CRAN (R 3.4.0)
    ##  parallel       3.4.4   2018-03-15 local         
    ##  phylolm      * 2.6     2018-05-31 CRAN (R 3.4.4)
    ##  phyr         * 0.1.5   2018-06-09 local         
    ##  pillar         1.2.3   2018-05-25 CRAN (R 3.4.4)
    ##  pkgconfig      2.0.1   2017-03-21 CRAN (R 3.4.0)
    ##  purrr        * 0.2.5   2018-05-29 CRAN (R 3.4.4)
    ##  R6             2.2.2   2017-06-17 CRAN (R 3.4.0)
    ##  Rcpp           0.12.17 2018-05-18 CRAN (R 3.4.4)
    ##  readr        * 1.1.1   2017-05-16 CRAN (R 3.4.0)
    ##  rlang          0.2.1   2018-05-30 CRAN (R 3.4.4)
    ##  rmarkdown      1.9     2018-03-01 CRAN (R 3.4.3)
    ##  rprojroot      1.3-2   2018-01-03 CRAN (R 3.4.3)
    ##  stats        * 3.4.4   2018-03-15 local         
    ##  stringi        1.2.2   2018-05-02 CRAN (R 3.4.4)
    ##  stringr        1.3.1   2018-05-10 CRAN (R 3.4.4)
    ##  tibble         1.4.2   2018-01-22 CRAN (R 3.4.3)
    ##  tidyr        * 0.8.1   2018-05-18 CRAN (R 3.4.4)
    ##  tidyselect     0.2.4   2018-02-26 CRAN (R 3.4.3)
    ##  tools          3.4.4   2018-03-15 local         
    ##  utf8           1.1.4   2018-05-24 CRAN (R 3.4.4)
    ##  utils        * 3.4.4   2018-03-15 local         
    ##  withr          2.1.2   2018-03-15 CRAN (R 3.4.4)
    ##  yaml           2.1.19  2018-05-01 CRAN (R 3.4.4)
