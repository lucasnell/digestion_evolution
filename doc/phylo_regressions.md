Phylogenetic linear regression
================
Lucas Nell
04 Dec 2017

-   [`source` the `R` directory](#source-the-r-directory)
-   [`SEF` on `Diet`](#sef-on-diet)
-   [`Absorption` on `Taxon`](#absorption-on-taxon)
-   [`Morphometrics` on `Taxon`](#morphometrics-on-taxon)
-   [`Morphometrics` on `Taxon`, separately by segment](#morphometrics-on-taxon-separately-by-segment)
-   [`Clearance` on `SEF`](#clearance-on-sef)
-   [`Clearance` on `log_enterocyte_density`](#clearance-on-log_enterocyte_density)
-   [`Absorption` on `log_total_enterocytes`](#absorption-on-log_total_enterocytes)
-   [Assembling all output into one object](#assembling-all-output-into-one-object)
-   [Session info](#session-info)

Loading packages:

``` r
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(phylolm)
    library(ape)
})
devtools::load_all('corphyloCpp')
```

    ## Loading corphyloCpp

`source` the `R` directory
==========================

The `R` directory provides functions to summarize `phylolm` objects, run a version of `ape::corphylo` with confidence interval output, and retrieve morphometric, clearance, and absorption data. See `tidy_csvs.md` for more info.

``` r
invisible(sapply(list.files('R', '*.R', full.names = TRUE), source))
```

The `get_tr` function in `R/get_data.R` reads the main phylogenetic tree, cleans species names, and removes unnecessary species from it for a given analysis set.

The function `ci` in `R/model_summaries.R` gets 95% CIs from a bootstrapped `phylolm` model object.

The function `ci_df` creates a tibble with 95% CIs for all parameters in a single model.

`SEF` on `Diet`
===============

Necessary data:

> The `spp_df` data frame is used for both this analysis and `Morphometrics` on `Taxon`. The `tr` tree is used for this analysis, `Morphometrics` on `Taxon`, and `Morphometrics` on `Taxon`, separately by segment.

``` r
spp_df <- get_df('spp')
tr <- get_tr('spp')
```

`phylolm` call and output:

``` r
set.seed(581120)
diet_fit <- phylolm(sef ~ diet, data = spp_df, phy = tr,
                    model = 'lambda', boot = 2000)
```

I'm saving output for this fit because I'll be using that for summarizing.

``` r
readr::write_rds(diet_fit, 'output/models_diet.rds')
```

Summary:

    ## dietOmnivorous: -2.545 (P = 0.123)
    ## dietProtein: -1.341 (P = 0.379)

`Absorption` on `Taxon`
=======================

"Absorption" here means `Fractional absorption / (total intestinal surface)`, where `total intestinal surface = NSA * SEF`

Necessary data:

``` r
absorp_df <- get_df('absorp')
absorp_tr <- get_tr('absorp')
```

`phylolm` call and output:

``` r
set.seed(454094511)
absorp_fit <- suppressWarnings(  # gives warning about lambda being very low
    phylolm(absorp ~ taxon + log_mass, data = absorp_df, phy = absorp_tr, 
            model = 'lambda', boot = 2000)
)
```

I'm saving output for this fit because I'll be using that for plotting.

``` r
readr::write_rds(absorp_fit, 'output/models_absorp.rds')
```

Summary:

    ## taxonBat: 0.004354 (P = 0)

`Morphometrics` on `Taxon`
==========================

List of `Morphometrics`:

-   Intestinal length
-   NSA
-   Villus surface area
-   Total number of enterocytes (log-transformed)
    -   Calculated as such: `log(NSA * enterocyte_density)`

> log(body mass) as covariate for all

These are the column names for the above parameters:

``` r
spp_ys <- c("intestinal_length", "nsa", "vill_surface_area", "log_total_enterocytes")
```

Necessary data: `spp_df` and `tr` are already created from fitting `sef ~ diet`.

`phylolm` call:

> The actual analyses (takes ~7.5 min, which is why I saved the output):

``` r
set.seed(88754829)
spp_fits <- lapply(
    spp_ys,
    function(y) {
        f <- paste(y, '~ taxon + log_mass')
        suppressWarnings(
            do.call("phylolm", list(as.formula(f), data = as.name("spp_df"),
                                    phy = as.name("tr"), model = 'lambda',
                                    boot = 2000))
        )
    })
names(spp_fits) <- spp_ys
readr::write_rds(spp_fits, 'output/models_spp.rds')
```

Loading the output and summarizing:

``` r
spp_fits <- readr::read_rds('output/models_spp.rds')
p <- 'taxonBat'
for (nn in names(spp_fits)) {
    cat(sprintf('Y = %s\n', nn))
    cat(sprintf('  %s: %.4g (P = %.4g)\n', p,
                    coef(spp_fits[[nn]])[p], 
                    pval(spp_fits[[nn]], p)))
}
```

    ## Y = intestinal_length
    ##   taxonBat: -11.67 (P = 0.024)
    ## Y = nsa
    ##   taxonBat: -7.001 (P = 0.007)
    ## Y = vill_surface_area
    ##   taxonBat: 33.71 (P = 0.339)
    ## Y = log_total_enterocytes
    ##   taxonBat: 0.2853 (P = 0.318)

`Morphometrics` on `Taxon`, separately by segment
=================================================

(Segment = proximal, medial, or distal)

List of `Y`s:

-   Intestinal diameter (log-transformed, body mass as covariate)
-   Villus height (body mass as covariate)
-   Villus width
-   Crypt width
-   Surface enlargement factor (SEF)
-   Enterocyte diameter
-   Enterocytes per cm^2 NSA (log-transformed)

Below are the column names for these parameters and all the segment types.

``` r
pos_ys <- c('log_intestinal_diameter', 'villus_height', 'villus_width', 
            'crypt_width', 'sef', 'enterocyte_diameter', 'log_enterocyte_density')
seg_types <- c('prox', 'med', 'dist')
```

Below is a data frame including whether or not to include `log_mass` as a covariate. This determination was based on whether `log_mass` had a significant effect when it was included in the model, where p-values were based on parametric bootstrapping (see `docs/include_mass.md`).

    ## # A tibble: 21 x 3
    ##      pos                       y include
    ##    <chr>                   <chr>   <lgl>
    ##  1  dist log_intestinal_diameter    TRUE
    ##  2  dist           villus_height    TRUE
    ##  3  dist            villus_width   FALSE
    ##  4  dist             crypt_width   FALSE
    ##  5  dist                     sef   FALSE
    ##  6  dist     enterocyte_diameter   FALSE
    ##  7  dist  log_enterocyte_density   FALSE
    ##  8   med log_intestinal_diameter    TRUE
    ##  9   med           villus_height    TRUE
    ## 10   med            villus_width   FALSE
    ## # ... with 11 more rows

`phylolm` call:

The actual analyses (takes ~21.7 min, which is why I saved the output):

``` r
set.seed(25413535)
pos_fits <- lapply(
    seg_types,
    function(pos) {
        # Assigning to obj named <pos>_df so that the call identifies the position
        assign(paste0(pos, '_df'), get_df('pos', .pos = pos))
        lapply(
            pos_ys,
            function(y) {
                f <- paste(y, ' ~ taxon')
                # Whether to include log_mass covariate
                imc <- {include_mass %>% filter(pos == pos, y == y)}$include
                if (imc[1]) f <- paste(f, '+ log_mass')
                arg_list <- list(
                    as.formula(f),
                    data = as.name(paste0(pos, "_df")),
                    phy = as.name("tr"), model = 'lambda')
                # Some models don't find the peak likelihood unless specifying a
                # starting value of 0.1.
                LL_nostart <- suppressWarnings(do.call("phylolm", arg_list))$logLik
                LL_wstart <- suppressWarnings(do.call(
                    "phylolm", c(arg_list, starting.value = 0.1)))$logLik
                if (LL_wstart > LL_nostart) {
                    arg_list <- c(arg_list, starting.value = 0.1)
                }
                arg_list <- c(arg_list, boot = 2000)
                # Now create the final phylolm object
                suppressWarnings(do.call("phylolm", arg_list))
            })
    })
names(pos_fits) <- seg_types
for (i in 1:length(pos_fits)) names(pos_fits[[i]]) <- pos_ys; rm(i)
readr::write_rds(pos_fits, 'output/models_pos.rds')
```

Loading the output and summarizing:

``` r
pos_fits <- readr::read_rds('output/models_pos.rds')
cbind(sapply(pos_fits$dist, pval))
```

    ##                          [,1]
    ## log_intestinal_diameter 0.912
    ## villus_height           0.000
    ## villus_width            0.554
    ## crypt_width             0.006
    ## sef                     0.000
    ## enterocyte_diameter     0.577
    ## log_enterocyte_density  0.001

``` r
cbind(sapply(pos_fits$med, pval))
```

    ##                          [,1]
    ## log_intestinal_diameter 0.959
    ## villus_height           0.001
    ## villus_width            0.085
    ## crypt_width             0.092
    ## sef                     0.000
    ## enterocyte_diameter     0.050
    ## log_enterocyte_density  0.397

``` r
cbind(sapply(pos_fits$prox, pval))
```

    ##                          [,1]
    ## log_intestinal_diameter 0.527
    ## villus_height           0.140
    ## villus_width            0.012
    ## crypt_width             0.310
    ## sef                     0.000
    ## enterocyte_diameter     0.336
    ## log_enterocyte_density  0.003

`Clearance` on `SEF`
====================

Clearance = "paracellular probe L-arabinose clearance"

Both are log-transformed.

From the original manuscript: &gt; ... we used reduced major axis regression (model II regression)... because both &gt; variables \[X and Y\] were subject to error

Instead of an RMA regression, I'll be using a modified version of `ape::corphylo` that can conduct parametric bootstrapping. P-values are calculated using these bootstrapping replicates.

``` r
clear_df <- get_df('clear')
clear_se_df <- get_df('clear', .stat = 'se')  # <-- contains standard errors
clear_tr <- get_tr('clear')

Xmat <- cp_mat(clear_df, c('log_sef', 'log_clear'))
MEmat <- cp_mat(clear_se_df, c('log_sef', 'log_clear'))

# Using p-values, this doesn't need Umat for either parameter

# corphylo_cpp run with bootstrapping (takes ~1 min)
set.seed(1844365955)
clear_sef <- corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, boot = 2000, n_cores = 4)
clear_sef
```

    ## Call to corphylo
    ## 
    ## logLik    AIC    BIC 
    ## -18.58  51.16  44.52 
    ## 
    ## correlation matrix:
    ##        1      2
    ## 1 1.0000 0.5765
    ## 2 0.5765 1.0000
    ## 
    ## from OU process:
    ##        d
    ## 1 1.1383
    ## 2 0.5887
    ## 
    ## coefficients:
    ##        Value Std.Error  Zscore  Pvalue    
    ## B1.0 2.41029   0.23771 10.1396 < 2e-16 ***
    ## B2.0 0.97926   0.48545  2.0172 0.04367 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Bootstrapped 95% CI:
    ##   r12       0.56 [    -0.14     0.947]
    ##   d1        1.01 [ 9.22e-12      3.66]
    ##   d2       0.635 [ 2.21e-14      2.83]

``` r
# P-value for correlation != 0
pval(clear_sef)[1]
```

    ## [1] 0.086

`Clearance` on `log_enterocyte_density`
=======================================

> One species (*Rattus norvegicus*) doesn't have `log_enterocyte_density` data, which is why I'm removing that row below.

``` r
Xmat <- cp_mat(clear_df, c('log_enterocyte_density', 'log_clear'))
Xmat <- Xmat[!is.na(rowSums(Xmat)),]

MEmat <- cp_mat(clear_se_df, c('log_enterocyte_density', 'log_clear'))
MEmat <- MEmat[!is.na(rowSums(MEmat)),]
# Using p-values, this doesn't need Umat for either parameter

clear_ed_tr <- ape::drop.tip(
    clear_tr, 
    tip = clear_tr$tip.label[!clear_tr$tip.label %in% rownames(Xmat)]
)

# Fit and bootstrap r (takes ~1 min)
set.seed(1442148819)
clear_ed <- corphylo_cpp(Xmat, phy = clear_ed_tr, SeM = MEmat, boot = 2000, n_cores = 4)

# P-value for correlation != 0
pval(clear_ed)[1]
```

    ## [1] 0.654

`Absorption` on `log_total_enterocytes`
=======================================

``` r
absorp_df <- get_df('absorp')
absorp_se_df <- get_df('absorp', .stat = 'se')  # <-- contains standard errors
absorp_tr <- get_tr('absorp')

Xmat <- cp_mat(absorp_df, c('absorp', 'log_total_enterocytes'))
MEmat <- cp_mat(absorp_se_df, c('absorp', 'log_total_enterocytes'))
# Using p-values, this doesn't need Umat for either parameter

# Fit and bootstrap
set.seed(2016097648)
absorp_te <- corphylo_cpp(Xmat, phy = absorp_tr, SeM = MEmat, boot = 2000, n_cores = 4)

# P-value for correlation != 0
pval(absorp_te)[1]
```

    ## [1] 0

Assembling all output into one object
=====================================

I ran `summ_df` on all models above. This function summarizes `phylolm` and `corphylo` objects.

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
        summ_df(clear_sef, .corr_pars = c('log_sef', 'log_clear')),
        summ_df(clear_ed, .corr_pars = c('log_enterocyte_density', 'log_clear')),
        summ_df(absorp_te, .corr_pars = c('absorp', 'log_total_enterocytes'))
        ))
```

I lastly write this summary to a csv file.

``` r
write_csv(mod_summaries, 'output/models_summaries.csv')
```

Session info
============

This outlines the package versions I used for these analyses.

    ## Session info -------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.4.2 (2017-09-28)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2017-12-04

    ## Packages -----------------------------------------------------------------

    ##  package     * version date       source        
    ##  ape         * 5.0     2017-10-30 CRAN (R 3.4.2)
    ##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.0)
    ##  backports     1.1.1   2017-09-25 CRAN (R 3.4.2)
    ##  base        * 3.4.2   2017-10-04 local         
    ##  bindr         0.1     2016-11-13 CRAN (R 3.4.0)
    ##  bindrcpp    * 0.2     2017-06-17 CRAN (R 3.4.0)
    ##  commonmark    1.4     2017-09-01 CRAN (R 3.4.1)
    ##  compiler      3.4.2   2017-10-04 local         
    ##  corphyloCpp * 1.0     <NA>       local         
    ##  datasets    * 3.4.2   2017-10-04 local         
    ##  devtools      1.13.3  2017-08-02 CRAN (R 3.4.1)
    ##  digest        0.6.12  2017-01-27 CRAN (R 3.4.0)
    ##  dplyr       * 0.7.4   2017-09-28 CRAN (R 3.4.2)
    ##  evaluate      0.10.1  2017-06-24 CRAN (R 3.4.1)
    ##  glue          1.2.0   2017-10-29 CRAN (R 3.4.2)
    ##  graphics    * 3.4.2   2017-10-04 local         
    ##  grDevices   * 3.4.2   2017-10-04 local         
    ##  grid          3.4.2   2017-10-04 local         
    ##  hms           0.3     2016-11-22 CRAN (R 3.4.0)
    ##  htmltools     0.3.6   2017-04-28 cran (@0.3.6) 
    ##  knitr         1.17    2017-08-10 CRAN (R 3.4.1)
    ##  lattice       0.20-35 2017-03-25 CRAN (R 3.4.2)
    ##  magrittr      1.5     2014-11-22 CRAN (R 3.4.0)
    ##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods     * 3.4.2   2017-10-04 local         
    ##  nlme          3.1-131 2017-02-06 CRAN (R 3.4.2)
    ##  parallel      3.4.2   2017-10-04 local         
    ##  phylolm     * 2.5     2016-10-17 CRAN (R 3.4.0)
    ##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.0)
    ##  purrr       * 0.2.4   2017-10-18 CRAN (R 3.4.2)
    ##  R6            2.2.2   2017-06-17 CRAN (R 3.4.0)
    ##  Rcpp          0.12.13 2017-09-28 CRAN (R 3.4.2)
    ##  readr       * 1.1.1   2017-05-16 CRAN (R 3.4.0)
    ##  rlang         0.1.4   2017-11-05 CRAN (R 3.4.2)
    ##  rmarkdown     1.6     2017-06-15 CRAN (R 3.4.0)
    ##  roxygen2      6.0.1   2017-02-06 CRAN (R 3.4.0)
    ##  rprojroot     1.2     2017-01-16 cran (@1.2)   
    ##  stats       * 3.4.2   2017-10-04 local         
    ##  stringi       1.1.5   2017-04-07 CRAN (R 3.4.0)
    ##  stringr       1.2.0   2017-02-18 CRAN (R 3.4.0)
    ##  tibble        1.3.4   2017-08-22 CRAN (R 3.4.1)
    ##  tidyr       * 0.7.2   2017-10-16 CRAN (R 3.4.2)
    ##  tools         3.4.2   2017-10-04 local         
    ##  utils       * 3.4.2   2017-10-04 local         
    ##  withr         2.1.0   2017-11-01 CRAN (R 3.4.2)
    ##  xml2          1.1.1   2017-01-24 CRAN (R 3.4.0)
    ##  yaml          2.1.14  2016-11-12 cran (@2.1.14)
