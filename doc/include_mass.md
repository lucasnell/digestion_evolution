Include log\_mass in phylogenetic linear regression
================
Lucas Nell
01 Dec 2017

-   [Loading packages:](#loading-packages)
-   [`SEF` on `Diet`](#sef-on-diet)
-   [`Absorption` on `Taxon`](#absorption-on-taxon)
-   [`Morphometrics` on `Taxon`](#morphometrics-on-taxon)
-   [`Morphometrics` on `Taxon`, separately by segment](#morphometrics-on-taxon-separately-by-segment)

This file determines whether to use `log_mass` in `phylolm` regressions or not. These analyses take about half an hour to run.

#### Loading packages:

``` r
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(phylolm)
    library(ape)
})
```

``` r
invisible(sapply(list.files('R', '*.R', full.names = TRUE), source))
```

`SEF` on `Diet`
===============

``` r
spp_df <- get_df('spp')
tr <- get_tr('spp')
```

`phylolm` call and output:

``` r
set.seed(29851644)
diet_fit <- phylolm(sef ~ diet + log_mass, data = spp_df, phy = tr,
                    model = 'lambda', boot = 2000)
pval(diet_fit, 'log_mass')
```

    ## [1] 0.369

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
set.seed(1092141389)
absorp_fit <- suppressWarnings(  # gives warning about lambda being very low
    phylolm(absorp ~ taxon + log_mass, data = absorp_df, phy = absorp_tr, 
            model = 'lambda', boot = 2000)
)
pval(absorp_fit, 'log_mass')
```

    ## [1] 0

`Morphometrics` on `Taxon`
==========================

``` r
spp_ys <- c("intestinal_length", "nsa", "vill_surface_area", "log_total_enterocytes")
```

``` r
set.seed(357885189)
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

cbind(sapply(spp_fits, pval, 'log_mass'))
```

    ##                       [,1]
    ## intestinal_length     0.04
    ## nsa                   0.00
    ## vill_surface_area     0.00
    ## log_total_enterocytes 0.00

`Morphometrics` on `Taxon`, separately by segment
=================================================

``` r
pos_ys <- c('log_intestinal_diameter', 'villus_height', 'villus_width', 
            'crypt_width', 'sef', 'enterocyte_diameter', 'log_enterocyte_density')
seg_types <- c('prox', 'med', 'dist')
```

`phylolm` call:

``` r
set.seed(632929430)
pos_fits <- lapply(
    seg_types,
    function(pos) {
        # Assigning to obj named <pos>_df so that the call identifies the position
        assign(paste0(pos, '_df'), get_df('pos', .pos = pos))
        lapply(
            pos_ys,
            function(y) {
                f <- paste(y, ' ~ taxon + log_mass')
                arg_list <- list(
                    as.formula(f),
                    data = as.name(paste0(pos, "_df")),
                    phy = as.name("tr"), model = 'lambda',
                    boot = 2000)
                # These models don't find the peak likelihood unless specifying a
                # starting value of 0.1.
                if ((y == "log_enterocyte_density" & pos == "med") |
                    (y == "crypt_width" & pos == "prox")) {
                    arg_list <- c(arg_list, starting.value = 0.1)
                }
                # Now call phylolm
                suppressWarnings(do.call("phylolm", arg_list))
            })
    })
names(pos_fits) <- seg_types
for (i in 1:length(pos_fits)) names(pos_fits[[i]]) <- pos_ys; rm(i)
cbind(sapply(pos_fits$dist, pval, parameter = 'log_mass'))
```

    ##                          [,1]
    ## log_intestinal_diameter 0.000
    ## villus_height           0.012
    ## villus_width            0.193
    ## crypt_width             0.500
    ## sef                     0.283
    ## enterocyte_diameter     0.782
    ## log_enterocyte_density  0.462

``` r
cbind(sapply(pos_fits$med, pval, parameter = 'log_mass'))
```

    ##                          [,1]
    ## log_intestinal_diameter 0.006
    ## villus_height           0.013
    ## villus_width            0.062
    ## crypt_width             0.379
    ## sef                     0.349
    ## enterocyte_diameter     0.621
    ## log_enterocyte_density  0.477

``` r
cbind(sapply(pos_fits$prox, pval, parameter = 'log_mass'))
```

    ##                          [,1]
    ## log_intestinal_diameter 0.000
    ## villus_height           0.001
    ## villus_width            0.009
    ## crypt_width             0.807
    ## sef                     0.024
    ## enterocyte_diameter     0.823
    ## log_enterocyte_density  0.203
