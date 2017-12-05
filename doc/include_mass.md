Include log\_mass in phylogenetic linear regression
================
Lucas Nell
04 Dec 2017

-   [Loading packages:](#loading-packages)
-   [`SEF` on `Diet`](#sef-on-diet)
-   [`Absorption` on `Taxon`](#absorption-on-taxon)
-   [`Morphometrics` on `Taxon`](#morphometrics-on-taxon)
-   [`Morphometrics` on `Taxon`, separately by segment](#morphometrics-on-taxon-separately-by-segment)
-   [`Clearance` on `SEF`](#clearance-on-sef)
-   [`Clearance` on `log_enterocyte_density`](#clearance-on-log_enterocyte_density)
-   [`Absorption` on `log_total_enterocytes`](#absorption-on-log_total_enterocytes)

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
devtools::load_all('corphyloCpp')
```

    ## Loading corphyloCpp

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
                    phy = as.name("tr"), model = 'lambda')
                # Some models don't find the peak likelihood unless specifying a
                # starting value of 0.1.
                LL_nostart <- suppressWarnings(do.call("phylolm", arg_list))$logLik
                LL_wstart <- suppressWarnings(do.call(
                    "phylolm", c(arg_list, starting.value = 0.1)))$logLik
                if (LL_wstart > LL_nostart) arg_list <- c(arg_list, starting.value = 0.1)
                # Adding number of bootstrap replicates
                arg_list <- c(arg_list, boot = 2000)
                # Now call phylolm
                suppressWarnings(do.call("phylolm", arg_list))
            })
    })
names(pos_fits) <- seg_types
for (i in 1:length(pos_fits)) names(pos_fits[[i]]) <- pos_ys; rm(i)
cbind(sapply(pos_fits$dist, pval, .parameters = 'log_mass'))
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
cbind(sapply(pos_fits$med, pval, .parameters = 'log_mass'))
```

    ##                          [,1]
    ## log_intestinal_diameter 0.006
    ## villus_height           0.013
    ## villus_width            0.062
    ## crypt_width             0.379
    ## sef                     0.349
    ## enterocyte_diameter     0.621
    ## log_enterocyte_density  0.208

``` r
cbind(sapply(pos_fits$prox, pval, .parameters = 'log_mass'))
```

    ##                          [,1]
    ## log_intestinal_diameter 0.000
    ## villus_height           0.001
    ## villus_width            0.009
    ## crypt_width             0.881
    ## sef                     0.024
    ## enterocyte_diameter     0.823
    ## log_enterocyte_density  0.203

`Clearance` on `SEF`
====================

> P-values in all following chunks are for whether the coefficient for the `U` matrix is not zero. There are two p-values because I'm including the `U` matrix separately for the first and second `X` matrix parameters.

``` r
clear_df <- get_df('clear')
clear_se_df <- get_df('clear', .stat = 'se')  # <-- contains standard errors
clear_tr <- get_tr('clear')

Xmat <- cp_mat(clear_df, c('log_sef', 'log_clear'))
MEmat <- cp_mat(clear_se_df, c('log_sef', 'log_clear'))

# For this comparison, I have to remove one row that doesn't have log_mass
Xmat <- Xmat[!is.na(clear_df$log_mass),]
MEmat <- MEmat[!is.na(clear_df$log_mass),]
Umat <- list( cbind(clear_df$log_mass[!is.na(clear_df$log_mass)]), NULL)
rownames(Umat[[1]]) <- rownames(Xmat)
Umat2 <- list(NULL, cbind(clear_df$log_mass[!is.na(clear_df$log_mass)]))
rownames(Umat2[[2]]) <- rownames(Xmat)

clear_tr <- ape::drop.tip(
    clear_tr,
    tip = clear_tr$tip.label[!clear_tr$tip.label %in% rownames(Xmat)])

# Function to retrieve the U coefficient from a corphylo object
get_U <- function(cp_obj) {
    rn <- rownames(cp_obj$B)[grepl('\\.1', rownames(cp_obj$B))]
    matrix(as.numeric(cp_obj$B[rn,]), nrow = 1)
}

# corphylo_cpp run with bootstrapping (takes ~1 min)
set.seed(1844365955)
clear_sef <- corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, U = Umat, 
                          boot = 2000, n_cores = 4, boot_out = get_U)
clear_sef2 <- corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, U = Umat2, 
                          boot = 2000, n_cores = 4, boot_out = get_U)

# P-value for coefficient != 0
pval(clear_sef); pval(clear_sef2)
```

    ## [1] 0.997

    ## [1] 0.308

`Clearance` on `log_enterocyte_density`
=======================================

``` r
Xmat <- cp_mat(clear_df, c('log_enterocyte_density', 'log_clear'))
Xmat <- Xmat[!is.na(rowSums(Xmat)),]

MEmat <- cp_mat(clear_se_df, c('log_enterocyte_density', 'log_clear'))
MEmat <- MEmat[!is.na(rowSums(MEmat)),]

# Fit and bootstrap r (takes ~1 min)
set.seed(1442148819)
clear_ed <- corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, U = Umat, 
                          boot = 2000, n_cores = 4, boot_out = get_U)
clear_ed2 <- corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, U = Umat2, 
                          boot = 2000, n_cores = 4, boot_out = get_U)

# P-value for coefficient != 0
pval(clear_ed); pval(clear_ed2)
```

    ## [1] 0.947

    ## [1] 0.228

`Absorption` on `log_total_enterocytes`
=======================================

``` r
absorp_df <- get_df('absorp')
absorp_se_df <- get_df('absorp', .stat = 'se')  # <-- contains standard errors
absorp_tr <- get_tr('absorp')


Xmat <- cp_mat(absorp_df, c('absorp', 'log_total_enterocytes'))
MEmat <- cp_mat(absorp_se_df, c('absorp', 'log_total_enterocytes'))

Umat <- list( cbind(absorp_df$log_mass[!is.na(absorp_df$log_mass)]), NULL)
rownames(Umat[[1]]) <- rownames(Xmat)
Umat2 <- list(NULL, cbind(absorp_df$log_mass[!is.na(absorp_df$log_mass)]))
rownames(Umat2[[2]]) <- rownames(Xmat)


# Fit and bootstrap
set.seed(2016097648)
absorp_te <- corphylo_cpp(X = Xmat, phy = absorp_tr, SeM = MEmat, U = Umat, 
                          boot = 2000, n_cores = 4, boot_out = get_U)
absorp_te2 <- corphylo_cpp(X = Xmat, phy = absorp_tr, SeM = MEmat, U = Umat2, 
                          boot = 2000, n_cores = 4, boot_out = get_U)

# P-value for coefficient != 0
pval(absorp_te); pval(absorp_te2)
```

    ## [1] 0.981

    ## [1] 0.42