Test whether to include log\_mass in phylogenetic linear regression
================
Lucas Nell
05 Dec 2017

-   [Necessary data:](#necessary-data)
-   [`SEF` on `Diet`](#sef-on-diet)
-   [`Absorption` on `Clade`](#absorption-on-clade)
-   [`Morphometrics` on `Clade`](#morphometrics-on-clade)
-   [`Morphometrics` on `Clade`, separately by segment](#morphometrics-on-clade-separately-by-segment)
-   [`Clearance` and `SEF`](#clearance-and-sef)
-   [`Clearance` and `log_enterocyte_density`](#clearance-and-log_enterocyte_density)
-   [`Absorption` and `log_total_enterocytes`](#absorption-and-log_total_enterocytes)
-   [Session info](#session-info)

This file determines whether to use `log_mass` in `phylolm` regressions or not. These analyses take about half an hour to run.

All p-values below are for whether the coefficient for log(mass) are not equal to zero.

Necessary data:
===============

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

`SEF` on `Diet`
===============

`phylolm` call and output:

``` r
set.seed(29851644)
diet_fit <- phylolm(log_sef ~ diet + log_mass, data = spp_df, phy = tr,
                    model = 'lambda', boot = 2000)
```

P-value:

    ## P = 0.501

`Absorption` on `Clade`
=======================

`phylolm` call and output:

``` r
set.seed(1092141389)
absorp_fit <- suppressWarnings(  # gives warning about lambda being very low
    phylolm(log_absorp ~ clade + log_mass, data = absorp_df, phy = absorp_tr, 
            model = 'lambda', boot = 2000)
)
```

P-value:

    ## P = 0

`Morphometrics` on `Clade`
==========================

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

`Morphometrics` on `Clade`, separately by segment
=================================================

`phylolm` calls:

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

`Clearance` and `SEF`
=====================

> P-values in all following sections are for whether the coefficient for the `U` matrix is not zero. There are two p-values because I'm including the `U` matrix separately for the first and second `X` matrix parameters.

``` r
# Function to retrieve the U coefficient(s) from a corphylo object
get_U <- function(cp_obj) {
    rn <- rownames(cp_obj$B)[grepl('\\.1', rownames(cp_obj$B))]
    uc <- matrix(as.numeric(cp_obj$B[rn,]), nrow = 1)
    colnames(uc) <- rn
    return(uc)
}

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

# corphylo_cpp run with bootstrapping (takes ~1 min)
set.seed(1844365955)
clear_sef <- corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, U = Umat, 
                          boot = 2000, n_cores = 4, boot_out = get_U)
clear_sef2 <- corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, U = Umat2, 
                          boot = 2000, n_cores = 4, boot_out = get_U)
```

P-values:

    ## P = 0.997

    ## P = 0.308

`Clearance` and `log_enterocyte_density`
========================================

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
```

P-values:

    ## P = 0.947

    ## P = 0.228

`Absorption` and `log_total_enterocytes`
========================================

``` r
Xmat <- cp_mat(absorp_df, c('log_absorp', 'log_total_enterocytes'))
MEmat <- cp_mat(absorp_se_df, c('log_absorp', 'log_total_enterocytes'))

Umat <- list(cbind(absorp_df$log_mass[!is.na(absorp_df$log_mass)]), NULL)
rownames(Umat[[1]]) <- rownames(Xmat)
Umat2 <- list(NULL, cbind(absorp_df$log_mass[!is.na(absorp_df$log_mass)]))
rownames(Umat2[[2]]) <- rownames(Xmat)

# Fit and bootstrap
set.seed(2016097648)
absorp_te <- corphylo_cpp(X = Xmat, phy = absorp_tr, SeM = MEmat, U = Umat, 
                          boot = 2000, n_cores = 4, boot_out = get_U)
# # This one gives numerical issues: non positive definite correlation matrix
# absorp_te2 <- corphylo_cpp(X = Xmat, phy = absorp_tr, SeM = MEmat, U = Umat2, 
#                           boot = 2000, n_cores = 4, boot_out = get_U)
```

P-value:

    ## P = 1

``` r
save(clear_sef, clear_sef2, clear_ed, clear_ed2, absorp_te, 
     file = 'output/inc_mass_corphylo.RData')
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
    ##  date     2017-12-05

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
