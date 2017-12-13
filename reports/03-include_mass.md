Test whether to include log\_mass in phylogenetic linear regression
================
Lucas Nell
13 Dec 2017

-   [Necessary data:](#necessary-data)
-   [`SEF` on `Diet`](#sef-on-diet)
-   [`Absorption` on `Clade`](#absorption-on-clade)
-   [`Morphometrics` on `Clade`](#morphometrics-on-clade)
-   [`Morphometrics` on `Clade`, separately by segment](#morphometrics-on-clade-separately-by-segment)
-   [`Clearance` and `SEF`](#clearance-and-sef)
-   [Session info](#session-info)

This file determines whether to use `log_mass` in `phylolm` regressions and `corphylo` correlations. These analyses take about half an hour to run.

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

    ## P = 0.000

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

Because there are so many models, I'm writing whether to include mass to a csv file:

``` r
include_df <- lapply(
    seg_types, 
    function(p) {
        data_frame(pos = p, y = names(pos_fits[[p]]),
                   include = sapply(pos_fits[[p]], pval, 
                                    .parameters = 'log_mass') < 0.05)
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

`Clearance` and `SEF`
=====================

In `corphylo`, you can input a `U` object if one or more of your variables of interest has a covariate that might be having confounding effects. `U` is a list of length `p`, where `p` is the number of variables you're interested in getting correlations between (in my case, `p = 2`). If `U[[i]]` is `NULL`, then variable `i` is considered to not have a covariate, while if `U[[i]]` is a matrix, then each column in that matrix is considered a covariate for variable `i`. Below, I'm trying out whether either clearance or SEF needs body mass as a covariate by including a 1-column matrix of body mass values in the `U` object, first for the position associated with SEF, then for the position associated with clearance.

P-values in this section are for whether the coefficient for the `U` matrix is not zero.

``` r
# Function to retrieve the U coefficient(s) from a corphylo object
get_U <- function(cp_obj) {
    rn <- rownames(cp_obj$B)[grepl('\\.1', rownames(cp_obj$B))]
    uc <- matrix(as.numeric(cp_obj$B[rn,]), nrow = 1)
    colnames(uc) <- rn
    return(uc)
}
# Matrix of mean SEF and clearance values by species
Xmat <- cp_mat(clear_df, c('log_sef', 'log_clear'))
# Matrix of standard error SEF and clearance values by species
MEmat <- cp_mat(clear_se_df, c('log_sef', 'log_clear'))

# For this comparison, I have to remove one row that doesn't have body mass
Xmat <- Xmat[!is.na(clear_df$log_mass),]
MEmat <- MEmat[!is.na(clear_df$log_mass),]

# Now creating two U objects, one for having it as a covariate for SEF, 
# then another for clearance
U_sef <- list( cbind(clear_df$log_mass[!is.na(clear_df$log_mass)]), NULL)
rownames(U_sef[[1]]) <- rownames(Xmat)
U_clear <- list(NULL, cbind(clear_df$log_mass[!is.na(clear_df$log_mass)]))
rownames(U_clear[[2]]) <- rownames(Xmat)

# Phylogenetic tree, removing species with no body mass data
clear_tr <- ape::drop.tip(
    clear_tr,
    tip = clear_tr$tip.label[!clear_tr$tip.label %in% rownames(Xmat)])

# corphylo_cpp run with bootstrapping (takes ~1 min each)
set.seed(1844365955)
clear_sef <- list(
    sef = corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, U = U_sef, 
                       boot = 2000, n_cores = 4, boot_out = get_U),
    clear = corphylo_cpp(X = Xmat, phy = clear_tr, SeM = MEmat, U = U_clear, 
                       boot = 2000, n_cores = 4, boot_out = get_U))
```

P-values for including body mass for SEF and clearance, respectively:

    ## P = 0.997

    ## P = 0.308

Writing this object to an `rds` file for later.

``` r
write_rds(clear_sef, 'output/inc_mass_corphylo.rds')
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
    ##  date     2017-12-13

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
