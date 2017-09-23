Phylogenetic linear regression
================
Lucas Nell
23 Sep 2017

-   [`source` the `R` directory](#source-the-r-directory)
-   [`SEF ~ Diet`](#sef-diet)
-   [`Absorption ~ Taxon`](#absorption-taxon)
-   [`Morphometrics ~ Taxon`](#morphometrics-taxon)
-   [`Morphometrics ~ Taxon`, separately by segment](#morphometrics-taxon-separately-by-segment)
-   [`Clearance ~ SEF`](#clearance-sef)
-   [Assembling all output into one object](#assembling-all-output-into-one-object)
-   [Session info](#session-info)

Loading packages:

``` r
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(phylolm)
    library(ape)
})
```

`source` the `R` directory
==========================

The `R` directory provides functions to summarize `phylolm` objects, run a version of `ape::corphylo` with confidence interval output, and retrieve morphometric, clearance, and absorption data. See `tidy_csvs.md` for more info.

``` r
invisible(sapply(list.files('R', '*.R', full.names = TRUE), source))
```

The `get_tr` function in `R/get_data.R` reads the main phylogenetic tree, cleans species names, and removes unnecessary species from it for a given analysis set.

The function `ci` in `R/model_summaries.R` gets 95% CIs from a bootstrapped `phylolm` model object.

The function `ci_df` creates a tibble with 95% CIs for all parameters in a single model.

`SEF ~ Diet`
============

Necessary data:

``` r
diet_df <- get_df('diet')
diet_tr <- get_tr('diet')
```

`phylolm` call and output:

``` r
set.seed(581120)
diet_fit <- phylolm(sef ~ diet, data = diet_df, phy = diet_tr, 
                    model = 'lambda', boot = 2000)
ci(diet_fit, c('dietOmnivorous', 'dietProtein'))
```

    ##       dietOmnivorous dietProtein
    ## 2.5%       -5.708251   -4.179319
    ## 97.5%       1.298453    2.299283

``` r
summary(diet_fit)
```

    ## 
    ## Call:
    ## phylolm(formula = sef ~ diet, data = diet_df, phy = diet_tr, 
    ##     model = "lambda", boot = 2000)
    ## 
    ##    AIC logLik 
    ##   87.2  -38.6 
    ## 
    ## Raw residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6736 -3.2315 -0.5839  3.4188  6.3726 
    ## 
    ## Mean tip height: 96.46239
    ## Parameter estimate(s) using ML:
    ## lambda : 0.718295
    ## sigma2: 0.118923 
    ## 
    ## Coefficients:
    ##                Estimate   StdErr  t.value lowerbootCI upperbootCI
    ## (Intercept)    13.47067  2.29622  5.86646     9.42748     17.4570
    ## dietOmnivorous -2.30099  1.95193 -1.17883    -5.70825      1.2985
    ## dietProtein    -1.05811  1.78688 -0.59216    -4.17932      2.2993
    ##                  p.value    
    ## (Intercept)    5.533e-05 ***
    ## dietOmnivorous    0.2596    
    ## dietProtein       0.5639    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Note: p-values are conditional on lambda=0.718295.
    ## 
    ## sigma2: 0.118923
    ##       bootstrap mean: 0.08165768 (on raw scale)
    ##                       0.06906805 (on log scale, then back transformed)
    ##       bootstrap 95% CI: (0.02363493,0.2130501)
    ## 
    ## lambda: 0.718295
    ##       bootstrap mean: 0.3968617 (on raw scale)
    ##       bootstrap 95% CI: (1e-07,1)
    ## 
    ## Parametric bootstrap results based on 2000 fitted replicates

`Absorption ~ Taxon`
====================

Necessary data:

``` r
absorp_df <- get_df('absorp')
absorp_tr <- get_tr('absorp')
```

`phylolm` call and output:

``` r
set.seed(454094511)
absorp_fit <- suppressWarnings(  # gives warning about lambda being very low
    phylolm(absorp ~ taxon, data = absorp_df, phy = absorp_tr, 
            model = 'lambda', boot = 2000)
)
ci(absorp_fit)
```

    ##       2.5%      97.5% 
    ## 0.03289512 0.06095654

``` r
summary(absorp_fit)
```

    ## 
    ## Call:
    ## phylolm(formula = absorp ~ taxon, data = absorp_df, phy = absorp_tr, 
    ##     model = "lambda", boot = 2000)
    ## 
    ##    AIC logLik 
    ## -42.73  25.37 
    ## 
    ## Raw residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -0.0186031 -0.0030835 -0.0006498  0.0030048  0.0191412 
    ## 
    ## Mean tip height: 96.46239
    ## Parameter estimate(s) using ML:
    ## lambda : 1e-07
    ## sigma2: 1.06921e-06 
    ## 
    ## Coefficients:
    ##              Estimate    StdErr   t.value lowerbootCI upperbootCI  p.value
    ## (Intercept) 0.0180469 0.0058634 3.0778850   0.0081828      0.0281 0.021722
    ## taxonBat    0.0466086 0.0082921 5.6208371   0.0328951      0.0610 0.001355
    ##               
    ## (Intercept) * 
    ## taxonBat    **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Note: p-values are conditional on lambda=1e-07.
    ## 
    ## sigma2: 1.06921e-06
    ##       bootstrap mean: 8.161444e-07 (on raw scale)
    ##                       6.878408e-07 (on log scale, then back transformed)
    ##       bootstrap 95% CI: (1.722724e-07,1.962124e-06)
    ## 
    ## lambda: 1e-07
    ##       bootstrap mean: 1e-07 (on raw scale)
    ##       bootstrap 95% CI: (1e-07,1e-07)
    ## 
    ## Parametric bootstrap results based on 2000 fitted replicates

I'm saving output for this fit because I'll be using that for plotting.

``` r
readr::write_rds(absorp_fit, 'output/models_absorp.rds')
```

`Morphometrics ~ Taxon`
=======================

List of `Morphometrics`:

-   Intestinal length / body mass^0.4
-   NSA / body mass^0.75
-   Villus surface area / body mass^0.75
-   Total number of enterocytes (log-transformed; log body mass as covariate)
    -   Calculated as such: `NSA * mean(<enterocyte density among segments>)`
-   Fractional absorption / (total intestinal surface / mass^0.75)
    -   total intestinal surface = `NSA * SEF`

These are the column names for the above parameters:

``` r
spp_ys <- c("int_length_mass", "nsa_mass", "vill_area_mass", "log_total_enterocytes")
```

Necessary data:

> The `tr` tree is used for both this analysis and the one separated also by segment.

``` r
spp_df <- get_df('spp')
tr <- get_tr('spp')
```

`phylolm` call:

> The actual analyses (takes ~6.5 min, which is why I saved the output):

``` r
set.seed(88754829)
spp_fits <- lapply(
    spp_ys,
    function(y) {
        f <- paste(y, ' ~ taxon',
                   ifelse(grepl('total_enterocytes', y), '+ log_mass', ''))
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
sapply(spp_fits, ci)
```

    ##       int_length_mass   nsa_mass vill_area_mass log_total_enterocytes
    ## 2.5%       -4.9166135 -0.7478906      -2.867647            -0.2823279
    ## 97.5%      -0.6899693 -0.1411825       4.987878             0.8362962

``` r
ci(spp_fits$log_total_enterocytes, parameter = 'log_mass')
```

    ##      2.5%     97.5% 
    ## 0.3237234 1.0874341

`Morphometrics ~ Taxon`, separately by segment
==============================================

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
                f <- paste(y, ' ~ taxon',
                           ifelse(grepl('intestinal_diameter|villus_height', y),
                                  '+ log_mass', ''))
                arg_list <- list(
                    as.formula(f),
                    data = as.name(paste0(pos, "_df")),
                    phy = as.name("tr"), model = 'lambda',
                    boot = 2000)
                # This model doesn't find the peak likelihood unless specifying a
                # starting value of 0.1.
                if (y == "log_enterocyte_density" & pos == "prox") {
                    arg_list <- c(arg_list, starting.value = 0.1)
                }
                # Now call phylolm
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
sapply(pos_fits$dist, ci)
```

    ##       log_intestinal_diameter villus_height villus_width  crypt_width
    ## 2.5%               -0.1690558     0.1292392  -0.02308230 -0.012647378
    ## 97.5%               0.1936576     0.2965390   0.00293386 -0.003572344
    ##            sef enterocyte_diameter log_enterocyte_density
    ## 2.5%  3.425866       -0.0015038563              0.3086293
    ## 97.5% 7.622589        0.0009059397              1.1674755

``` r
sapply(pos_fits$med, ci)
```

    ##       log_intestinal_diameter villus_height villus_width  crypt_width
    ## 2.5%               -0.1984294    0.05459952  -0.02886595 -0.013205706
    ## 97.5%               0.2049461    0.25264785  -0.00524069 -0.001736775
    ##            sef enterocyte_diameter log_enterocyte_density
    ## 2.5%  2.589389       -2.339798e-03             -0.9042551
    ## 97.5% 7.925245       -1.515785e-05              2.2500809

``` r
sapply(pos_fits$prox, ci)
```

    ##       log_intestinal_diameter villus_height villus_width crypt_width
    ## 2.5%               -0.1338032   -0.02618653  -0.04749609 -0.03151140
    ## 97.5%               0.2469784    0.20931842  -0.01601584  0.01036664
    ##            sef enterocyte_diameter log_enterocyte_density
    ## 2.5%  1.552915       -0.0017673248              0.1512608
    ## 97.5% 7.471908        0.0005472164              1.0464692

``` r
sapply(pos_fits$dist[c('log_intestinal_diameter', 'villus_height')], 
       ci, parameter = 'log_mass')
```

    ##       log_intestinal_diameter villus_height
    ## 2.5%                0.1336546    0.01753835
    ## 97.5%               0.4005606    0.12853220

``` r
sapply(pos_fits$med[c('log_intestinal_diameter', 'villus_height')], 
       ci, parameter = 'log_mass')
```

    ##       log_intestinal_diameter villus_height
    ## 2.5%               0.05162472    0.01717659
    ## 97.5%              0.32610715    0.15016828

``` r
sapply(pos_fits$prox[c('log_intestinal_diameter', 'villus_height')], 
       ci, parameter = 'log_mass')
```

    ##       log_intestinal_diameter villus_height
    ## 2.5%                0.1625201    0.07717142
    ## 97.5%               0.4200719    0.23395167

`Clearance ~ SEF`
=================

Clearance = "paracellular probe L-arabinose clearance"

Both are log-transformed.

From the original manuscript: &gt; ... we used reduced major axis regression (model II regression)... because both &gt; variables \[X and Y\] were subject to error

Instead of an RMA regression, I'll be using a modified version of `ape::corphylo` to estimate Pearson correlation coefficients. Confidence intervals are calculated using Fisher information.

``` r
clear_df <- get_df('clear')
clear_se_df <- get_df('clear', .stat = 'se')  # <-- contains standard errors
clear_tr <- get_tr('clear')

Xmat <- cbind(clear_df$log_sef, clear_df$log_clear)
rownames(Xmat) <- rownames(clear_df)

MEmat <- cbind(clear_se_df$log_sef, clear_se_df$log_clear)
rownames(MEmat) <- clear_se_df$species

clear_cor <- corp(Xmat, phy = clear_tr, SeM = MEmat)

# Correlation with 95% CI
clear_cor['r',]
```

    ##       value    upper      lower
    ## r 0.5772012 1.066615 0.08778754

Assembling all output into one object
=====================================

I ran `ci_df` on all models above. I also added the output from the correlation between `log_sef` and `log_clear` manually at the end, which I have to do because it isn't a `phylolm` object.

``` r
mod_summaries <- bind_rows(
    list(
        ci_df(diet_fit),
        ci_df(absorp_fit),
        bind_rows(lapply(spp_fits, ci_df)),
        bind_rows(
            lapply(names(pos_fits), function(p) {
                bind_rows(lapply(pos_fits[[p]], ci_df, .pos = p))
            })),
        bind_cols(X = 'log_clear', Y = 'log_sef', clear_cor['r',])
        ))
knitr::kable(mod_summaries, format = 'markdown')
```

| Y                         | X              | pos  |       value|       lower|       upper|
|:--------------------------|:---------------|:-----|-----------:|-----------:|-----------:|
| sef                       | dietOmnivorous | NA   |  -2.3009901|  -5.7082512|   1.2984528|
| sef                       | dietProtein    | NA   |  -1.0581114|  -4.1793189|   2.2992833|
| absorp                    | taxonBat       | NA   |   0.0466086|   0.0328951|   0.0609565|
| int\_length\_mass         | taxonBat       | NA   |  -2.8583328|  -4.9166135|  -0.6899693|
| nsa\_mass                 | taxonBat       | NA   |  -0.4367081|  -0.7478906|  -0.1411825|
| vill\_area\_mass          | taxonBat       | NA   |   0.9414326|  -2.8676465|   4.9878776|
| log\_total\_enterocytes   | taxonBat       | NA   |   0.2853143|  -0.2823279|   0.8362962|
| log\_total\_enterocytes   | log\_mass      | NA   |   0.7165538|   0.3237234|   1.0874341|
| log\_intestinal\_diameter | taxonBat       | prox |   0.0612330|  -0.1338032|   0.2469784|
| log\_intestinal\_diameter | log\_mass      | prox |   0.2919863|   0.1625201|   0.4200719|
| villus\_height            | taxonBat       | prox |   0.0887347|  -0.0261865|   0.2093184|
| villus\_height            | log\_mass      | prox |   0.1560026|   0.0771714|   0.2339517|
| villus\_width             | taxonBat       | prox |  -0.0316416|  -0.0474961|  -0.0160158|
| crypt\_width              | taxonBat       | prox |  -0.0104877|  -0.0315114|   0.0103666|
| sef                       | taxonBat       | prox |   4.4720094|   1.5529154|   7.4719083|
| enterocyte\_diameter      | taxonBat       | prox |  -0.0006223|  -0.0017673|   0.0005472|
| log\_enterocyte\_density  | taxonBat       | prox |   0.5985526|   0.1512608|   1.0464692|
| log\_intestinal\_diameter | taxonBat       | med  |   0.0079669|  -0.1984294|   0.2049461|
| log\_intestinal\_diameter | log\_mass      | med  |   0.1927109|   0.0516247|   0.3261071|
| villus\_height            | taxonBat       | med  |   0.1526178|   0.0545995|   0.2526478|
| villus\_height            | log\_mass      | med  |   0.0843465|   0.0171766|   0.1501683|
| villus\_width             | taxonBat       | med  |  -0.0166444|  -0.0288660|  -0.0052407|
| crypt\_width              | taxonBat       | med  |  -0.0073315|  -0.0132057|  -0.0017368|
| sef                       | taxonBat       | med  |   5.2512236|   2.5893885|   7.9252450|
| enterocyte\_diameter      | taxonBat       | med  |  -0.0011836|  -0.0023398|  -0.0000152|
| log\_enterocyte\_density  | taxonBat       | med  |   0.6254613|  -0.9042551|   2.2500809|
| log\_intestinal\_diameter | taxonBat       | dist |   0.0109267|  -0.1690558|   0.1936576|
| log\_intestinal\_diameter | log\_mass      | dist |   0.2637575|   0.1336546|   0.4005606|
| villus\_height            | taxonBat       | dist |   0.2094520|   0.1292392|   0.2965390|
| villus\_height            | log\_mass      | dist |   0.0737162|   0.0175383|   0.1285322|
| villus\_width             | taxonBat       | dist |  -0.0095327|  -0.0230823|   0.0029339|
| crypt\_width              | taxonBat       | dist |  -0.0082673|  -0.0126474|  -0.0035723|
| sef                       | taxonBat       | dist |   5.5291580|   3.4258656|   7.6225894|
| enterocyte\_diameter      | taxonBat       | dist |  -0.0002774|  -0.0015039|   0.0009059|
| log\_enterocyte\_density  | taxonBat       | dist |   0.7402773|   0.3086293|   1.1674755|
| log\_sef                  | log\_clear     | NA   |   0.5772012|   0.0877875|   1.0666149|

I lastly write this summary to a csv file.

``` r
write_csv(mod_summaries, 'output/model_summaries.csv')
```

Session info
============

This outlines the package versions I used for these analyses.

    ## Session info -------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.4.1 (2017-06-30)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2017-09-23

    ## Packages -----------------------------------------------------------------

    ##  package    * version date       source        
    ##  ape        * 4.1     2017-02-14 CRAN (R 3.4.0)
    ##  assertthat   0.2.0   2017-04-11 CRAN (R 3.4.0)
    ##  backports    1.1.0   2017-05-22 CRAN (R 3.4.0)
    ##  base       * 3.4.1   2017-07-07 local         
    ##  bindr        0.1     2016-11-13 CRAN (R 3.4.0)
    ##  bindrcpp   * 0.2     2017-06-17 CRAN (R 3.4.0)
    ##  compiler     3.4.1   2017-07-07 local         
    ##  datasets   * 3.4.1   2017-07-07 local         
    ##  devtools     1.13.3  2017-08-02 CRAN (R 3.4.1)
    ##  digest       0.6.12  2017-01-27 CRAN (R 3.4.0)
    ##  dplyr      * 0.7.3   2017-09-09 CRAN (R 3.4.1)
    ##  evaluate     0.10.1  2017-06-24 CRAN (R 3.4.1)
    ##  glue         1.1.1   2017-06-21 CRAN (R 3.4.0)
    ##  graphics   * 3.4.1   2017-07-07 local         
    ##  grDevices  * 3.4.1   2017-07-07 local         
    ##  grid         3.4.1   2017-07-07 local         
    ##  highr        0.6     2016-05-09 cran (@0.6)   
    ##  hms          0.3     2016-11-22 CRAN (R 3.4.0)
    ##  htmltools    0.3.6   2017-04-28 cran (@0.3.6) 
    ##  knitr        1.17    2017-08-10 CRAN (R 3.4.1)
    ##  lattice      0.20-35 2017-03-25 CRAN (R 3.4.1)
    ##  magrittr     1.5     2014-11-22 CRAN (R 3.4.0)
    ##  memoise      1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods    * 3.4.1   2017-07-07 local         
    ##  nlme         3.1-131 2017-02-06 CRAN (R 3.4.1)
    ##  parallel     3.4.1   2017-07-07 local         
    ##  phylolm    * 2.5     2016-10-17 CRAN (R 3.4.0)
    ##  pkgconfig    2.0.1   2017-03-21 CRAN (R 3.4.0)
    ##  purrr        0.2.3   2017-08-02 CRAN (R 3.4.1)
    ##  R6           2.2.2   2017-06-17 CRAN (R 3.4.0)
    ##  Rcpp         0.12.12 2017-07-15 CRAN (R 3.4.1)
    ##  readr      * 1.1.1   2017-05-16 CRAN (R 3.4.0)
    ##  rlang        0.1.2   2017-08-09 CRAN (R 3.4.1)
    ##  rmarkdown    1.6     2017-06-15 CRAN (R 3.4.0)
    ##  rprojroot    1.2     2017-01-16 cran (@1.2)   
    ##  stats      * 3.4.1   2017-07-07 local         
    ##  stringi      1.1.5   2017-04-07 CRAN (R 3.4.0)
    ##  stringr      1.2.0   2017-02-18 CRAN (R 3.4.0)
    ##  tibble       1.3.4   2017-08-22 CRAN (R 3.4.1)
    ##  tidyr      * 0.7.1   2017-09-01 CRAN (R 3.4.1)
    ##  tools        3.4.1   2017-07-07 local         
    ##  utils      * 3.4.1   2017-07-07 local         
    ##  withr        2.0.0   2017-07-28 CRAN (R 3.4.1)
    ##  yaml         2.1.14  2016-11-12 cran (@2.1.14)
