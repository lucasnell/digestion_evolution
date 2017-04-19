Phylogenetic linear regression
================
Lucas Nell
19 Apr 2017

-   [Morphometric measurements](#morphometric-measurements)
-   [Phylogenetic tree](#phylogenetic-tree)
    -   [Visualizing tree](#visualizing-tree)
-   [Fitting phylogenetic linear regression models](#fitting-phylogenetic-linear-regression-models)
-   [Model output](#model-output)
    -   [Summaries](#summaries)
    -   [P-values](#p-values)
-   [Session info](#session-info)

Loading packages:

``` r
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(magrittr)
    library(phylolm)
    library(ape)
    library(ggplot2)
    library(ggtree)
})
```

Morphometric measurements
=========================

Cleaning `./data/morphometrics.csv` file for use and providing a useful function to retrieve columns from it

``` r
source('tidy_csv.R')
str(morph_df)
```

    ## Classes 'tbl_df', 'tbl' and 'data.frame':    1479 obs. of  7 variables:
    ##  $ diet   : chr  "Omnivorous" "Omnivorous" "Omnivorous" "Omnivorous" ...
    ##  $ taxon  : chr  "Rodent" "Rodent" "Rodent" "Rodent" ...
    ##  $ species: chr  "Akodon montensis" "Akodon montensis" "Akodon montensis" "Akodon montensis" ...
    ##  $ id     : chr  "Ak1" "Ak2" "Ak3" "Ak4" ...
    ##  $ measure: chr  "crypt width" "crypt width" "crypt width" "crypt width" ...
    ##  $ pos    : chr  "prox" "prox" "prox" "prox" ...
    ##  $ value  : num  0.0386 0.0492 0.0373 0.037 0.0236 ...

All measures found in `morph_df`:

-   `crypt width`
-   `enterocyte density`
-   `enterocyte width`
-   `intestinal diameter`
-   `intestinal length`
-   `mass`
-   `nsa`
-   `sef`
-   `villa surface area`
-   `villus height`
-   `villus width`

I am only using three: `nsa`, `sef`, and `mass`. Now I create a data frame with just these columns and their log-transformed versions. Species names are row names because `phylolm` requires that.

``` r
sp_df <- prep_df(measures = c('nsa', 'sef', 'mass'))
str(sp_df)
```

    ## 'data.frame':    18 obs. of  9 variables:
    ##  $ diet    : chr  "Herbivorous" "Herbivorous" "Protein" "Protein" ...
    ##  $ taxon   : chr  "Bat" "Bat" "Bat" "Bat" ...
    ##  $ species : chr  "Artibeus lituratus" "Carollia perspicillata" "Desmodus rotundus" "Eumops glaucinus" ...
    ##  $ nsa_log : num  3.33 1.86 2.33 2.12 1.74 ...
    ##  $ sef_log : num  2.82 2.83 2.51 2.8 2.52 ...
    ##  $ mass_log: num  4.24 2.85 3.65 3.53 2.59 ...
    ##  $ nsa     : num  28.2 6.48 10.43 8.32 5.79 ...
    ##  $ sef     : num  17.2 17.6 12.4 16.9 12.6 ...
    ##  $ mass    : num  69.6 17.9 38.5 34.1 13.5 ...

Phylogenetic tree
=================

Reading phylogenetic tree, cleaning species names, and removing unnecessary species from it

``` r
tr <- read.tree('./data/tree.nwk')
tr$tip.label <- gsub('_', ' ', tr$tip.label)
tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% (morph_df$species %>% unique)])
tr
```

    ## 
    ## Phylogenetic tree with 18 tips and 17 internal nodes.
    ## 
    ## Tip labels:
    ##  Carollia perspicillata, Artibeus lituratus, Desmodus rotundus, Eptesicus fuscus, Myotis lucifugus, Molossus rufus, ...
    ## Node labels:
    ##  , '8', '13', '14', '22', '11', ...
    ## 
    ## Rooted; includes branch lengths.

Visualizing tree
----------------

Here is the phylogenetic tree with log(NSA) as tip color and log(SEF) as tip size.

![](phylo_regr_files/figure-markdown_github/phylo_plot-1.png)

Fitting phylogenetic linear regression models
=============================================

Below fits phylogenetic linear regression models using `phylolm::phylolm`. For both `nsa` and `sef` (log-transformed), I fit models using log-transformed mass and taxon (a factor based on whether that species is a rodent or bat) as covariates. (I tried including the interaction between mass and taxon, but it increased the AIC in all models.)

I fit two types of phylogenetic-covariance models for both `nsa` and `sef` regression models: "the Ornstein-Uhlenbeck model with an ancestral state to be estimated at the root (OUfixedRoot) ... \[and\] Pagel's lambda model." (see `phylolm` documentation) As you can see from the results below, the covariance model had little effect on our conclusions.

I also ran 2,000 parametric bootstrap replicates to estimate model parameters.

The code below takes ~10 minutes to run.

``` r
set.seed(352)
nsa_fits <- lapply(c('lambda', 'OUfixedRoot'), 
                  function(m) {
                      phylolm(nsa_log ~ mass_log + taxon, data = sp_df, phy = tr,
                              model = m, boot = 2000, 
                              upper.bound = ifelse(m == 'lambda', 1.2, Inf))})
names(nsa_fits) <- c('lambda', 'ou')
sef_fits <- lapply(c('lambda', 'OUfixedRoot'), 
                   function(m) {
                       phylolm(sef_log ~ mass_log + taxon, data = sp_df, phy = tr,
                               model = m, boot = 2000,
                               upper.bound = ifelse(m == 'lambda', 1.2, Inf))})
names(sef_fits) <- c('lambda', 'ou')
save(nsa_fits, sef_fits, file = './data/model_fits.RData', compress = FALSE)
```

Model output
============

Summaries
---------

### `nsa`

*Pagel's lambda*

    ## 
    ## Call:
    ## phylolm(formula = nsa_log ~ mass_log + taxon, data = sp_df, phy = tr, 
    ##     model = m, upper.bound = ifelse(m == "lambda", 1.2, Inf), 
    ##     boot = 2000)
    ## 
    ##    AIC logLik 
    ## 10.826 -0.413 
    ## 
    ## Raw residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.46973 -0.16573 -0.01707  0.14016  0.50960 
    ## 
    ## Mean tip height: 96.46239
    ## Parameter estimate(s) using ML:
    ## lambda : 1e-07
    ## sigma2: 0.0006354724 
    ## 
    ## Coefficients:
    ##             Estimate   StdErr  t.value lowerbootCI upperbootCI  p.value
    ## (Intercept)  0.46555  0.33661  1.38305    -0.14208      1.0466 0.186891
    ## mass_log     0.55616  0.10451  5.32140     0.37210      0.7421 8.55e-05
    ## taxonRodent  0.53124  0.15099  3.51833     0.26848      0.8001 0.003105
    ##                
    ## (Intercept)    
    ## mass_log    ***
    ## taxonRodent ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Note: p-values are conditional on lambda=1e-07.
    ## 
    ## sigma2: 0.0006354724
    ##       bootstrap mean: 0.0005418981 (on raw scale)
    ##                       0.0005066157 (on log scale, then back transformed)
    ##       bootstrap 95% CI: (0.0002412875,0.001007715)
    ## 
    ## lambda: 1e-07
    ##       bootstrap mean: 0.003704019 (on raw scale)
    ##       bootstrap 95% CI: (1e-07,1e-07)
    ## 
    ## Parametric bootstrap results based on 2000 fitted replicates

*Ornstein-Uhlenbeck*

    ## 
    ## Call:
    ## phylolm(formula = nsa_log ~ mass_log + taxon, data = sp_df, phy = tr, 
    ##     model = m, upper.bound = ifelse(m == "lambda", 1.2, Inf), 
    ##     boot = 2000)
    ## 
    ##     AIC  logLik 
    ## 10.5299 -0.2649 
    ## 
    ## Raw residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.46184 -0.18192 -0.02546  0.15220  0.47670 
    ## 
    ## Mean tip height: 96.46239
    ## Parameter estimate(s) using ML:
    ## alpha: 0.1327371
    ## sigma2: 0.01640746 
    ## 
    ## Coefficients:
    ##             Estimate   StdErr  t.value lowerbootCI upperbootCI   p.value
    ## (Intercept)  0.40545  0.32698  1.24002    -0.16159      0.9678  0.234018
    ## mass_log     0.57810  0.10074  5.73861     0.39858      0.7580 3.919e-05
    ## taxonRodent  0.50513  0.15504  3.25813     0.22918      0.7919  0.005294
    ##                
    ## (Intercept)    
    ## mass_log    ***
    ## taxonRodent ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Note: p-values are conditional on alpha=0.1327371.
    ## 
    ## sigma2: 0.01640746
    ##       bootstrap mean: 2.680095 (on raw scale)
    ##                       0.05420355 (on log scale, then back transformed)
    ##       bootstrap 95% CI: (0.002771107,3.214494)
    ## 
    ## alpha: 0.1327371
    ##       bootstrap mean: 32.21303 (on raw scale)
    ##                       0.567499 (on log scale, then back transformed)
    ##       bootstrap 95% CI: (0.03313966,33.07642)
    ## 
    ## Parametric bootstrap results based on 2000 fitted replicates

### `sef`

*Pagel's lambda*

    ## 
    ## Call:
    ## phylolm(formula = sef_log ~ mass_log + taxon, data = sp_df, phy = tr, 
    ##     model = m, upper.bound = ifelse(m == "lambda", 1.2, Inf), 
    ##     boot = 2000)
    ## 
    ##    AIC logLik 
    ##  1.593  4.204 
    ## 
    ## Raw residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.28739 -0.21506  0.05644  0.10672  0.32298 
    ## 
    ## Mean tip height: 96.46239
    ## Parameter estimate(s) using ML:
    ## lambda : 1e-07
    ## sigma2: 0.000380469 
    ## 
    ## Coefficients:
    ##              Estimate    StdErr   t.value lowerbootCI upperbootCI
    ## (Intercept)  2.289672  0.260456  8.791011    1.816941      2.7541
    ## mass_log     0.113560  0.080869  1.404242   -0.032224      0.2587
    ## taxonRodent -0.519852  0.116833 -4.449527   -0.733102     -0.3149
    ##               p.value    
    ## (Intercept) 2.646e-07 ***
    ## mass_log    0.1806133    
    ## taxonRodent 0.0004682 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Note: p-values are conditional on lambda=1e-07.
    ## 
    ## sigma2: 0.000380469
    ##       bootstrap mean: 0.0003149273 (on raw scale)
    ##                       0.0002940711 (on log scale, then back transformed)
    ##       bootstrap 95% CI: (0.0001286653,0.000581292)
    ## 
    ## lambda: 1e-07
    ##       bootstrap mean: 0.002138407 (on raw scale)
    ##       bootstrap 95% CI: (1e-07,1e-07)
    ## 
    ## Parametric bootstrap results based on 2000 fitted replicates

*Ornstein-Uhlenbeck*

    ## 
    ## Call:
    ## phylolm(formula = sef_log ~ mass_log + taxon, data = sp_df, phy = tr, 
    ##     model = m, upper.bound = ifelse(m == "lambda", 1.2, Inf), 
    ##     boot = 2000)
    ## 
    ##    AIC logLik 
    ##  1.593  4.204 
    ## 
    ## Raw residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.28739 -0.21506  0.05644  0.10672  0.32298 
    ## 
    ## Mean tip height: 96.46239
    ## Parameter estimate(s) using ML:
    ## alpha: 4.128654
    ## sigma2: 0.303051 
    ## 
    ## Coefficients:
    ##              Estimate    StdErr   t.value lowerbootCI upperbootCI
    ## (Intercept)  2.289672  0.260456  8.791011    1.839643      2.7707
    ## mass_log     0.113560  0.080869  1.404242   -0.031936      0.2581
    ## taxonRodent -0.519852  0.116833 -4.449527   -0.716238     -0.3072
    ##               p.value    
    ## (Intercept) 2.646e-07 ***
    ## mass_log    0.1806133    
    ## taxonRodent 0.0004682 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Note: p-values are conditional on alpha=4.128654.
    ## 
    ## sigma2: 0.303051
    ##       bootstrap mean: 7.850735 (on raw scale)
    ##                       0.08702424 (on log scale, then back transformed)
    ##       bootstrap 95% CI: (0.001814778,8.089423)
    ## 
    ## alpha: 4.128654
    ##       bootstrap mean: 105.323 (on raw scale)
    ##                       1.539463 (on log scale, then back transformed)
    ##       bootstrap 95% CI: (0.04052872,126.7784)
    ## 
    ## Parametric bootstrap results based on 2000 fitted replicates

P-values
--------

These are p-value based on bootstrap replicates for whether the coefficient for the taxon covariate is significantly different from zero.

### `nsa`

    ## P for Pagel's lambda     = 0.001

    ## P for Ornstein-Uhlenbeck = 0.001

### `sef`

    ## P for Pagel's lambda     = 0.000

    ## P for Ornstein-Uhlenbeck = 0.000

Session info
============

This outlines the package versions I used for these analyses.

``` r
devtools::session_info()
```

    ## Session info --------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.3.3 (2017-03-06)
    ##  system   x86_64, darwin13.4.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2017-04-19

    ## Packages ------------------------------------------------------------------

    ##  package    * version date       source        
    ##  ape        * 4.1     2017-02-14 CRAN (R 3.3.2)
    ##  assertthat   0.2.0   2017-04-11 CRAN (R 3.3.2)
    ##  backports    1.0.5   2017-01-18 CRAN (R 3.3.2)
    ##  colorspace   1.3-2   2016-12-14 CRAN (R 3.3.2)
    ##  DBI          0.6-1   2017-04-01 CRAN (R 3.3.2)
    ##  devtools     1.12.0  2016-06-24 CRAN (R 3.3.0)
    ##  digest       0.6.12  2017-01-27 CRAN (R 3.3.2)
    ##  dplyr      * 0.5.0   2016-06-24 CRAN (R 3.3.0)
    ##  evaluate     0.10    2016-10-11 CRAN (R 3.3.0)
    ##  ggplot2    * 2.2.1   2016-12-30 CRAN (R 3.3.2)
    ##  ggtree     * 1.6.11  2017-03-15 Bioconductor  
    ##  gtable       0.2.0   2016-02-26 CRAN (R 3.3.0)
    ##  hms          0.3     2016-11-22 CRAN (R 3.3.2)
    ##  htmltools    0.3.5   2016-03-21 CRAN (R 3.3.0)
    ##  jsonlite     1.4     2017-04-08 CRAN (R 3.3.2)
    ##  knitr        1.15.1  2016-11-22 CRAN (R 3.3.2)
    ##  labeling     0.3     2014-08-23 CRAN (R 3.3.0)
    ##  lattice      0.20-35 2017-03-25 CRAN (R 3.3.2)
    ##  lazyeval     0.2.0   2016-06-12 CRAN (R 3.3.0)
    ##  magrittr   * 1.5     2014-11-22 CRAN (R 3.3.0)
    ##  memoise      1.0.0   2016-01-29 CRAN (R 3.3.0)
    ##  munsell      0.4.3   2016-02-13 CRAN (R 3.3.0)
    ##  nlme         3.1-131 2017-02-06 CRAN (R 3.3.3)
    ##  phylolm    * 2.5     2016-10-17 CRAN (R 3.3.0)
    ##  plyr         1.8.4   2016-06-08 CRAN (R 3.3.0)
    ##  R6           2.2.0   2016-10-05 CRAN (R 3.3.0)
    ##  Rcpp         0.12.10 2017-03-19 CRAN (R 3.3.2)
    ##  readr      * 1.1.0   2017-03-22 CRAN (R 3.3.2)
    ##  rmarkdown    1.4     2017-03-24 CRAN (R 3.3.2)
    ##  rprojroot    1.2     2017-01-16 CRAN (R 3.3.2)
    ##  scales       0.4.1   2016-11-09 CRAN (R 3.3.2)
    ##  stringi      1.1.5   2017-04-07 CRAN (R 3.3.2)
    ##  stringr      1.2.0   2017-02-18 CRAN (R 3.3.2)
    ##  tibble       1.3.0   2017-04-01 CRAN (R 3.3.2)
    ##  tidyr      * 0.6.1   2017-01-10 CRAN (R 3.3.2)
    ##  withr        1.0.2   2016-06-20 CRAN (R 3.3.0)
    ##  yaml         2.1.14  2016-11-12 CRAN (R 3.3.2)
