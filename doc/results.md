phylolm results
================
Lucas Nell
25 Oct 2017

-   [Coefficient estimates](#coefficient-estimates)
-   [Phylogenetic signal](#phylogenetic-signal)
-   [SEF on Taxon and log(Mass)](#sef-on-taxon-and-logmass)
-   [Percent differences](#percent-differences)

This script displays results from analyses using [`phylolm`](https://doi.org/10.1093/sysbio/syu005). Column `estimate` is the coefficient estimate for the specified `X`, while `lower` and `upper` are lower and upper bounds of the 95% confidence interval for the coefficient estimate obtained by parametric bootstrapping. Multiple `X` parameters for a given `Y` indicate that both `X`s were included in the regression with `Y`, *not* that separate regressions were performed for each `X`.

> This page uses the variable name `clade` for the dummy variable indicating whether a species is a bat, but in the rest of the scripts in this repo, this variable is referred to as `taxon`.

Coefficient estimates
=====================

Positions combined
------------------

| Y                       | X              | estimate | lower   | upper   |
|:------------------------|:---------------|:---------|:--------|:--------|
| absorp                  | clade          | 0.04661  | 0.0329  | 0.06096 |
| int\_length\_mass       | clade          | -2.858   | -4.917  | -0.69   |
| log\_sef                | log\_clear     | 0.5772   | 0.08779 | 1.067   |
| log\_total\_enterocytes | clade          | 0.2853   | -0.2823 | 0.8363  |
| log\_total\_enterocytes | log\_mass      | 0.7166   | 0.3237  | 1.087   |
| nsa\_mass               | clade          | -0.4367  | -0.7479 | -0.1412 |
| sef                     | dietOmnivorous | -2.545   | -5.745  | 0.6651  |
| sef                     | dietProtein    | -1.341   | -4.513  | 1.515   |
| vill\_area\_mass        | clade          | 0.9414   | -2.868  | 4.988   |

Proximal
--------

| Y                         | X         | estimate   | lower     | upper     |
|:--------------------------|:----------|:-----------|:----------|:----------|
| crypt\_width              | clade     | -0.01099   | -0.01631  | -0.005623 |
| enterocyte\_diameter      | clade     | -0.0006223 | -0.001767 | 0.0005472 |
| log\_enterocyte\_density  | clade     | 0.5986     | 0.1513    | 1.048     |
| log\_intestinal\_diameter | clade     | 0.06123    | -0.1338   | 0.247     |
| log\_intestinal\_diameter | log\_mass | 0.292      | 0.1625    | 0.4201    |
| sef                       | clade     | 4.472      | 1.553     | 7.472     |
| villus\_height            | clade     | 0.08873    | -0.02619  | 0.2093    |
| villus\_height            | log\_mass | 0.156      | 0.07717   | 0.234     |
| villus\_width             | clade     | -0.03164   | -0.0475   | -0.01602  |

Medial
------

| Y                         | X         | estimate  | lower    | upper      |
|:--------------------------|:----------|:----------|:---------|:-----------|
| crypt\_width              | clade     | -0.007331 | -0.01321 | -0.001737  |
| enterocyte\_diameter      | clade     | -0.001184 | -0.00234 | -1.516e-05 |
| log\_enterocyte\_density  | clade     | 0.7739    | 0.366    | 1.171      |
| log\_intestinal\_diameter | clade     | 0.007967  | -0.1984  | 0.2049     |
| log\_intestinal\_diameter | log\_mass | 0.1927    | 0.05162  | 0.3261     |
| sef                       | clade     | 5.251     | 2.589    | 7.925      |
| villus\_height            | clade     | 0.1526    | 0.0546   | 0.2526     |
| villus\_height            | log\_mass | 0.08435   | 0.01718  | 0.1502     |
| villus\_width             | clade     | -0.01664  | -0.02887 | -0.005241  |

Distal
------

| Y                         | X         | estimate   | lower     | upper     |
|:--------------------------|:----------|:-----------|:----------|:----------|
| crypt\_width              | clade     | -0.008267  | -0.01265  | -0.003572 |
| enterocyte\_diameter      | clade     | -0.0002774 | -0.001504 | 0.0009059 |
| log\_enterocyte\_density  | clade     | 0.7403     | 0.3086    | 1.167     |
| log\_intestinal\_diameter | clade     | 0.01093    | -0.1691   | 0.1937    |
| log\_intestinal\_diameter | log\_mass | 0.2638     | 0.1337    | 0.4006    |
| sef                       | clade     | 5.529      | 3.426     | 7.623     |
| villus\_height            | clade     | 0.2095     | 0.1292    | 0.2965    |
| villus\_height            | log\_mass | 0.07372    | 0.01754   | 0.1285    |
| villus\_width             | clade     | -0.009533  | -0.02308  | 0.002934  |

Phylogenetic signal
===================

Positions combined
------------------

| Y                       | X                 | lambda                 | sigma2                             |
|:------------------------|:------------------|:-----------------------|:-----------------------------------|
| absorp                  | clade             | 1e-07 \[1e-07, 1e-07\] | 1.069e-06 \[1.723e-07, 1.962e-06\] |
| int\_length\_mass       | clade             | 1e-07 \[1e-07, 1e-07\] | 0.05219 \[0.02052, 0.08702\]       |
| log\_total\_enterocytes | clade + log\_mass | 1e-07 \[1e-07, 1e-07\] | 0.002756 \[0.0009931, 0.004587\]   |
| nsa\_mass               | clade             | 1e-07 \[1e-07, 1e-07\] | 0.001104 \[0.0004277, 0.001804\]   |
| sef                     | diet              | 0.7184 \[1e-07, 1\]    | 0.1088 \[0.02197, 0.2035\]         |
| vill\_area\_mass        | clade             | 1e-07 \[1e-07, 1e-07\] | 0.1973 \[0.07853, 0.3488\]         |

Proximal
--------

| Y                         | X                 | lambda                 | sigma2                             |
|:--------------------------|:------------------|:-----------------------|:-----------------------------------|
| crypt\_width              | clade             | 1e-07 \[1e-07, 1e-07\] | 3.595e-07 \[1.395e-07, 5.682e-07\] |
| enterocyte\_diameter      | clade             | 1e-07 \[1e-07, 1e-07\] | 1.623e-08 \[6.276e-09, 2.725e-08\] |
| log\_enterocyte\_density  | clade             | 1e-07 \[1e-07, 1e-07\] | 0.002594 \[0.001042, 0.004351\]    |
| log\_intestinal\_diameter | clade + log\_mass | 1e-07 \[1e-07, 1e-07\] | 0.0002917 \[0.000108, 0.00048\]    |
| sef                       | clade             | 1e-07 \[1e-07, 1e-07\] | 0.1055 \[0.04028, 0.1791\]         |
| villus\_height            | clade + log\_mass | 1e-07 \[1e-07, 1e-07\] | 0.0001122 \[3.926e-05, 0.000175\]  |
| villus\_width             | clade             | 1e-07 \[1e-07, 1e-07\] | 2.884e-06 \[1.11e-06, 4.73e-06\]   |

Medial
------

| Y                         | X                 | lambda                 | sigma2                             |
|:--------------------------|:------------------|:-----------------------|:-----------------------------------|
| crypt\_width              | clade             | 1e-07 \[1e-07, 1e-07\] | 4.01e-07 \[1.563e-07, 6.709e-07\]  |
| enterocyte\_diameter      | clade             | 1e-07 \[1e-07, 1e-07\] | 1.554e-08 \[5.808e-09, 2.576e-08\] |
| log\_enterocyte\_density  | clade             | 1e-07 \[1e-07, 1e-07\] | 0.002 \[0.0007841, 0.003197\]      |
| log\_intestinal\_diameter | clade + log\_mass | 1e-07 \[1e-07, 1e-07\] | 0.0003425 \[0.0001136, 0.000534\]  |
| sef                       | clade             | 1e-07 \[1e-07, 1e-07\] | 0.0858 \[0.0322, 0.1428\]          |
| villus\_height            | clade + log\_mass | 1e-07 \[1e-07, 1e-07\] | 8.31e-05 \[2.9e-05, 0.000133\]     |
| villus\_width             | clade             | 1e-07 \[1e-07, 1e-07\] | 1.632e-06 \[6.245e-07, 2.725e-06\] |

Distal
------

| Y                         | X                 | lambda                   | sigma2                             |
|:--------------------------|:------------------|:-------------------------|:-----------------------------------|
| crypt\_width              | clade             | 1e-07 \[1e-07, 1e-07\]   | 2.751e-07 \[1.048e-07, 4.518e-07\] |
| enterocyte\_diameter      | clade             | 1e-07 \[1e-07, 1e-07\]   | 1.635e-08 \[6.41e-09, 2.745e-08\]  |
| log\_enterocyte\_density  | clade             | 1e-07 \[1e-07, 1e-07\]   | 0.002139 \[0.0008007, 0.003781\]   |
| log\_intestinal\_diameter | clade + log\_mass | 1e-07 \[1e-07, 1e-07\]   | 0.00031 \[0.0001039, 0.0005158\]   |
| sef                       | clade             | 1e-07 \[1e-07, 1e-07\]   | 0.0576 \[0.02207, 0.0991\]         |
| villus\_height            | clade + log\_mass | 1e-07 \[1e-07, 1e-07\]   | 6.046e-05 \[2.079e-05, 9.808e-05\] |
| villus\_width             | clade             | 1e-07 \[1e-07, 0.01453\] | 1.973e-06 \[7.84e-07, 3.498e-06\]  |

SEF on Taxon and log(Mass)
==========================

This is to determine whether there's an effect of body mass on SEF. It appears there is not.

``` r
invisible(sapply(list.files('R', '*.R', full.names = TRUE), source))
tr <- get_tr('spp')
spp_df <- get_df('spp') %>% rename(clade = taxon)
set.seed(940318092)
mod <- suppressWarnings(
    phylolm(sef ~ clade + log_mass, data = spp_df, phy = tr, model = "lambda", 
            boot = 2000))
# Confidence intervals for coefficient estimates using parametric bootstrapping:
t(ci(mod, c('log_mass', 'cladeBat')))
```

    ##                2.5%    97.5%
    ## log_mass -0.2788004 3.237297
    ## cladeBat  3.4640604 8.652407

Percent differences
===================

This is for the percent differences in the discussion and abstract.

``` r
perc_diff <- function(.m, log_trans = FALSE) {
    .c <- as.numeric({.m %>% summary %>% coef}[,'Estimate'])
    if (log_trans) {
        prop_diff <- {exp(.c[2] + .c[1]) - exp(.c[1])} / exp(.c[1])
    } else {
        prop_diff <- .c[2] / .c[1]
    }
    return(prop_diff * 100)
}
```

NSA corrected for mass
----------------------

``` r
cat(sprintf('%.4g %%\n', perc_diff(models$spp$nsa_mass)))
```

    ## -32.37 %

SEF
---

``` r
cat(sprintf('%.4g %%\n', mean(c(
    perc_diff(models$pos$prox$sef),
    perc_diff(models$pos$med$sef),
    perc_diff(models$pos$dist$sef)
))))
```

    ## 57.39 %

Enterocyte density
------------------

``` r
cat(sprintf('%.4g %%\n', mean(c(
    perc_diff(models$pos$prox$log_enterocyte_density, TRUE), 
    perc_diff(models$pos$med$log_enterocyte_density, TRUE), 
    perc_diff(models$pos$dist$log_enterocyte_density, TRUE)
))))
```

    ## 102.8 %
