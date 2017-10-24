phylolm results
================
Lucas Nell
24 Oct 2017

-   [Positions combined](#positions-combined)
-   [Proximal](#proximal)
-   [Medial](#medial)
-   [Distal](#distal)
-   [SEF on Taxon and log(Mass)](#sef-on-taxon-and-logmass)

This script displays results from [`phylolm`](https://doi.org/10.1093/sysbio/syu005). Column `estimate` is the coefficient estimate for the specified `X`, while `lower` and `upper` are lower and upper bounds of the 95% confidence interval for the coefficient estimate obtained by parametric bootstrapping. Multiple `X` parameters for a given `Y` indicate that both `X`s were included in the regression with `Y`, *not* that separate regressions were performed for each `X`.

Positions combined
==================

| Y                       | X              |    estimate|       lower|       upper|
|:------------------------|:---------------|-----------:|-----------:|-----------:|
| absorp                  | taxon          |   0.0466086|   0.0328951|   0.0609565|
| int\_length\_mass       | taxon          |  -2.8583328|  -4.9166135|  -0.6899693|
| log\_sef                | log\_clear     |   0.5772012|   0.0877875|   1.0666149|
| log\_total\_enterocytes | log\_mass      |   0.7165538|   0.3237234|   1.0874341|
| log\_total\_enterocytes | taxon          |   0.2853143|  -0.2823279|   0.8362962|
| nsa\_mass               | taxon          |  -0.4367081|  -0.7478906|  -0.1411825|
| sef                     | dietOmnivorous |  -2.5446199|  -5.7446261|   0.6650802|
| sef                     | dietProtein    |  -1.3408502|  -4.5131826|   1.5150948|
| vill\_area\_mass        | taxon          |   0.9414326|  -2.8676465|   4.9878776|

Proximal
========

| Y                         | X         |    estimate|       lower|       upper|
|:--------------------------|:----------|-----------:|-----------:|-----------:|
| crypt\_width              | taxon     |  -0.0109879|  -0.0163081|  -0.0056231|
| enterocyte\_diameter      | taxon     |  -0.0006223|  -0.0017673|   0.0005472|
| log\_enterocyte\_density  | taxon     |   0.5985526|   0.1512608|   1.0477713|
| log\_intestinal\_diameter | log\_mass |   0.2919863|   0.1625201|   0.4200719|
| log\_intestinal\_diameter | taxon     |   0.0612330|  -0.1338032|   0.2469784|
| sef                       | taxon     |   4.4720094|   1.5529154|   7.4719083|
| villus\_height            | log\_mass |   0.1560026|   0.0771714|   0.2339517|
| villus\_height            | taxon     |   0.0887347|  -0.0261865|   0.2093184|
| villus\_width             | taxon     |  -0.0316416|  -0.0474961|  -0.0160158|

Medial
======

| Y                         | X         |    estimate|       lower|       upper|
|:--------------------------|:----------|-----------:|-----------:|-----------:|
| crypt\_width              | taxon     |  -0.0073315|  -0.0132057|  -0.0017368|
| enterocyte\_diameter      | taxon     |  -0.0011836|  -0.0023398|  -0.0000152|
| log\_enterocyte\_density  | taxon     |   0.7738868|   0.3659795|   1.1710685|
| log\_intestinal\_diameter | log\_mass |   0.1927109|   0.0516247|   0.3261071|
| log\_intestinal\_diameter | taxon     |   0.0079669|  -0.1984294|   0.2049461|
| sef                       | taxon     |   5.2512236|   2.5893885|   7.9252450|
| villus\_height            | log\_mass |   0.0843465|   0.0171766|   0.1501683|
| villus\_height            | taxon     |   0.1526178|   0.0545995|   0.2526478|
| villus\_width             | taxon     |  -0.0166444|  -0.0288660|  -0.0052407|

Distal
======

| Y                         | X         |    estimate|       lower|       upper|
|:--------------------------|:----------|-----------:|-----------:|-----------:|
| crypt\_width              | taxon     |  -0.0082673|  -0.0126474|  -0.0035723|
| enterocyte\_diameter      | taxon     |  -0.0002774|  -0.0015039|   0.0009059|
| log\_enterocyte\_density  | taxon     |   0.7402773|   0.3086293|   1.1674755|
| log\_intestinal\_diameter | log\_mass |   0.2637575|   0.1336546|   0.4005606|
| log\_intestinal\_diameter | taxon     |   0.0109267|  -0.1690558|   0.1936576|
| sef                       | taxon     |   5.5291580|   3.4258656|   7.6225894|
| villus\_height            | log\_mass |   0.0737162|   0.0175383|   0.1285322|
| villus\_height            | taxon     |   0.2094520|   0.1292392|   0.2965390|
| villus\_width             | taxon     |  -0.0095327|  -0.0230823|   0.0029339|

SEF on Taxon and log(Mass)
==========================

This is to determine whether there's an effect of body mass on SEF. It appears there is not.

``` r
suppressPackageStartupMessages(library(phylolm))
invisible(sapply(list.files('R', '*.R', full.names = TRUE), source))
tr <- get_tr('spp')
spp_df <- get_df('spp')
set.seed(940318092)
mod <- suppressWarnings(
    phylolm(sef ~ taxon + log_mass, data = spp_df, phy = tr, model = "lambda", 
            boot = 2000))
ci(mod, 'log_mass')
```

    ##       2.5%      97.5% 
    ## -0.2788004  3.2372969
