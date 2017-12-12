Aggregate from individual to species-level data
================
Lucas Nell
12 Dec 2017

-   [Summary of output](#summary-of-output)
-   [Morphometric data aggregated by species](#morphometric-data-aggregated-by-species)
-   [Morphometric data by species AND position](#morphometric-data-by-species-and-position)
-   [Clearance data](#clearance-data)
-   [Absorption data](#absorption-data)
-   [Saving these data frames](#saving-these-data-frames)
-   [Session info](#session-info)

This converts the cleaned Excel-file-derived data frames (with data for individual animals) into CSV files with means and standard errors by species and/or intestinal segment. It needs to be run only once, and not at all if you have the csv files already.

Sourcing the file `doc/01-xl_to_dfs.R` creates three data frames: `morph_df`, `clear_df`, and `absorp_df`. Each is associated with one or more data set(s) output from this file.

``` r
source("reports/01-xl_to_dfs.R")
```

Summary of output
=================

I created two data sets from morphometric data in `morph_df`:

1.  Measurements separated by species (found in `output/tidy_spp.csv`)
2.  Measurements separated by species and intestinal segment (`output/tidy_pos.csv`)

I also created one data set each from `clear_df` and `absorp_df`. These csvs are found in `output/tidy_clear.csv` and `output/tidy_absorp.csv`, respectively.

Each data set is associated with a set of analyses and only needs certain measurements, so only those are included.

The csv files are used in the function `get_df` in the `R/get_data.R` file to retrieve a data frame for a given analysis set. **The analysis sets are as follows:**

1.  `'spp'`: Measurements separated by species.
2.  `'pos'`: Measurements by species and position. (You have to also provide the position for this analysis set.)
3.  `'clear'`: Clearance data by species. (Uses a different set of individuals entirely from the 1–2.)
4.  `'absorp'`: Absorption data by species. (Uses a different set of individuals entirely from the 1–3.)

For the second analysis set, I need to do the analyses separately for each position because modelling within-species and within-individual variance due to position rather than process error would be difficult and not likely possible with this small dataset.

Morphometric data aggregated by species
=======================================

For the following measures, we need a single morphometric value per species:

-   Intestinal length
-   NSA
-   Villus surface area
-   Total number of enterocytes
-   Calculated as such: `log(NSA * enterocyte_density)`

> All are log-transformed to improve linearity

We need the following columns from `morph_df` to compute these values:

``` r
spp_measures <- c('mass', 'intestinal_length', 'nsa', 'vill_surface_area',
                  'enterocyte_density', 'sef')
```

Next I manipulated `morph_df` as such to get means for each parameter that I will be using for my analyses.

``` r
spp_df <- morph_df %>% 
    # Selecting measurement columns, plus the identifying columns
    select_(.dots = append(list('diet', 'clade', 'species', 'id', 'pos'), 
                           spp_measures)) %>% 
    # Removing all rows with all NAs in measures columns
    filter(Reduce(`+`, lapply(.[,spp_measures], is.na)) < length(spp_measures)) %>% 
    # Add nsa to all positions' estimates (for total_enterocytes below)
    group_by(clade, diet, species, id) %>% 
    mutate(nsa = ifelse(is.na(nsa), nsa[!is.na(nsa)], nsa)) %>% 
    ungroup %>% 
    # Doing the calculations / transformations now, before taking any means
    mutate(log_total_enterocytes = log(enterocyte_density * nsa),
           log_mass = log(mass),
           log_intestinal_length = log(intestinal_length), 
           log_vill_surface_area = log(vill_surface_area), 
           log_nsa = log(nsa),
           log_sef = log(sef)) %>% 
    select(clade, diet, species, id,
           log_sef, log_intestinal_length, log_nsa, log_vill_surface_area, 
           log_total_enterocytes, log_mass) %>% 
    # Taking mean by individual (i.e., across segments)
    group_by(clade, diet, species, id) %>% 
    summarize_all(mean, na.rm = TRUE) %>% 
    ungroup %>% 
    # Now taking mean and SE by species
    select(-id) %>% 
    group_by(clade, diet, species) %>% 
    summarize_all(funs(mean, se), na.rm = TRUE) %>% 
    # Adjusting species with NA for their SE (i.e., those spp with only 1 individual)
    mutate_at(vars(ends_with('_se')), na_se) %>%
    ungroup %>%
    arrange(clade, diet, species)
```

Lastly I add a column named `stat` to distinguish between mean and SE values, making filtering for a given statistic easier.

``` r
spp_df <- spp_df %>%
    gather('measure', 'value', -clade, -diet, -species) %>%
    mutate(stat = sapply(measure, function(s) tail(strsplit(s, '_')[[1]], 1)),
           measure = sapply(measure,
                            function(s) {
                                paste0(head(strsplit(s, '_')[[1]], -1), collapse = "_")
                            })) %>%
    spread(measure, value)
```

Morphometric data by species AND position
=========================================

For the following measures, we need 3 morphometric values per species, one per intestinal segment:

-   Intestinal diameter (log-transformed)
-   Villus height (log-transformed)
-   Villus width
-   Crypt width
-   Surface enlargement factor (SEF) (log-transformed)
-   Enterocyte diameter
-   Enterocytes per cm^2 NSA (log-transformed)

We need the following columns from `morph_df` to compute these values:

``` r
pos_measures <- c('mass', 'intestinal_diameter', 'villus_height',  'villus_width',
                  'crypt_width', 'sef', 'enterocyte_diameter', 'enterocyte_density')
```

Next I manipulated `morph_df` as such to get means and standard errors for each parameter that I will be analyzing.

``` r
pos_df <- morph_df %>% 
    select_(.dots = append(list('clade', 'diet', 'species', 'pos', 'id'), 
                           as.list(pos_measures))) %>% 
    # Add mass to all positions' estimates
    group_by(clade, diet, species, id) %>% 
    mutate(mass = ifelse(is.na(mass), mass[!is.na(mass)], mass)) %>% 
    ungroup %>% 
    # Now removing rows with pos == NA bc they don't have the other measurements
    filter(!is.na(pos)) %>% 
    # Taking the log now, before taking any means
    mutate_(.dots = setNames(as.list(sprintf('%s(%s)', 'log', pos_measures)), 
                             paste0('log_', pos_measures))) %>% 
    # Grouping by, then taking mean and SE of all measurement columns and transformed-
    # measurement columns
    group_by(clade, diet, species, pos) %>% 
    summarize_at(.vars = c(pos_measures, paste0('log_', pos_measures)), funs(mean, se), 
                 na.rm = TRUE) %>% 
    ungroup %>% 
    # Adjusting species with NA for their SE (i.e., those spp with only 1 individual)
    mutate_at(vars(ends_with('_se')), na_se) %>%
    # Now only outputting columns that are necessary
    select_(
        .dots = as.list(c('clade', 'diet', 'species', 'pos',
          paste0(rep(c('log_intestinal_diameter', 'log_villus_height', 'villus_width',
                       'crypt_width', 'log_sef', 'enterocyte_diameter',
                       'log_enterocyte_density', 'log_mass'), each = 2),
                 c('_mean', '_se'))))
    ) %>%
    arrange(clade, diet, species, pos)
```

Lastly I add a column named `stat` to distinguish between mean and SE values, making filtering for a given statistic easier.

``` r
pos_df <- pos_df %>%
    gather('measure', 'value', -clade, -diet, -species, -pos) %>%
    mutate(stat = sapply(measure, function(s) tail(strsplit(s, '_')[[1]], 1)),
           measure = sapply(measure,
                            function(s) {
                                paste0(head(strsplit(s, '_')[[1]], -1), collapse = "_")
                            })) %>%
    spread(measure, value)
```

Clearance data
==============

``` r
clear_df <- clear_df %>%
    mutate(
        # They lumped herbivores and omnivores together as "carb eater <clade>"
        diet = ifelse(diet == "Protein", diet, "Carb"),
        # Averaging SEF by individual, on log scale
        log_sef = (log(prox) + log(med) + log(dist)) / 3,
        # Taking log of clearance before any means are calculated
        # Some clearances were negative, so there would be NaNs produced
        # The next two lines are to avoid the warning message
        log_clear = ifelse(clear < 0, NA, clear),
        log_clear = log(log_clear)
    ) %>% 
    # These are no longer necessary
    select(-prox, -med, -dist, -clear) %>% 
    group_by(diet, clade, species) %>% 
    summarize_at(vars(log_sef, log_clear), funs(mean, se), na.rm = TRUE) %>% 
    ungroup
```

I'm going to add in columns `log_enterocyte_density` and `log_mass` from the morphometric data frame.

``` r
# Data frame for retrieving
retr_df <- morph_df %>% 
    select(clade, diet, species, pos, id, mass, enterocyte_density) %>% 
    # Taking the log now, before taking any means
    mutate(log_mass = log(mass), log_enterocyte_density = log(enterocyte_density)) %>% 
    # Grouping by, then taking mean by individual (i.e., among segments)
    group_by(species, id) %>% 
    summarize_at(vars(log_mass, log_enterocyte_density), mean, na.rm = TRUE) %>% 
    ungroup %>% 
    # Now getting mean and SE by species
    group_by(species) %>% 
    summarize_at(vars(log_mass, log_enterocyte_density), funs(mean, se), na.rm = TRUE) %>% 
    ungroup %>% 
    # We don't need SE for log_mass
    rename(log_mass = log_mass_mean) %>% 
    select(-log_mass_se) %>% 
    # Adjusting species with NA for their SE (i.e., those spp with only 1 individual)
    mutate_at(vars(ends_with('_se')), na_se)

clear_df <- left_join(clear_df, retr_df, by = 'species')
```

I now add a column named `stat` to distinguish between mean and SE values, making filtering for a given statistic easier.

``` r
clear_df <- clear_df %>%
    gather('measure', 'value', -diet, -clade, -species, -log_mass) %>%
    mutate(stat = sapply(measure, function(s) tail(strsplit(s, '_')[[1]], 1)),
           measure = sapply(measure,
                            function(s) {
                                paste0(head(strsplit(s, '_')[[1]], -1), collapse = "_")
                            })) %>%
    spread(measure, value)
```

Absorption data
===============

The absorption data is a little weird in that different individuals were used for the various measurements necessary to get the final parameter values for each species. Because of this, I have to combine these measurements a bit differently than before.

The final parameter, `absorp`, equals the following:

`(gavage / injection) / (nsa * sef)`

This can be manipulated to the following:

`gavage * (1 / injection) * { 1 / (nsa * sef) }`

For three species, the individuals used for measurements for this calculation are split into three groups, corresponding to the three portions of the equation above. One set of individuals were used for `gavage`, another for `injection`, and a third for `nsa` and `sef`.

These are the above-mentioned three species:

``` r
sep_absorps <- c('Myotis lucifugus', 'Tadarida brasiliensis', 'Akodon montensis')
```

For the rest of the species, `gavage` and `injection` were measured together, but `nsa` and `sef` were measured in a different set of individuals.

This is important because `mean(X/Y) != mean(X) / mean(Y)` and `var(X/Y) != var(X) / var(Y)`. But because `mean(X*Y) = mean(X) * mean(Y)` and `var(XY) = mean(X)^2 * var(Y) + mean(Y)^2 * var(X) + var(X) * var(Y)`, I can inverse things before taking means or variances to allow me to combine these values.

So before taking means and variances, I have inversed `gavage` and `(nsa * sef)`, thereby allowing me to use the equations above. For the three species with `gavage` and `injection` measured separately, I had to combine variances for `gavage` and `injection` first, then combine the variance of `gavage * injection` with `1 / (nsa * sef)`.

The term `1 / (nsa * sef)` is called `ns` below.

Lastly, I want to log-transform `absorp`. I do that as follows:

    If X is lognormally distributed with mean M and variance V, X = exp(Y),
    then the mean m and variance v of Y are
    v = log(V/M^2 +1)
    m = log(M)-v/2

``` r
absorp_df <- absorp_df %>%
    mutate(
        # Averaging SEF by individual (i.e., across segments)
        sef = (prox + med + dist) / 3,
        ns = 1 / (nsa * sef),
        log_mass = log(mass)
    ) %>% 
    group_by(diet, clade, species) %>% 
    summarize(
        
        # I'm being conservative and taking the minimum sample size among
        # any measurements.
        # The ifelse statement is to include injection in the sample size vector
        # to choose the minimum from for species that had injection measured separately.
        # I don't want to include injection for other species, because its sample size 
        # is zero.
        n = ifelse(species[1] %in% sep_absorps,
                   min(length(gavage[!is.na(gavage)]), 
                       length(injection[!is.na(injection)]),
                       length(ns[!is.na(ns)])),
                   min(length(gavage[!is.na(gavage)]), 
                       length(ns[!is.na(ns)]))),
        
        gavage_v = var(gavage, na.rm = TRUE),
        gavage = mean(gavage, na.rm = TRUE),
        
        # For non-sep_absorps species, I'm setting injection mean to 1 and variance to 0
        # bc the final value is already in the gavage column, and doing this makes the
        # equation for combined varianced above simplify to = V(X), where X is gavage.
        # It also makes the mean end up equalling just the gavage mean.
        inv_injection_v = ifelse(species[1] %in% sep_absorps,
                                 var(1 / injection, na.rm = TRUE), 0),
        inv_injection = ifelse(species[1] %in% sep_absorps,
                               mean(1 / injection, na.rm = TRUE), 1),
        
        # "gi" stands for gavage and injection
        gi_v = gavage^2 * inv_injection_v + inv_injection^2 * gavage_v +
            gavage_v * inv_injection_v,
        gi = gavage * inv_injection,
        
        ns_v = var(ns, na.rm = TRUE),
        ns = mean(ns, na.rm = TRUE),
        
        absorp_v = ns^2 * gi_v + gi^2 * ns_v + ns_v * gi_v,
        absorp_mean = gavage * inv_injection * ns,
        
        # Now I log-transform absorp, first calculating the variance bc it's necessary 
        # for both mean and se
        log_absorp_v = log(absorp_v / absorp_mean^2 + 1),
        log_absorp_se = sqrt(log_absorp_v) / sqrt(n),
        log_absorp_mean = log(absorp_mean) - log_absorp_v / 2,
        
        # I don't need SE for log_mass because I'm assuming it's without error
        log_mass = mean(log_mass, na.rm = TRUE)
    ) %>% 
    ungroup %>% 
    select(clade, species, log_absorp_mean, log_absorp_se, log_mass)


rm(sep_absorps)
```

I'm going to add in the column `log_total_enterocytes` and `log_enterocyte_density` from the morphometric data frame.

``` r
# I'm getting `log_total_enterocytes` from `spp_df` and `log_enterocyte_density`
# from `retr_df` created above
absorp_df <- left_join(absorp_df,
          spp_df %>% filter(stat == 'mean') %>% select(species, log_total_enterocytes),
          by = 'species') %>% 
    left_join(spp_df %>% filter(stat == 'se') %>% select(species, log_total_enterocytes),
          by = 'species', suffix = c('_mean', '_se')) %>% 
    left_join(retr_df %>% select(-log_mass), by = 'species')
```

As above, I add a column named `stat` to distinguish between mean and SE values, making filtering for a given statistic easier.

``` r
absorp_df <- absorp_df %>%
    gather('measure', 'value', -clade, -species, -log_mass) %>%
    mutate(stat = sapply(measure, function(s) tail(strsplit(s, '_')[[1]], 1)),
           measure = sapply(measure,
                            function(s) {
                                paste0(head(strsplit(s, '_')[[1]], -1), collapse = "_")
                            })) %>%
    spread(measure, value)
```

Saving these data frames
========================

I'm saving them as csv files so they can be quickly loaded when retrieving a data frame.

``` r
write_csv(pos_df, 'output/tidy_pos.csv')
write_csv(spp_df, 'output/tidy_spp.csv')
write_csv(clear_df, 'output/tidy_clear.csv')
write_csv(absorp_df, 'output/tidy_absorp.csv')
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
    ##  date     2017-12-12

    ## Packages -----------------------------------------------------------------

    ##  package    * version date       source        
    ##  assertthat   0.2.0   2017-04-11 CRAN (R 3.4.0)
    ##  backports    1.1.1   2017-09-25 CRAN (R 3.4.2)
    ##  base       * 3.4.2   2017-10-04 local         
    ##  bindr        0.1     2016-11-13 CRAN (R 3.4.0)
    ##  bindrcpp   * 0.2     2017-06-17 CRAN (R 3.4.0)
    ##  cellranger   1.1.0   2016-07-27 CRAN (R 3.4.0)
    ##  compiler     3.4.2   2017-10-04 local         
    ##  datasets   * 3.4.2   2017-10-04 local         
    ##  devtools     1.13.3  2017-08-02 CRAN (R 3.4.1)
    ##  digest       0.6.12  2017-01-27 CRAN (R 3.4.0)
    ##  dplyr      * 0.7.4   2017-09-28 CRAN (R 3.4.2)
    ##  evaluate     0.10.1  2017-06-24 CRAN (R 3.4.1)
    ##  glue         1.2.0   2017-10-29 CRAN (R 3.4.2)
    ##  graphics   * 3.4.2   2017-10-04 local         
    ##  grDevices  * 3.4.2   2017-10-04 local         
    ##  hms          0.3     2016-11-22 CRAN (R 3.4.0)
    ##  htmltools    0.3.6   2017-04-28 cran (@0.3.6) 
    ##  knitr        1.17    2017-08-10 CRAN (R 3.4.1)
    ##  magrittr   * 1.5     2014-11-22 CRAN (R 3.4.0)
    ##  memoise      1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods    * 3.4.2   2017-10-04 local         
    ##  pkgconfig    2.0.1   2017-03-21 CRAN (R 3.4.0)
    ##  purrr      * 0.2.4   2017-10-18 CRAN (R 3.4.2)
    ##  R6           2.2.2   2017-06-17 CRAN (R 3.4.0)
    ##  Rcpp         0.12.13 2017-09-28 CRAN (R 3.4.2)
    ##  readr      * 1.1.1   2017-05-16 CRAN (R 3.4.0)
    ##  readxl     * 1.0.0   2017-04-18 CRAN (R 3.4.0)
    ##  rematch      1.0.1   2016-04-21 CRAN (R 3.4.0)
    ##  rlang        0.1.4   2017-11-05 CRAN (R 3.4.2)
    ##  rmarkdown    1.6     2017-06-15 CRAN (R 3.4.0)
    ##  rprojroot    1.2     2017-01-16 cran (@1.2)   
    ##  stats      * 3.4.2   2017-10-04 local         
    ##  stringi      1.1.5   2017-04-07 CRAN (R 3.4.0)
    ##  stringr      1.2.0   2017-02-18 CRAN (R 3.4.0)
    ##  tibble       1.3.4   2017-08-22 CRAN (R 3.4.1)
    ##  tidyr      * 0.7.2   2017-10-16 CRAN (R 3.4.2)
    ##  tidyselect   0.2.3   2017-11-06 CRAN (R 3.4.2)
    ##  tools        3.4.2   2017-10-04 local         
    ##  utils      * 3.4.2   2017-10-04 local         
    ##  withr        2.1.0   2017-11-01 CRAN (R 3.4.2)
    ##  yaml         2.1.14  2016-11-12 cran (@2.1.14)
