Clean CSVs into simpler data frames
================
Lucas Nell
11 Sep 2017

-   [Summary of output](#summary-of-output)
-   [Getting started](#getting-started)
-   [Full morphometric data frame](#full-morphometric-data-frame)
-   [Morphometric data aggregated by species](#morphometric-data-aggregated-by-species)
-   [Morphometric data by species AND position](#morphometric-data-by-species-and-position)
-   [Clearance data](#clearance-data)
-   [Absorption data](#absorption-data)
-   [Saving these data frames](#saving-these-data-frames)
-   [Session info](#session-info)

This reads the cleaned CSV files with data for individual animals and simplifies them into CSV files with summary statistics by species and/or intestinal segment. It needs to be run only once, and not at all if you have the csv files already.

Summary of output
=================

I created two data sets from morphometric data in `output/clean_morph_data.csv`:

1.  Measurements separated by species (found in `output/spp_df.csv`)
2.  Measurements separated by species and intestinal segment (`output/pos_df.csv`)

I also created one data set each from `output/clean_clearance_data.csv` and `output/clean_absorption_data.csv`. These csvs are found in `output/clear_df.csv` and `output/absorp_df.csv`, respectively.

Each data set is associated with a set of analyses and only needs certain measurements, so only those are included.

The csv files are used in the function `get_df` in the `R/get_data.R` file to retrieve a data frame for a given analysis set. The analysis sets are as follows:

1.  `'spp'`: Measurements separated by species.
2.  `'diet'`: Measurements by species and with non-`NA` diet data.
3.  `'pos'`: Measurements by species and position. (You have to also provide the position for this analysis set.)
4.  `'clear'`: Clearance data by species. (Uses a different set of individuals entirely from the 1–3.)
5.  `'absorp'`: Absorption data by species. (Uses a different set of individuals entirely from the 1–4.)

For the third analysis set, I need to do the analyses separately for each position because modelling within-species and within-individual variance due to position rather than measurement error or process error would be difficult and not likely possible with this small dataset.

Getting started
===============

**Load packages:**

``` r
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(magrittr)
})
```

Functions for standard error calculations

``` r
# Calculating standard error
# (The ... is added for compatibility with using alongside mean.)
se <- function(.x, ...) {
    .z <- .x[!is.na(.x)]
    return(sd(.z) / sqrt(length(.z)))
}
# Replacing NAs in a vector of standard errors with the average standard error in 
# that vector.
na_se <- function(.x) {
    .z <- mean(.x, na.rm = TRUE)
    .x[is.na(.x)] <- .z
    return(.x)
}
```

Full morphometric data frame
============================

This data frame of morphometric measurements will be useful for the next two sections.

``` r
morph_df <- read_csv('output/clean_morph_data.csv', col_types = 'cccccddd') %>%
    # Oligoryzomys seems to be the more standard spelling
    mutate(species = ifelse(species == 'Olygoryzomys nigripes', 
                            'Oligoryzomys nigripes', species))

# Number of individuals
N <- morph_df$id %>% unique %>% length

# Measures with no position (i.e., NA in pos column instead of prox, med, or dist)
no_pos <- morph_df %>% 
    filter(is.na(dist)) %>% 
    group_by(measure) %>% 
    summarize(total = n()) %>% 
    filter(total == N) %>% 
    select(measure) %>% 
    unlist %>% 
    paste

# Gathering into 'tall' format, fixing position column, removing spaces from 
# measure column, and changing the "enterocyte_width" measure to "enterocyte_diameter".
morph_df <- morph_df %>% 
    gather(pos, value, prox:dist, na.rm = TRUE) %>% 
    mutate(pos = ifelse(measure %in% no_pos, NA, pos),
           measure = gsub('enterocyte_width', 'enterocyte_diameter', 
                          gsub(' ', '_', measure)))

# These objects are no longer necessary
rm(no_pos, N)
```

Morphometric data aggregated by species
=======================================

For the following measures, we need a single morphometric value per species:

-   Intestinal length / body mass^0.4
-   NSA / body mass^0.75
-   Villus surface area / body mass^0.75
-   Total number of enterocytes (log-transformed; log body mass as covariate)
-   Calculated as such: `NSA * mean(<enterocyte density among segments>)`
-   Fractional absorption / (total intestinal surface / mass^0.75)
-   total intestinal surface = `NSA * SEF`

We need the following columns from `morph_df` to compute these values:

``` r
spp_measures <- c('mass', 'intestinal_length', 'nsa', 'villa_surface_area',
                  'enterocyte_density', 'sef')
```

Next I manipulated `morph_df` as such to get means and standard errors for each parameter that I will be analyzing.

``` r
spp_df <- morph_df %>%
    # Changing from tall to wide format
    spread(measure, value) %>% 
    # Selecting measurement columns, plus the identifying columns
    select_(.dots = append(list('diet', 'taxon', 'species', 'id', 'pos'), 
                           spp_measures)) %>% 
    # Removing all rows with all NAs in measures columns
    filter(Reduce(`+`, lapply(.[,spp_measures], is.na)) < length(spp_measures)) %>% 
    # Add nsa to all positions' estimates (for total_enterocytes and total_surface below)
    group_by(taxon, diet, species, id) %>% 
    mutate(nsa = ifelse(is.na(nsa), nsa[!is.na(nsa)], nsa)) %>% 
    ungroup %>% 
    # Doing the calculations now, before taking any means
    mutate(int_length_mass = intestinal_length / mass^0.4,
           nsa_mass = nsa / mass^0.75,
           vill_area_mass = villa_surface_area / mass^0.75,
           log_total_enterocytes = log(enterocyte_density * nsa),
           log_mass = log(mass)) %>% 
    select(taxon, diet, species, id,
           sef, int_length_mass, nsa_mass, vill_area_mass, 
           log_total_enterocytes, log_mass) %>% 
    # Taking mean by individual (i.e., across segments)
    group_by(taxon, diet, species, id) %>% 
    summarize_all(mean, na.rm = TRUE) %>% 
    ungroup %>% 
    # Now taking mean and SE by species
    select(-id) %>% 
    group_by(taxon, diet, species) %>% 
    summarize_all(funs(mean, se), na.rm = TRUE) %>% 
    ungroup %>% 
    # Adjusting species with NA for their SE (i.e., those spp with only 1 individual)
    mutate_at(vars(ends_with('_se')), na_se) %>%
    arrange(taxon, diet, species)
```

Lastly I add a column named `stat` to distinguish between mean and SE values, making filtering for a given statistic easier.

``` r
spp_df <- spp_df %>% 
    gather('measure', 'value', -taxon, -diet, -species) %>% 
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

-   Intestinal diameter (log-transformed, body mass as covariate)
-   Villus height (body mass as covariate)
-   Villus width
-   Crypt width
-   Surface enlargement factor (SEF)
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
    # Changing from tall to wide format
    spread(measure, value) %>% 
    select_(.dots = append(list('taxon', 'diet', 'species', 'pos', 'id'), 
                           as.list(pos_measures))) %>% 
    # Add mass to all positions' estimates
    group_by(taxon, diet, species, id) %>% 
    mutate(mass = ifelse(is.na(mass), mass[!is.na(mass)], mass)) %>% 
    ungroup %>% 
    # Now removing rows with pos == NA bc they don't have the other measurements
    filter(!is.na(pos)) %>% 
    # Taking the log now, before taking any means
    mutate_(.dots = setNames(as.list(sprintf('%s(%s)', 'log', pos_measures)), 
                             paste0('log_', pos_measures))) %>% 
    # Grouping by, then taking mean and SE of all measurement columns and transformed-
    # measurement columns
    group_by(taxon, diet, species, pos) %>% 
    summarize_at(.vars = c(pos_measures, paste0('log_', pos_measures)), funs(mean, se), 
                 na.rm = TRUE) %>% 
    ungroup %>% 
    # Now only outputting columns that are necessary
    select_(
        .dots = as.list(c('taxon', 'diet', 'species', 'pos', 
          paste0(rep(c('log_intestinal_diameter', 'villus_height', 'villus_width',
                       'crypt_width', 'sef', 'enterocyte_diameter',
                       'log_enterocyte_density', 'log_mass'), each = 2), 
                 c('_mean', '_se'))))
    ) %>% 
    arrange(taxon, diet, species, pos)
```

Lastly I add a column named `stat` to distinguish between mean and SE values, making filtering for a given statistic easier.

``` r
pos_df <- pos_df %>% 
    gather('measure', 'value', -taxon, -diet, -species, -pos) %>% 
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
clear_df <- read_csv('output/clean_clearance_data.csv', col_types = 'ccccdddd') %>%
    mutate(
        # They lumped herbivores and omnivores together as "carb eater <taxon>"
        diet = ifelse(diet == "Protein", diet, "Carb"),
        # Averaging SEF by individual, on log scale
        log_sef = (log(prox) + log(med) + log(dist)) / 3,
        # Taking log of clearance before any means are calculated
        # Some clearances were negative, so there would be NaNs produced
        # The below two lines are to avoid the warning message
        log_clear = ifelse(clear < 0, NA, clear),
        log_clear = log(log_clear)
    ) %>% 
    # These are no longer necessary
    select(-prox, -med, -dist, -clear) %>% 
    group_by(diet, taxon, species) %>% 
    summarize_at(vars(log_sef, log_clear), funs(mean, se), na.rm = TRUE) %>% 
    ungroup
```

As above, I add a column named `stat` to distinguish between mean and SE values, making filtering for a given statistic easier.

``` r
clear_df <- clear_df %>% 
    gather('measure', 'value', -taxon, -diet, -species) %>% 
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

`(gavage / injection) / { (nsa * sef) / (mass^0.75) }`

This can be manipulated to the following:

`gavage * (1 / injection) * { (mass^0.75) / (nsa * sef) }`

For three species, the individuals used for measurements for this calculation are split into three groups, corresponding to the three portions of the equation above. One set of individuals were used for `gavage`, another for `injection`, and a third for `nsa`, `sef`, and `mass`.

These are the above-mentioned three species:

``` r
sep_absorps <- c('Myotis lucifugus', 'Tadarida brasiliensis', 'Akodon montensis')
```

For the rest of the species, `gavage` and `injection` were measured together, but `nsa`, `sef`, and `mass` were measured in a different set of individuals.

This is important because `mean(X/Y) != mean(X) / mean(Y)` and `var(X/Y) != var(X) / var(Y)`. But because `mean(X*Y) = mean(X) * mean(Y)` and `var(XY) = mean(X)^2 * var(Y) + mean(Y)^2 * var(X) + var(X) * var(Y)`, I can inverse things before taking means or variances to allow me to combine these values.

So before taking means or variances, I have inversed `gavage` and `(nsa * sef) / (mass^0.75)`, thus allowing me to use the equations above. For the three species with `gavage` and `injection` measured separately, I had to combine variances for `gavage` and `injection` first, then combine the variance of `gavage * injection` with `(mass^0.75) / (nsa * sef)`.

The term `(mass^0.75) / (nsa * sef)` is called `mns` below.

``` r
absorp_df <- read_csv('output/clean_absorption_data.csv', col_types = 'ccccddddddd') %>%
    mutate(
        # Averaging SEF by individual (i.e., across segments)
        sef = (prox + med + dist) / 3,
        mns = (mass^0.75) / (nsa * sef)
    ) %>% 
    group_by(diet, taxon, species) %>% 
    summarize(
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
        
        # "gi_" stands for gavage and injection
        gi_v = gavage^2 * inv_injection_v + inv_injection^2 * gavage_v +
            gavage_v * inv_injection_v,
        gi = gavage * inv_injection,
        
        mns_v = var(mns, na.rm = TRUE),
        mns = mean(mns, na.rm = TRUE),
        
        # I added the sqrt bc I want the output to be the standard deviation (the closest
        # to standard error I've figure out to get)
        absorp_se = sqrt(mns^2 * gi_v + gi^2 * mns_v + mns_v * gi_v),
        absorp_mean = gavage * inv_injection * mns
    ) %>% 
    ungroup %>% 
    select(taxon, species, absorp_mean, absorp_se)

rm(sep_absorps)
```

As above, I add a column named `stat` to distinguish between mean and SE values, making filtering for a given statistic easier.

``` r
absorp_df <- absorp_df %>% 
    gather('measure', 'value', -taxon, -species) %>% 
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
write_csv(pos_df, 'output/pos_df.csv')
write_csv(spp_df, 'output/spp_df.csv')
write_csv(clear_df, 'output/clear_df.csv')
write_csv(absorp_df, 'output/absorp_df.csv')
```

Session info
============

This outlines the package versions I used for this script.

    ## Session info -------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.4.1 (2017-06-30)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2017-09-11

    ## Packages -----------------------------------------------------------------

    ##  package    * version date       source        
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
    ##  hms          0.3     2016-11-22 CRAN (R 3.4.0)
    ##  htmltools    0.3.6   2017-04-28 cran (@0.3.6) 
    ##  knitr        1.17    2017-08-10 CRAN (R 3.4.1)
    ##  magrittr   * 1.5     2014-11-22 CRAN (R 3.4.0)
    ##  memoise      1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods    * 3.4.1   2017-07-07 local         
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
    ##  tidyselect   0.2.0   2017-08-30 CRAN (R 3.4.1)
    ##  tools        3.4.1   2017-07-07 local         
    ##  utils      * 3.4.1   2017-07-07 local         
    ##  withr        2.0.0   2017-07-28 CRAN (R 3.4.1)
    ##  yaml         2.1.14  2016-11-12 cran (@2.1.14)
