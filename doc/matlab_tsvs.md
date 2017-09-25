Create MATLAB TSVs
================
Lucas Nell
25 Sep 2017

-   [Instructions for MATLAB script](#instructions-for-matlab-script)
-   [By species](#by-species)
-   [Absorption](#absorption)
-   [By species and segment](#by-species-and-segment)

This script creates TSVs ("tab-delimited values") for use in the MATLAB code.

**Load libraries:**

``` r
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(phylolm)
    library(ape)
    library(ggplot2)
})
```

The `R/get_data.R` file provides functions to retrieve morphometric, clearance, and absorption data. See `tidy_csvs.md` for more info.

``` r
invisible(source("R/get_data.R"))
```

Instructions for MATLAB script
==============================

All input files should be tab-delimited text files.

-   Data file: 1 row for each species
-   Measurement error file: standard error for each cell in data file
-   MATLAB file to run: `physigv2.m`
-   Run it, and dialogs will show up
-   `taxon` measurement error = 0
-   Use a variance-covariance matrix (using `ape::vcov`) for the phylogenetic tree file
-   You can use likelihoods for p-values: run it with and without, `2*logLik` is chi-square
-   Run it with ML, not REML

By species
==========

Retrieve data, change factors to binary variable(s), and remove unnecessary parameters:

``` r
spp_mean <- get_df('spp', 'mean') %>% 
    as_data_frame %>% 
    mutate(taxonBat = ifelse(taxon == 'Bat', 1L, 0L), 
           dietOmnivorous = ifelse(diet == 'Omnivorous', 1L, 0L),
           dietProtein = ifelse(diet == 'Protein', 1L, 0L)) %>% 
    select(species, taxonBat, dietOmnivorous, dietProtein, log_mass, everything(), 
           -diet, -taxon) %>% 
    arrange(species) %>% 
    select(-species)
spp_se <- get_df('spp', 'se') %>%
    as_data_frame %>% 
    # no measurement error in diet or taxon:
    mutate(taxonBat = 0L, dietOmnivorous = 0L, dietProtein = 0L) %>%
    select(species, taxonBat, dietOmnivorous, dietProtein, log_mass, everything(), 
           -diet, -taxon) %>% 
    arrange(species) %>% 
    select(-species)
```

Write files:

``` r
write_tsv(spp_mean, 'matlab/data/spp_mean.txt')
write_tsv(spp_se, 'matlab/data/spp_se.txt')
```

Phylogenetic tree for these data
--------------------------------

``` r
tr <- get_tr('spp')
tr_df <- vcv(tr)[order(tr$tip.label),order(tr$tip.label)] %>%
    as_data_frame
write_tsv(tr_df, 'matlab/data/spp_vcv.txt', col_names = FALSE)
```

Absorption
==========

Retrieve data, change factors to binary variable(s), and remove unnecessary parameters:

``` r
absorp_mean <- get_df('absorp', 'mean') %>% 
    as_data_frame %>% 
    mutate(taxonBat = ifelse(taxon == 'Bat', 1L, 0L)) %>% 
    select(species, taxonBat, everything(), -taxon) %>% 
    arrange(species) %>% 
    select(-species)
absorp_se <- get_df('absorp', 'se') %>%
    as_data_frame %>% 
    mutate(taxonBat = 0L) %>%  # <- no measurement error in taxon
    select(species, taxonBat, everything(), -taxon) %>% 
    arrange(species) %>% 
    select(-species)
```

Write files:

``` r
write_tsv(absorp_mean,  'matlab/data/absorp_mean.txt')
write_tsv(absorp_se, 'matlab/data/absorp_se.txt')
```

Phylogenetic tree for these data
--------------------------------

``` r
tr <- get_tr('absorp')
tr_df <- vcv(tr)[order(tr$tip.label),order(tr$tip.label)] %>%
    as_data_frame
write_tsv(tr_df, 'matlab/data/absorp_vcv.txt', col_names = FALSE)
```

By species and segment
======================

Retrieve data, change factors to binary variable(s), remove unnecessary parameters, and write output for each segment type:

``` r
for (s in c('prox', 'med', 'dist')) {
    pos_mean <- get_df('pos', 'mean', s) %>%
        as_data_frame %>% 
        mutate(taxonBat = ifelse(taxon == 'Bat', 1L, 0L)) %>% 
        select(species, taxonBat, log_mass, everything(), -diet, -taxon) %>% 
        arrange(species) %>% 
        select(-species)
    pos_se <- get_df('pos', 'se', s) %>%
        as_data_frame %>% 
        mutate(taxonBat = 0L) %>%  # <-- no measurement error in taxon
        select(species, taxonBat, log_mass, everything(), -diet, -taxon) %>% 
        arrange(species) %>% 
        select(-species)
    write_tsv(pos_mean, paste0('matlab/data/pos_', s, '_mean.txt'))
    write_tsv(pos_se, paste0('matlab/data/pos_', s, '_se.txt'))
}; rm(s, pos_mean, pos_se)
```

Phylogenetic tree for these data
--------------------------------

The same tree is used for these as for analyses by species only, so the file `matlab/data/spp_vcv.txt` can be used.
