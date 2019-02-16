Convert raw Excel file into simpler data frames
================
Lucas Nell
16 Feb 2019

This reads the initial Excel file and simplifies it into data frames. It
needs to be sourced every time you run `doc/02-aggregate.Rmd`.

**Load packages:**

``` r
suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(readr)
    library(tidyr)
    library(purrr)
})
```

Function to get integer index from Excel alphabetic index

``` r
.i <- function(xlcol) {
    xlcol <- tolower(xlcol)
    rcol <- sum(sapply(1:nchar(xlcol), 
                       function (i) {
                           ind_i <- which(letters == substr(xlcol, i, i))
                           ind <- length(letters)^(nchar(xlcol) - i) * ind_i
                           return(ind)
                       }))
    return(rcol)
}
.i('Z')
```

    ## [1] 26

# Morphometric data

Read from Excel file.

``` r
xl <- read_excel('data/raw_data.xlsx', col_names = paste0("X__", 1:88), 
                 col_types = rep('text', 88), sheet = 1)
```

Initial manipulation to retrieve desired info and correct ambiguity /
“untidy-ness”.

``` r
# Table of abbreviated species names
spp_df_ <- xl[69:86,2:3] %>% 
    rename(abbrev = X__2, full = X__3) %>% 
    arrange(full)
# Resolving ambiguous and untidy abbreviated names
spp_df_$abbrev[spp_df_$abbrev == 'Ds1'] <- 'Ds'
spp_df_$abbrev[spp_df_$full == 'Eumops glaucinus'] <- 'Eg'
spp_df_$abbrev[spp_df_$full == 'Microtus pennsylvanicus'] <- 'Microtus'
spp_df_$abbrev[spp_df_$full == 'Molossus molossus'] <- 'Mmo'


# Start of final morphometrics data frame, starting with diet, clade, name, and 
# villus heights
initial_df <- xl[c(4:34,36:64),1:6] %>%  # (<-- row 35 is all NAs)
    set_names(c('diet', 'clade', 'id', 'prox', 'mid', 'dist')) %>% 
    mutate_at(vars(prox, mid, dist),
              function(x) as.numeric(ifelse(x == 'None', NA, x))) %>%
    mutate(measure = 'villus height')

# Resolving ambiguous ids
initial_df$id[initial_df$clade == 'Bat' & initial_df$id == 'E'] <- 'Eg'
initial_df$id[initial_df$clade == 'Bat' & grepl('Mm', initial_df$id)] <- paste0('Mmo', 1:3)
initial_df$id[grepl('My', initial_df$id)] <- paste0('My', 1:3)

# Assigning species names
initial_df <- initial_df %>%
    mutate(species = sapply(id, function(.id) {
        spp_df_$full[spp_df_$abbrev == gsub('[0-9]', '', .id)]
    })) %>% 
    select(diet, clade, species, id, measure, everything())
```

I’m next going to make a new data frame of column names along with their
column location in the Excel file and the number of additional columns
to take (after the one specified). It’s assumed that all these are
located in rows `c(4:34,36:64)`, like the villus height data. Also, `n`
should only be 2 or 0, for measures with and without separate proximal,
middle, and distal measurements, respectively. Lastly, although I am
extracting columns from the Excel file, I want the new csv file to be in
“long” format, so I’m binding them by rows.

``` r
new_cols <- read_csv('name,col1,n
"villus width",M,2
"crypt width",U,2
"enterocyte density",AC,2
sef,AM,2
"intestinal diameter",AU,2
"intestinal length",BB,0
mass,BH,0
"vill surface area",BV,0
"enterocyte width",CF,2
nsa,BE,0')

new_cols <- new_cols %>% 
    mutate(col1 = as.integer(sapply(col1, .i) - 1))


# Function to extract info from the xl data frame based on a row from new_cols
xl_extr <- function(.x) {
    .df <- xl[c(4:34,36:64), .x$col1:(.x$col1+.x$n)]
    if (.x$n == 0) {
        lhs <- sprintf('X__%i', .x$col1 + 1:2)
        .df <- .df %>%
            mutate(!!lhs[1] := NA_real_, !!lhs[2] := NA_real_)
    }
    
    .df <- .df %>%
        set_names(c('prox', 'mid', 'dist')) %>% 
        mutate_all(function(x) as.numeric(ifelse(x == 'None', NA, x))) %>% 
        mutate(measure = .x$name)
    # Now adding columns common to all dataframes
    out_df <- bind_cols(initial_df %>% select(diet, clade, species, id), .df) %>% 
        select(diet, clade, species, id, measure, everything())
    
    return(out_df)
}

morph_df <- new_cols %>%
    split(.$name) %>% 
    map(xl_extr) %>%
    bind_rows() %>%
    bind_rows(initial_df) %>%
    arrange(measure, species, id)


# Some fixes to errors in the initial Excel file
morph_df <- morph_df %>%
    # These values were input with an extra zero in the Excel file
    mutate(dist = ifelse(species == 'Microtus pennsylvanicus' &
                              measure == 'enterocyte width', 
                         dist * 10, dist),
           mid = ifelse(species == 'Microtus pennsylvanicus' &
                             measure == 'enterocyte width', 
                        mid * 10, mid)) %>% 
    # These species' diets weren't included
    mutate(diet = ifelse(species == 'Microtus pennsylvanicus', 'Herbivorous',
                         ifelse(species == 'Eptesicus fuscus', 'Protein', diet))) %>% 
    # Oligoryzomys seems to be the more standard spelling
    mutate(species = ifelse(species == 'Olygoryzomys nigripes', 
                        'Oligoryzomys nigripes', species))

# Number of individuals
N <- morph_df$id %>% unique %>% length

# Measures with no position (i.e., NA in pos column instead of prox, mid, or dist)
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
                          gsub(' ', '_', measure))) %>%
    # Changing back from tall to wide format
    spread(measure, value)

# These objects are no longer necessary
rm(no_pos, N)
```

# Clearance data

For “L-arabinose clearance (μl min^-1)” vs SEF (Figure 7A)

Necessary functions:

``` r
# Replacing abbreviated names for full species namees
full_spp <- function(species) {
    spp <- unique(morph_df$species)
    spp_split <- sapply(strsplit(spp, '\\s+'), function(x) tolower(x[2]))
    species_split <- sapply(strsplit(species, '\\s+'), function(x) tolower(x[2]))
    species_out <- sapply(species_split, function(s) spp[spp_split == s][1], 
                          USE.NAMES = FALSE)
    # This is the brown rat.
    species_out[is.na(species_out) & species == 'R. norvegicus'] <- 'Rattus norvegicus'
    return(species_out)
}
# Abbreviate ids once I have the full species info
abbrev_id <- function(ids) {
    abbrev_ids <- sapply(strsplit(ids, '\\s+'), 
                        function(x) {
                            paste(substr(x[1],1,1), 
                                  substr(x[2], 1, 4),
                                  x[3],
                                  sep = "_")
                        })
    return(abbrev_ids)
}
# Next two functions: Find clade and diet from a full species name
find_clade <- function(species) {
    taxa <- sapply(species, function(s) morph_df$clade[morph_df$species == s][1])
    # In case you didn't know, rats are rodents.
    taxa[is.na(taxa) & species == 'Rattus norvegicus'] <- 'Rodent'
    return(taxa)
}
find_diet <- function(species) {
    diets <- sapply(species, function(s) morph_df$diet[morph_df$species == s][1])
    # In case you didn't know, rats are rodents.
    diets[is.na(diets) & species == 'Rattus norvegicus'] <- 'Omnivorous'
    return(diets)
}
# Boolean vector for whether items in (character) vector are NAs or could be coerced 
# into a numeric
is_num <- function(x) {
    sapply(x, function(.x) is.na(.x) | !is.na(suppressWarnings(as.numeric(.x))))
}
```

Reading the
data.

``` r
xl_sef <- read_excel('data/raw_data.xlsx', sheet = 2, range = "B20:E138",
                     col_types = rep('text', 4),
                     col_names = c('id', 'prox', 'mid', 'dist'))


xl_clear <- read_excel('data/raw_data.xlsx', col_names = c('id', 'clear'), 
                       col_types = rep('text', 2), sheet = 2, range = "I20:J140")
```

Manipulating for SEF and clearance, then combining them into one data
frame.

``` r
sef_df <- xl_sef %>%
    filter(!is.na(id) | !is.na(prox) | !is.na(mid) | !is.na(dist),
           is_num(prox), is_num(mid), is_num(dist), 
           id != "SEF") %>% 
    mutate_at(vars(prox, mid, dist), as.numeric) %>% 
    mutate(
        # this one was input incorrectly
        id = gsub("perspicillatac", "perspicillata", id),
        species = sapply(strsplit(id, "\\s+"), 
                         function(x) paste(head(x, -1), collapse = " ")),
        species = full_spp(species),
        # now changing ids to abbreviations since I now have full species names
        id = abbrev_id(id),
        # Now for clade, then diet
        diet = find_diet(species),
        clade = find_clade(species)
    ) %>% 
    select(diet, clade, species, id, prox, mid, dist)



clear_df <- xl_clear %>% 
    filter(!is.na(id), !is.na(clear), is_num(clear)) %>% 
    mutate(
        clear = as.numeric(clear),
        # this one was input incorrectly
        id = gsub("perspicillatac", "perspicillata", id),
        species = sapply(strsplit(id, " "), 
                         function(x) paste(x[1:2], collapse = " ")),
        species = full_spp(species),
        # now changing ids to abbreviations since I now have full species names
        id = abbrev_id(id),
        # Now for clade, then diet
        diet = find_diet(species),
        clade = find_clade(species)
        ) %>% 
    select(diet, clade, species, id, clear)


# Combine
clear_df <- bind_rows(sef_df, clear_df) %>% 
    arrange(clade, diet, species, id)
```

# Absorption data

For “Fractional absorption / total intestinal surface (cm^-2)” vs clade
(Figure 7B)

This uses many of the same functions that the clearance data did.

Reading
data:

``` r
xl_abs <- read_excel('data/raw_data.xlsx', sheet = 3, range = "C19:M131",
                     col_types = rep('text', 11),
                     col_names = c('id', 'gavage', 'injection', 
                                   'id2', 'prox', 'mid', 'dist', 'animal_avg', 
                                   'sp_avg', 'nsa', 'mass'))
```

Manipulating for the two main sections of `xl_abs`, then combining them
into one data frame.

``` r
abs_df <- xl_abs %>%
    select(id, gavage, injection) %>% 
    filter(!is.na(id), !is.na(gavage) | !is.na(injection),
           is_num(gavage), is_num(injection), 
           !grepl(id, pattern = 'average|absorption', ignore.case = TRUE)) %>% 
    mutate_at(vars(gavage, injection), as.numeric) %>% 
    mutate(
        # this one was input incorrectly
        id = gsub("Microtus", "M. pennsylvanicus", id),
        species = sapply(strsplit(id, "\\s+"), 
                         function(x) paste(head(x, -1), collapse = " ")),
        species = full_spp(species),
        # now changing ids to abbreviations since I now have full species names
        id = abbrev_id(id),
        # Now for clade, then diet
        diet = find_diet(species),
        clade = find_clade(species)
    ) %>% 
    select(diet, clade, species, id, gavage, injection)


abs_df2 <- xl_abs %>%
    select(id2, prox, mid, dist, nsa, mass) %>% 
    filter(!is.na(id2), !is.na(prox) , !is.na(mid),
           !is.na(dist), !is.na(nsa), !is.na(mass)) %>% 
    mutate_at(vars(prox, mid, dist, nsa, mass), as.numeric) %>% 
    mutate(
        # this one was input incorrectly
        id2 = gsub("Microtus", "M. pennsylvanicus ", id2),
        species = sapply(strsplit(id2, "\\s+"), 
                         function(x) paste(head(x, -1), collapse = " ")),
        species = full_spp(species),
        # now changing ids to abbreviations since I now have full species names
        id = abbrev_id(id2),
        # Now for clade, then diet
        diet = find_diet(species),
        clade = find_clade(species)
    ) %>% 
    select(diet, clade, species, id, prox, mid, dist, nsa, mass)

absorp_df <- bind_rows(abs_df, abs_df2) %>% 
    arrange(clade, diet, species, id)
```

# Remove excess objects

Only the objects `morph_df`, `clear_df`, and `absorp_df` are needed
downstream.

``` r
rm(list = c(".i", "abbrev_id", "abs_df",
            "abs_df2", "find_clade", "find_diet", "full_spp",
            "initial_df", "is_num", "new_cols", "sef_df",
            "spp_df_", "xl", "xl_abs", "xl_clear",
            "xl_extr", "xl_sef"))
```

# Session info

This outlines the package versions I used for this
    script.

    ## ─ Session info ──────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 3.5.2 (2018-12-20)
    ##  os       macOS Mojave 10.14.3        
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  ctype    en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2019-02-16                  
    ## 
    ## ─ Packages ──────────────────────────────────────────────────────────────
    ##  package     * version date       lib source        
    ##  assertthat    0.2.0   2017-04-11 [1] CRAN (R 3.5.0)
    ##  backports     1.1.3   2018-12-14 [1] CRAN (R 3.5.0)
    ##  callr         3.1.1   2018-12-21 [1] CRAN (R 3.5.0)
    ##  cellranger    1.1.0   2016-07-27 [1] CRAN (R 3.5.0)
    ##  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.0)
    ##  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
    ##  desc          1.2.0   2018-05-01 [1] CRAN (R 3.5.0)
    ##  devtools      2.0.1   2018-10-26 [1] CRAN (R 3.5.1)
    ##  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.0)
    ##  dplyr       * 0.8.0.1 2019-02-15 [1] CRAN (R 3.5.2)
    ##  evaluate      0.13    2019-02-12 [1] CRAN (R 3.5.2)
    ##  fs            1.2.6   2018-08-23 [1] CRAN (R 3.5.0)
    ##  glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.0)
    ##  highr         0.7     2018-06-09 [1] CRAN (R 3.5.0)
    ##  hms           0.4.2   2018-03-10 [1] CRAN (R 3.5.0)
    ##  htmltools     0.3.6   2017-04-28 [1] CRAN (R 3.5.0)
    ##  knitr         1.21    2018-12-10 [1] CRAN (R 3.5.2)
    ##  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
    ##  memoise       1.1.0   2017-04-21 [1] CRAN (R 3.5.0)
    ##  pillar        1.3.1   2018-12-15 [1] CRAN (R 3.5.0)
    ##  pkgbuild      1.0.2   2018-10-16 [1] CRAN (R 3.5.0)
    ##  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.0)
    ##  pkgload       1.0.2   2018-10-29 [1] CRAN (R 3.5.0)
    ##  prettyunits   1.0.2   2015-07-13 [1] CRAN (R 3.5.0)
    ##  processx      3.2.1   2018-12-05 [1] CRAN (R 3.5.0)
    ##  ps            1.3.0   2018-12-21 [1] CRAN (R 3.5.0)
    ##  purrr       * 0.3.0   2019-01-27 [1] CRAN (R 3.5.2)
    ##  R6            2.4.0   2019-02-14 [1] CRAN (R 3.5.2)
    ##  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.0)
    ##  readr       * 1.3.1   2018-12-21 [1] CRAN (R 3.5.0)
    ##  readxl      * 1.3.0   2019-02-15 [1] CRAN (R 3.5.2)
    ##  rematch       1.0.1   2016-04-21 [1] CRAN (R 3.5.0)
    ##  remotes       2.0.2   2018-10-30 [1] CRAN (R 3.5.0)
    ##  rlang         0.3.1   2019-01-08 [1] CRAN (R 3.5.2)
    ##  rmarkdown     1.11    2018-12-08 [1] CRAN (R 3.5.0)
    ##  rprojroot     1.3-2   2018-01-03 [1] CRAN (R 3.5.0)
    ##  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.0)
    ##  stringi       1.3.1   2019-02-13 [1] CRAN (R 3.5.2)
    ##  stringr       1.4.0   2019-02-10 [1] CRAN (R 3.5.2)
    ##  testthat      2.0.1   2018-10-13 [1] CRAN (R 3.5.0)
    ##  tibble        2.0.1   2019-01-12 [1] CRAN (R 3.5.2)
    ##  tidyr       * 0.8.2   2018-10-28 [1] CRAN (R 3.5.0)
    ##  tidyselect    0.2.5   2018-10-11 [1] CRAN (R 3.5.0)
    ##  usethis       1.4.0   2018-08-14 [1] CRAN (R 3.5.0)
    ##  withr         2.1.2   2018-03-15 [1] CRAN (R 3.5.0)
    ##  xfun          0.4     2018-10-23 [1] CRAN (R 3.5.0)
    ##  yaml          2.2.0   2018-07-25 [1] CRAN (R 3.5.0)
    ## 
    ## [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library
