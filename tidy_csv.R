# 
# This cleans csv files for use.
# 

# I'm now making sure all necessary packages are loaded:
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(magrittr)
})


morph_df <- read_csv('./data/clean_data.csv', col_types = 'cccccddd') %>%
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

# This object is no longer necessary
rm(no_pos)

# These values were input with an extra zero in the Excel file
morph_df <- morph_df %>% 
    mutate(value = ifelse(species == 'Microtus pennsylvanicus' &
                          measure == 'enterocyte_diameter' &
                          pos %in% c('dist', 'med'), value * 10, value))




# --------------------
# Clearance and absorption data
# --------------------

clear_df <- read_csv('data/clean_clearance_data.csv', col_types = 'cccdd') %>%
    # In the plot, they lumped herbivores and omnivores together as "carb eater <taxon>"
    mutate(diet = ifelse(diet == "Protein", diet, "Carb")) %>% 
    as.data.frame %>% 
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
           diet = factor(diet, levels = c("Carb", "Protein")))
row.names(clear_df) <- paste(clear_df$species)

absorp_df <- read_csv('data/clean_absorption_data.csv', col_types = 'ccd') %>%
    as.data.frame %>% 
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')))
row.names(absorp_df) <- paste(absorp_df$species)
