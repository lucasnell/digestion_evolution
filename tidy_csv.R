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


# ======================================================================================
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==

# Full morphometric data frame

# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==
# ======================================================================================


morph_df <- read_csv('data/clean_morph_data.csv', col_types = 'cccccddd') %>%
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



# ======================================================================================
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==

# Morphometric data aggregated by species

# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==
# ======================================================================================


spp_measures <- c('mass', 'intestinal_length', 'nsa', 'villa_surface_area',
                  'enterocyte_density', 'sef')

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
           total_enterocytes = enterocyte_density * nsa,
           log_total_enterocytes = log(enterocyte_density * nsa),
           total_surface = nsa * sef,
           log_mass = log(mass)) %>% 
    select(taxon, diet, species, pos,
           int_length_mass, nsa_mass, vill_area_mass, 
           total_enterocytes, log_total_enterocytes, total_surface, sef, 
           mass, log_mass) %>% 
    # Taking mean by position (dist, med, prox, or NA)
    group_by(taxon, diet, species, pos) %>% 
    summarize_all(mean, na.rm = TRUE) %>% 
    ungroup %>% 
    # Now taking mean by species
    group_by(taxon, diet, species) %>% 
    summarize_at(.vars = vars(int_length_mass, nsa_mass, vill_area_mass, 
                              total_enterocytes, log_total_enterocytes, total_surface, 
                              sef, mass, log_mass), 
                 mean, na.rm = TRUE) %>% 
    ungroup %>% 
    # Converting taxon and diet to factors
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
           diet = factor(diet, levels = c("Herbivorous", "Omnivorous", "Protein"))) %>% 
    arrange(taxon, diet, species) %>% 
    # To change row names, it can't be a tibble, so I'm reverting back to normal
    # data frame
    as.data.frame

# phylolm requires that the rownames match the species names
rownames(spp_df) <- spp_df$species


# Removing NA from data frame for diet
diet_df <- spp_df %>% filter(!is.na(diet))
rownames(diet_df) <- diet_df$species




# ======================================================================================
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==

# Morphometric data aggregated by species AND position

# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==
# ======================================================================================


pos_measures <- c('mass', 'intestinal_diameter', 'villus_height',  'villus_width',
                  'crypt_width', 'sef', 'enterocyte_diameter', 'enterocyte_density')

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
    # Grouping by, then taking mean of all measurement columns and transformed-
    # measurement columns
    group_by(taxon, diet, species, pos) %>% 
    summarize_at(.vars = c(pos_measures, paste0('log_', pos_measures)), mean, 
                 na.rm = TRUE) %>% 
    ungroup %>% 
    # Converting taxon and diet to factors
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
           diet = factor(diet, levels = c("Herbivorous", "Omnivorous", "Protein"))) %>% 
    arrange(taxon, diet, species, pos)


for (p in unique(pos_df$pos)) {
    out_df <- pos_df %>% 
        filter(pos == p) %>% 
        select(-pos) %>% 
        as.data.frame
    rownames(out_df) <- out_df$species
    assign(paste0(p, '_df'), out_df)
}; rm(p, out_df)





# ======================================================================================
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==

# Clearance and absorption data

# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==
# ======================================================================================

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
