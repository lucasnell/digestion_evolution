# 
# This is the same as tidy_csv.R, except that it returns data frames of standard errors.
# The data frames are kept as tibbles with no row names because they will not
# be used directly in phylolm.
# 

# I'm now making sure all necessary packages are loaded:
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(magrittr)
})


# Calculating standard error
se <- function(.x) {
    .z <- .x[!is.na(.x)]
    return(sd(.z) / sqrt(length(.z)))
}
# Replacing NAs in a vector of standard errors with the average standard error in 
# that vector
na_se <- function(.x) {
    .z <- mean(.x, na.rm = TRUE)
    .x[is.na(.x)] <- .z
    return(.x)
}

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

# These objects are no longer necessary
rm(no_pos, N)



# ======================================================================================
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==

# Morphometric data aggregated by species

# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==
# ======================================================================================


spp_measures <- c('mass', 'intestinal_length', 'nsa', 'villa_surface_area',
                  'enterocyte_density', 'sef')

spp_se_df <- morph_df %>%
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
    # Doing the calculations now, before taking any means / standard errors
    mutate(int_length_mass = intestinal_length / mass^0.4,
           nsa_mass = nsa / mass^0.75,
           vill_area_mass = villa_surface_area / mass^0.75,
           total_enterocytes = enterocyte_density * nsa,
           log_total_enterocytes = log(enterocyte_density * nsa),
           total_surface = nsa * sef,
           log_mass = log(mass)) %>% 
    select(taxon, diet, species, id,
           int_length_mass, nsa_mass, vill_area_mass, 
           total_enterocytes, log_total_enterocytes, total_surface, sef, 
           mass, log_mass) %>% 
    # Taking mean by individual (i.e., across segments)
    group_by(taxon, diet, species, id) %>% 
    summarize_all(mean, na.rm = TRUE) %>% 
    ungroup %>% 
    # Now taking standard error by species
    group_by(taxon, diet, species) %>% 
    summarize_at(.vars = vars(int_length_mass, nsa_mass, vill_area_mass, 
                              total_enterocytes, log_total_enterocytes, total_surface, 
                              sef, mass, log_mass), 
                 se) %>% 
    ungroup %>% 
    # Adjusting species with NA for their SE (i.e., those spp with only 1 individual)
    mutate_at(.vars = vars(int_length_mass, nsa_mass, vill_area_mass, 
                           total_enterocytes, log_total_enterocytes, total_surface, 
                           sef, mass, log_mass), 
              na_se) %>% 
    # Converting taxon and diet to factors
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
           diet = factor(diet, levels = c("Herbivorous", "Omnivorous", "Protein"))) %>% 
    arrange(taxon, diet, species)

# Removing NA from data frame for diet
diet_se_df <- spp_se_df %>% filter(!is.na(diet))





# ======================================================================================
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==

# Morphometric data aggregated by species AND position

# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==
# ======================================================================================


pos_measures <- c('mass', 'intestinal_diameter', 'villus_height',  'villus_width',
                  'crypt_width', 'sef', 'enterocyte_diameter', 'enterocyte_density')

pos_se_df <- morph_df %>%
    # Changing from tall to wide format
    spread(measure, value) %>% 
    select_(.dots = c(list('taxon', 'diet', 'species', 'pos', 'id'), 
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
    # Grouping by, then taking standard error of all measurement columns and transformed-
    # measurement columns
    group_by(taxon, diet, species, pos) %>% 
    summarize_at(.vars = c(pos_measures, paste0('log_', pos_measures)), se) %>% 
    ungroup %>% 
    # Adjusting species with NA for their SE (i.e., those spp with only 1 individual)
    mutate_at(.vars = c(pos_measures, paste0('log_', pos_measures)), na_se) %>% 
    ungroup %>% 
    # Converting taxon and diet to factors
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
           diet = factor(diet, levels = c("Herbivorous", "Omnivorous", "Protein"))) %>% 
    arrange(taxon, diet, species, pos)



for (p in unique(pos_se_df$pos)) {
    out_se_df <- pos_se_df %>% 
        filter(pos == p) %>% 
        select(-pos)
    assign(paste0(p, '_se_df'), out_se_df)
}; rm(p, out_se_df)
rm(pos_se_df)




# ======================================================================================
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==

# Clearance data

# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==
# ======================================================================================

# To remove warning about logging negative clearance values 
wog <- getOption('warn')
options(warn = -1)

clear_se_df <- read_csv('data/clean_clearance_data.csv', col_types = 'ccccdddd') %>%
    mutate(
        # They lumped herbivores and omnivores together as "carb eater <taxon>"
        diet = ifelse(diet == "Protein", diet, "Carb"),
        # Averaging SEF by individual, both on "identity" and log scale
        sef = (prox + med + dist) / 3,
        log_sef = (log(prox) + log(med) + log(dist)) / 3,
        # Taking log of clearance before any means are calculated
        # Some clearances were negative, so there will be NaNs produced
        log_clear = log(clear)
    ) %>% 
    # SE by species
    group_by(diet, taxon, species) %>% 
    summarize_at(vars(sef, log_sef, clear, log_clear), se) %>% 
    ungroup %>%
    # This dataset has no species with only 1 individual, so no need to replace NAs
    arrange(diet, taxon, species) %>% 
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
           diet = factor(diet, levels = c('Carb', 'Protein')))

# Setting warning setting back to what it was
options(warn = wog)
rm(wog)



# ======================================================================================
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==

# Absorption data

# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==
# ======================================================================================


# These species had gavage and injection done separately in different individuals
sep_absorps <- c('Myotis lucifugus', 'Tadarida brasiliensis', 'Akodon montensis')



absorp_se_df <- read_csv('data/clean_absorption_data.csv', col_types = 'ccccddddddd') %>%
    mutate(
        # Averaging SEF by individual
        sef = (prox + med + dist) / 3,
        # I'm inversing this parameter because...
        # E(X*Y) = E(X) * E(Y), but E(X/Y) != E(X) / E(Y)
        # And the final parameter (`absorp`) equals the following:
        # (gavage / injection) / { (nsa * sef) / (mass^0.75) }
        rhs = 1 / {(nsa * sef) / (mass^0.75)}
    ) %>% 
    group_by(diet, taxon, species) %>%
    summarize(
        # # The below n values for for absorp2 below. This parameter is not used currently.
        # n1 = ifelse(species[1] %in% sep_absorps, sum(!is.na(injection)), 0),
        # n2 = sum(!is.na(gavage)),
        # n3 = sum(!is.na(rhs)),
        
        # For sep_absorps species, I'm inversing injection here for the same reason
        # as for rhs above.
        # For non-sep_absorps species, I'm setting injection mean to 1 and variance to 0
        # bc the final value is already in the gavage column, and doing this makes the
        # equation below simplify to = V(Y), where Y is gavage
        inv_injection_v = ifelse(species[1] %in% sep_absorps,
                               var(1 / injection, na.rm = TRUE), 0),
        inv_injection = ifelse(species[1] %in% sep_absorps,
                               mean(1 / injection, na.rm = TRUE), 1),
        gavage_v = var(gavage, na.rm = TRUE),
        gavage = mean(gavage, na.rm = TRUE),
        # Variances for X*Y follow the following function:
        # V(X*Y) = E(X)^2 * V(Y) + E(Y)^2 * V(X) + V(X) * V(Y)
        # So I'm doing this for the above two parameters first, then doing it again
        # when combining it with `rhs` below
        # "ig_" stands for injection and gavage
        ig_v = inv_injection^2 * gavage_v + gavage^2 * inv_injection_v +
            inv_injection_v * gavage_v,
        ig = inv_injection * gavage,
        rhs_v = var(rhs, na.rm = TRUE),
        rhs = mean(rhs, na.rm = TRUE),
        # I added the sqrt bc I want the output to be the standard deviation
        absorp = sqrt(rhs^2 * ig_v + ig^2 * rhs_v + rhs_v * ig_v)#,
        # I'm just going with SD for now. The below equation doesn't make sense to me.
        # # Could this possibly work for the SE, or should I just use SD?
        # absorp2 = absorp / sqrt(n1^2 + n2^2 + n3^2)
    ) %>%
    ungroup %>% 
    select(taxon, species, absorp) %>%  # , absorp2) %>%
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')))

rm(sep_absorps)
