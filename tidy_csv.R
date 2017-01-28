library(readr)
library(dplyr)
library(tidyr)

morph_df <- read_csv('morphometrics.csv', col_types = 'cccccddd') %>% 
    # Oligoryzomys seems to be the more standard spelling
    mutate(species = ifelse(species == 'Olygoryzomys nigripes', 
                            'Oligoryzomys nigripes', species))




N <- morph_df$individual %>% unique %>% length

# Measures with no position
no_pos <- morph_df %>% 
    filter(is.na(dist)) %>% 
    group_by(measure) %>% 
    summarize(total = n()) %>% 
    filter(total == N) %>% 
    select(measure) %>% 
    unlist %>% 
    paste

# Gathering into 'tall' format and fixing position column
morph_df <- morph_df %>% 
    gather(pos, value, prox:dist, na.rm = TRUE) %>% 
    mutate(pos = ifelse(measure %in% no_pos, NA, pos))

# This object is no longer necessary
rm(no_pos)


# ============
# Formulas from Excel sheet
# ============

# # Total enterocytes
# { prox_enterocytes * (NSA / 3) } +{ med_enterocytes * (NSA / 3) } + { dist_enterocytes * (NSA / 3) }
# 
# # Intestinal length/bodymass^0,4
# intestinal_length / body_mass^0.4
# 
# # Intestinal diameter/body mass^0,75
# # This was done for proximal, medial, and distal intestinal diameters separately
# intestinal_diameter / body_mass^0.75
# 
# # NSA/Body mass0,75
# nsa / body_mass^0.75
# 
# # Mucosa / Body mass0,75
# total_villa_surface_area / body_mass^0.75

