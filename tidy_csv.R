library(readr)
library(dplyr)
library(tidyr)

morph_df <- read_csv('morphometrics.csv', col_types = 'cccccddd')


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

