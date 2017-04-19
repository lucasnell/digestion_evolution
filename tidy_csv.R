
# This cleans the csv file for use and provides a useful function to retrieve columns 
# from it

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

# Gathering into 'tall' format and fixing position column
morph_df <- morph_df %>% 
    gather(pos, value, prox:dist, na.rm = TRUE) %>% 
    mutate(pos = ifelse(measure %in% no_pos, NA, pos))

# This object is no longer necessary
rm(no_pos)



# Function to return 'wide' version of an input data frame of morphometric measurements 
# (defaults to morph_df), given a character vector of measures
prep_df <- function(measures, input_df = morph_df, by_sp = TRUE, 
                    trans_fun = 'log'){
    
    # List of measures wrapped in ticks (``) so they'll play nice with dplyr and tidyr
    # as column names even if they have spaces
    meas_list <- lapply(measures, function(s) sprintf('`%s`', s))
    # Replacing all spaces with underscores for better naming in the output data frame
    meas_clean <- gsub(' ', '_', measures)
    # Names of transformed columns, also with underscores rather than spaces
    trans_meas <- paste(meas_clean, trans_fun, sep = '_')
    
    new_df <- input_df %>% 
        # Changing from tall to wide format
        spread(measure, value) %>% 
        # Selecting measurement columns, plus the identifying columns
        select_(.dots = append(list('diet', 'taxon', 'species', 'id'), 
                               meas_list)) %>% 
        # Removing all rows with all NAs in measures columns
        filter(Reduce(`+`, lapply(.[,measures], is.na)) < length(measures)) %>% 
        # Replacing spaces in measures-column names with underscores
        rename_(.dots = setNames(as.list(sprintf('`%s`', measures)), meas_clean)) %>% 
        # Doing the transformation now, before taking any means
        mutate_(.dots = setNames(as.list(sprintf('%s(%s)', trans_fun, meas_clean)), 
                                 trans_meas)) %>% 
        # Taking mean by sample
        group_by(diet, taxon, species, id) %>% 
        summarize_all(mean, na.rm = TRUE) %>% 
        ungroup
    
    if (by_sp) {
        new_df <- new_df %>% 
            # Grouping by, then taking mean of all measurement columns and transformed-
            # measurement columns
            group_by(diet, taxon, species) %>% 
            summarize_at(.cols = c(trans_meas, meas_clean), mean) %>% 
            ungroup %>% 
            arrange(taxon, diet, species) %>% 
            # To change row names, it can't be a tibble, so I'm reverting back to normal
            # data frame
            as.data.frame
        # phylolm requires that the rownames match the species names
        rownames(new_df) <- new_df$species
    } else {
        new_df <- new_df %>% select(species, everything()) %>% 
            as.data.frame
    }
    
    return(new_df)
}


