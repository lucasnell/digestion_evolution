
# This reads the initial Excel file and simplifies it into a csv file

library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)


xl <- read_excel('./data/raw_data.xlsx', col_names = FALSE, col_types = rep('text', 88))

# Table of abbreviated species names
spp_df <- xl[69:86,2:3] %>% 
    rename(abbrev = X__2, full = X__3) %>% 
    arrange(full)
# Resolving ambiguous and untidy abbreviated names
spp_df$abbrev[spp_df$abbrev == 'Ds1'] <- 'Ds'
spp_df$abbrev[spp_df$full == 'Eumops glaucinus'] <- 'Eg'
spp_df$abbrev[spp_df$full == 'Microtus pennsylvanicus'] <- 'Microtus'
spp_df$abbrev[spp_df$full == 'Molossus molossus'] <- 'Mmo'



# Start of final morphometrics data frame, starting with diet, taxon, name, and 
# villus heights
initial_df <- xl[c(4:34,36:64),1:6] %>%  # (<-- row 35 is all NAs)
    rename_(.dots = setNames(paste0('X__', 1:6), 
                             c('diet', 'taxon', 'id', 
                               'prox', 'med', 'dist'))) %>% 
    mutate_at(vars(prox, med, dist), 
              function(x) as.numeric(ifelse(x == 'None', NA, x))) %>% 
    mutate(measure = 'villus height')

# Resolving ambiguous ids
initial_df$id[initial_df$taxon == 'Bat' & initial_df$id == 'E'] <- 'Eg'
initial_df$id[initial_df$taxon == 'Bat' & grepl('Mm', initial_df$id)] <- paste0('Mmo', 1:3)
initial_df$id[grepl('My', initial_df$id)] <- paste0('My', 1:3)

# Assigning species names
initial_df <- initial_df %>% 
    mutate(species = sapply(id, function(.id) {
        spp_df$full[spp_df$abbrev == gsub('[0-9]', '', .id)]
    })) %>% 
    select(diet, taxon, species, id, measure, everything())



# Function to get integer index from Excel alphabetic index
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


# I'm going to make a new data frame of column names along with their column location
# in the Excel file and the number of additional columns to take (after the one 
# specified).
# It's assumed that all these are located in rows c(4:34,36:64), like the villus height
# data.
# Also, n should only be 2 or 0, for measures with and without separate proximal, 
# medial, and distal measurements, respectively.
# Lastly, although I am extracting columns from the Excel file, I want the new csv file
# to be in "long" format, so I'm binding them by rows.

new_cols <- read_csv('name,col1,n
"villus width",M,2
"crypt width",U,2
"enterocyte density",AC,2
sef,AM,2
"intestinal diameter",AU,2
"intestinal length",BB,0
mass,BH,0
"villa surface area",BV,0
"enterocyte width",CF,2
nsa,BE,0')

new_cols <- new_cols %>% 
    mutate(col1 = as.integer(sapply(col1, .i) - 1))


# Function to extract info from the xl data frame based on a row from new_cols
xl_extr <- function(.x) {
    .df <- xl[c(4:34,36:64), .x$col1:(.x$col1+.x$n)]
    if (.x$n == 0) {
        .df <- .df %>% 
            mutate_(.dots = setNames(list(as.numeric(NA), as.numeric(NA)), 
                                     paste0('X__', .x$col1 + 1:2)))
    }
    
    .df <- .df %>%
        rename_(.dots = setNames(paste0('X__', .x$col1:(.x$col1 + 2)), 
                                 c('prox', 'med', 'dist'))) %>% 
        mutate_all(function(x) as.numeric(ifelse(x == 'None', NA, x))) %>% 
        mutate(measure = .x$name)
    # Now adding columns common to all dataframes
    out_df <- bind_cols(initial_df %>% select(diet, taxon, species, id), .df) %>% 
        select(diet, taxon, species, id, measure, everything())
    
    return(out_df)
}

morph_df <- new_cols %>% 
    split(.$name) %>% 
    map(xl_extr) %>% 
    bind_rows %>% 
    bind_rows(initial_df) %>% 
    arrange(measure, species, id)


write_csv(morph_df, './data/clean_data.csv')
