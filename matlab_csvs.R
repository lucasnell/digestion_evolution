suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(phylolm)
    library(ape)
    library(ggplot2)
})

source('tidy_csv.R')
source('tidy_csv_se.R')

spp_df <- prep_df(measures = c('nsa', 'sef', 'mass'), by_sp = FALSE) %>% 
    as_data_frame %>% 
    mutate(taxon = ifelse(taxon == 'Rodent', 1, 0), species = factor(species))

spp_df

# data file: row for each species
# me file: se for each cell in data file
# file: physigv2.m file
# run it, and dialogs will show up
# taxon me = 0
# vcov txt file for phylo file
# use likelihoods for p-values: run it with and with, 2*logLik is chi-square
# run it with ml, not reml


sp_se <- spp_df %>%
    group_by(species) %>% 
    summarize_at(vars(nsa_log, sef_log, mass_log, taxon), function(x) sd(x) / length(x)) %>% 
    mutate_at(vars(nsa_log, sef_log, mass_log, taxon),
              function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)) %>% 
    select(-species)

# write_tsv(sp_se, 'data/matlab/se.txt', col_names = FALSE)

sp_mean <- spp_df %>%
    group_by(species) %>% 
    summarize_at(vars(nsa_log, sef_log, mass_log, taxon), 'mean') %>% 
    select(-species)

# write_tsv(sp_mean, 'data/matlab/data.txt', col_names = FALSE)



tr <- read.tree('./data/tree.nwk')
tr$tip.label <- gsub('_', ' ', tr$tip.label)
tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% (morph_df$species %>% unique)])
tr

length(tr$tip.label) == length(unique(spp_df$species))


tr_df <- vcv(tr)[sort(tr$tip.label),sort(tr$tip.label)] %>% 
    as_data_frame %>% 
    magrittr::set_colnames(paste0('sp', 1:ncol(.)))
tr_df

# write_tsv(tr_df, 'data/matlab/vcv.txt', col_names = FALSE)
