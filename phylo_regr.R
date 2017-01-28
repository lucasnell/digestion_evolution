# Regression with phylogenetic covariance

source('tidy_csv.R')


library(phylolm)
library(ape)
library(nlme)
library(ggplot2)
# "... the best options are ape (that uses gls in nlme), phylolm, or phytools." (Tony)


tr <- read.tree('tree.nwk')
tr$tip.label <- gsub('_', ' ', tr$tip.label)
tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% (morph_df$species %>% unique)])


# Return 'wide' version of morph_df, given a character vector of measures
prep_df <- function(measures, df = morph_df, by_sp = TRUE){
    
    meas_list <- lapply(measures, function(s) sprintf('`%s`', s))
    meas_clean <- gsub(' ', '_', measures)
    
    new_df <- df %>% 
        spread(measure, value) %>% 
        select_(.dots = append(list('diet', 'taxa', 'species', 'individual'), 
                               meas_list)) %>% 
        # Removing all rows with all NAs in measures columns
        filter(Reduce(`+`, lapply(.[,measures], is.na)) < length(measures)) %>% 
        # Replacing spaces in measures-column names with underscores
        rename_(.dots = setNames(as.list(sprintf('`%s`', measures)), meas_clean))
    
    if (by_sp) {
        new_df <- new_df %>% 
            group_by(diet, taxa, species) %>% 
            summarize_at(.cols = meas_clean, mean) %>% 
            ungroup
    }
    new_df <- new_df %>% 
        arrange(taxa, diet, species) %>% 
        # To change row names, it can't be a tibble
        as.data.frame
    rownames(new_df) <- new_df$species
    
    return(new_df)
}




nsa_bm_df <- prep_df(c('nsa', 'body mass'))

phy_mod <- gls(log(nsa) ~ log(body_mass) + taxa, data = nsa_bm_df, 
               correlation = corBrownian(phy = tr))

phy_mod2 <- phylolm(log(nsa) ~ log(body_mass) + taxa, data = nsa_bm_df, phy = tr, 
                    model = 'BM', boot = 1000)

summary(phy_mod)
summary(phy_mod2)

plot(phy_mod)
plot(phy_mod2)


nsa_bm_df %>% 
    ggplot(aes(log(body_mass), log(nsa), color = taxa)) +
    geom_point() +
    geom_smooth(method = 'lm')


