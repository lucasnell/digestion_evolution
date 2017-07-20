#' ---
#' title: "Phylogenetic linear regression"
#' author: "Lucas Nell"
#' date: "`r format(Sys.Date(), '%d %b %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#' ---
#' 
#+ setup, include = FALSE, cache = FALSE
knitr::opts_chunk$set(echo = TRUE)
#' 
#' 
#' 
#' 
#' Loading packages:
#' 
#+ load_packages
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(magrittr)
    library(phylolm)
    library(ape)
    library(ggplot2)
})
#' 
#' 
#' 
#' # Morphometric measurements
#' 
#' Cleaning `./data/morphometrics.csv` file for use and providing a useful 
#' function to retrieve columns from it
#' 
#+ source_tidy_csv
source('tidy_csv.R')
str(morph_df)
#' 
#' 
#' 
#' 
#' All measures found in `morph_df`:
#' 
#' - `crypt width`
#' - `enterocyte density`
#' - `enterocyte width`
#' - `intestinal diameter`
#' - `intestinal length`
#' - `mass`
#' - `nsa`
#' - `sef`
#' - `villa surface area`
#' - `villus height`
#' - `villus width`
#' 
#' 
#' 
#' I am only using three: `nsa`, `sef`, and `mass`. Now I create a data frame with 
#' just these columns and their log-transformed versions. Species names are row names
#' because `phylolm` requires that.
#' 
#' 
#' 
#+ make_sp_df
sp_df <- prep_df(measures = c('nsa', 'sef', 'mass'))
str(sp_df)
#' 
#' 
#' 
#' 
#' 
#' 
#' # Phylogenetic tree
#' 
#' Reading phylogenetic tree, cleaning species names, and removing unnecessary species 
#' from it
#' 
#+ make_tr
tr <- read.tree('./data/tree.nwk')
tr$tip.label <- gsub('_', ' ', tr$tip.label)
tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% (morph_df$species %>% unique)])
tr

all_df <- prep_df(measures = c('nsa', 'sef', 'mass'), by_sp = FALSE) %>% 
    as_data_frame

library(Rphylopars)

library(phytools)


# set.seed(9)
# test <- all_df %>% 
#     select(species, taxon) %>% 
#     mutate(mass_log = log(rlnorm(nrow(test), mean(all_df$mass_log), sd(all_df$mass_log)))) %>% 
#     mutate(nsa_log = ifelse(taxon == 'Rodent', 1, 0) * 0.5 + mass_log * 0.6,
#            taxon = factor(taxon)) %>% 
#     select(species, nsa_log, mass_log, taxon)

# z <- phylopars.lm(nsa_log ~ mass_log + taxon, 
#                   trait_data = as.data.frame(test), 
#                   tree = tr, model = 'BM')
# 
# z <- phylopars.lm(nsa_log ~ mass_log + taxon, 
#                   trait_data = as.data.frame(
#                       all_df[,c('species', 'nsa_log', 'mass_log', 'taxon')]
#                   ), 
#                   tree = tr, model = 'BM')
# summary(z)

# data file: row for each species
# me file: se for each cell in data file
# file: physigv2.m file
# run it, and dialogs will show up
# taxon me = 0
# vcov txt file for phylo file
# use likelihoods for p-values: run it with and with, 2*logLik is chi-square
# run it with ml, not reml



X <- all_df$mass_log
names(X) <- all_df$species
y <- all_df$nsa_log
names(y) <- all_df$species
pgls.Ives(tr, X, y)

sp_se

pgls.SEy(nsa_log ~ mass_log + taxon, data = sp_df, 
         corClass=corBrownian, tree=tr,
         se=NULL, method="ML", interval=c(0,1000))

m <- phylolm(nsa_log ~ mass_log + taxon, data = sp_df, phy = tr,
             model = 'lambda', upper.bound = 1.2)
summary(m)





#' 
#' 
#' 
#' 
#' 
#' # Fitting phylogenetic linear regression models
#' 
#' Below fits phylogenetic linear regression models using `phylolm::phylolm`. For both
#' `nsa` and `sef` (log-transformed), I fit models using log-transformed mass and
#' taxon (a factor based on whether that species is a rodent or bat) as covariates. 
#' (I tried including the interaction between mass and taxon, but it increased the 
#' AIC in all models.)
#' 
#' I fit two types of phylogenetic-covariance models for both `nsa` and `sef` regression 
#' models: "the Ornstein-Uhlenbeck model with an ancestral state to be estimated at the 
#' root (OUfixedRoot) ... [and] Pagel's lambda model." (see `phylolm` documentation)
#' As you can see from the results below, the covariance model had little effect on our
#' conclusions.
#' The Ornstein-Uhlenbeck model gave much less precise estimates of the 
#' phylogenetic signal (see the bootstrapped confidence intervals of the alpha 
#' parameter below).
#' This is why I only mentioned Pagel's lambda in the manuscript.
#' 
#' I ran 2,000 parametric bootstrap replicates to estimate model parameters.
#' 
#' The code below takes ~10 minutes to run.
#' 
#' 
#+ fit_mods, eval = FALSE
set.seed(352)
nsa_fits <- lapply(c('lambda', 'OUfixedRoot'), 
                  function(m) {
                      phylolm(nsa_log ~ mass_log + taxon, data = sp_df, phy = tr,
                              model = m, boot = 2000, 
                              upper.bound = ifelse(m == 'lambda', 1.2, Inf))})
names(nsa_fits) <- c('lambda', 'ou')
sef_fits <- lapply(c('lambda', 'OUfixedRoot'), 
                   function(m) {
                       phylolm(sef_log ~ mass_log + taxon, data = sp_df, phy = tr,
                               model = m, boot = 2000,
                               upper.bound = ifelse(m == 'lambda', 1.2, Inf))})
names(sef_fits) <- c('lambda', 'ou')
save(nsa_fits, sef_fits, file = './data/model_fits.RData', compress = FALSE)



#+ load_fits, echo = FALSE
load('./data/model_fits.RData')

#' 
#' # Model output
#' 
#' ## Summaries
#' 
#' ### `nsa`
#' 
#' *Pagel's lambda*
#' 
#+ lambda_nsa, echo = FALSE
summary(nsa_fits[['lambda']])
#' 
#' 
#' 
#' *Ornstein-Uhlenbeck*
#' 
#+ ua_nsa, echo = FALSE
summary(nsa_fits[['ou']])



#' 
#' ### `sef`
#' 
#' *Pagel's lambda*
#' 
#+ lambda_sef, echo = FALSE
summary(sef_fits[['lambda']])
#' 
#' 
#' 
#' 
#' *Ornstein-Uhlenbeck*
#' 
#+ ua_sef, echo = FALSE
summary(sef_fits[['ou']])


#' 
#' ## P-values
#' 
#' These are p-value based on bootstrap replicates for whether the coefficient for the 
#' each covariate is significantly different from zero.
#' 
#' ### `nsa`
#' 
#+ nsa_ps, echo = FALSE
p_nsa_df <- data_frame(
    model = c("Pagel's lambda", "Ornstein-Uhlenbeck"),
    taxon = sapply(nsa_fits, 
                   function(m) mean(m$bootstrap[,'taxonRodent'] < 0) * 2),
    mass = sapply(nsa_fits, function(m) mean(m$bootstrap[,'mass_log'] < 0) * 2)
)
knitr::kable(p_nsa_df)
#' 
#' 
#' ### `sef`
#' 
#+ sef_ps, echo = FALSE
p_sef_df <- data_frame(
    model = c("Pagel's lambda", "Ornstein-Uhlenbeck"),
    taxon = sapply(sef_fits, 
                   function(m) mean(m$bootstrap[,'taxonRodent'] > 0) * 2),
    mass = sapply(sef_fits, function(m) mean(m$bootstrap[,'mass_log'] < 0) * 2)
)
knitr::kable(p_sef_df)
#' 
#' 
#' 
#' 
#' 
#' 
#' # Session info
#' 
#' This outlines the package versions I used for these analyses.
#' 

devtools::session_info()

