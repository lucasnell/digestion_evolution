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
    library(ggtree)
})
#' 
#' 
#' 
#' # Morphometric measurements
#' 
#' Cleaning `./rd_files/morphometrics.csv` file for use and providing a useful 
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
tr <- read.tree('tree.nwk')
tr$tip.label <- gsub('_', ' ', tr$tip.label)
tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% (morph_df$species %>% unique)])
tr
#' 
#' 
#' 
#' 
#' ## Visualizing tree
#' 
#' Here is the phylogenetic tree with log(NSA) as tip color and log(SEF) as tip size.
#' 
#' 
#+ phylo_plot, fig.width=8, fig.height=6, echo = FALSE
x_end <- 32
gg_tr <- ggtree(tr)
gg_tr$data$x <- gg_tr$data$x - max(gg_tr$data$x)
gg_tr %<+% {sp_df %>% select(species, everything())} +
    theme_tree2(axis.title.x = element_text(size = 14),
                legend.position = c(0.25, 0.75), legend.box = 'horizontal',
                legend.background = element_rect(color = NA, fill = NA)) +
    geom_tiplab(aes(x = x + 3), size = 3, fontface = 'bold.italic') +
    geom_tippoint(aes(x = x + 1, size = sef_log, color = nsa_log)) +
    geom_text(data = data_frame(x = rep(x_end + 1, 2), y = c(14, 5), 
                                label = c('Rodents', 'Bats')), 
              aes(label = label), angle = -90, vjust = 0, size = 6) +
    geom_segment(data = data_frame(x = rep(x_end, 2), xend = rep(x_end, 2), 
                                   y = c(1, 10), yend = c(9, 18)), 
              aes(xend = xend, yend = yend)) +
    scale_x_continuous('Time (mya)', limits = c(-100, x_end + 2),
                       breaks = seq(-100, 0, 25), labels = seq(100, 0, -25)) +
    scale_color_gradient('log(NSA)', low = 'darkblue', high = 'cadetblue1',
                         guide = guide_colorbar(direction = "horizontal",
                                                title.position = 'top', 
                                                title.hjust = 0.5)) +
    scale_size_continuous('log(SEF)', 
                         guide = guide_legend(direction = "horizontal",
                                              title.position = 'top', 
                                              title.hjust = 0.5))
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
#' 
#' I also ran 2,000 parametric bootstrap replicates to estimate model parameters.
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
save(nsa_fits, sef_fits, file = 'model_fits.RData', compress = FALSE)



#+ load_fits, echo = FALSE
load('model_fits.RData')

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
#' *Ornstein-Uhlenbeck*
#' 
#+ ua_sef, echo = FALSE
summary(sef_fits[['ou']])


#' 
#' ## P-values
#' 
#' These are p-value based on bootstrap replicates for whether the coefficient for the 
#' taxon covariate is significantly different from zero.
#' 
#' ### `nsa`
#' 
#+ nsa_ps, echo = FALSE
p_nsa_lambda <- mean(nsa_fits[['lambda']]$bootstrap[,'taxonRodent'] < 0) * 2
p_nsa_ou <- mean(nsa_fits[['ou']]$bootstrap[,'taxonRodent'] < 0) * 2
cat("P for Pagel's lambda     =", format(p_nsa_lambda, nsmall = 3))
cat("P for Ornstein-Uhlenbeck =", format(p_nsa_ou, nsmall = 3))
#' 
#' 
#' ### `sef`
#' 
#+ sef_ps, echo = FALSE
p_sef_lambda <- mean(sef_fits[['lambda']]$bootstrap[,'taxonRodent'] > 0) * 2
p_sef_ou <- mean(sef_fits[['ou']]$bootstrap[,'taxonRodent'] > 0) * 2
cat("P for Pagel's lambda     =", format(p_sef_lambda, nsmall = 3))
cat("P for Ornstein-Uhlenbeck =", format(p_sef_ou, nsmall = 3))




#' 
#' # Session info
#' 
#' This outlines the package versions I used for these analyses.
#' 

devtools::session_info()

