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
})
#' 
#' 
#' 
#' # Reading csv of morphometric measurements
#' 
#' Reading and cleaning `./data/morphometrics.csv` file for use.
#' 
#+ source_tidy_csv
source('tidy_csv.R')
morph_df
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
#' 
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
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Creating data frames of measurements
#' 
#' 
#' Now I create two data frames with the different measurements as columns, plus separate
#' columns for their log-transformed versions.
#' 
#' The first one summarizes by species only (and takes mean over positions), 
#' while the second one keeps the positions (distal, medial, proximal) separate.
#' Each data frame only needs certain measurements, so only those are included.
#' 
#' ## Species data frame
#' 
#+ make_spp_df
spp_measures <- c('mass',
                  'intestinal_length',
                  'nsa',
                  'villa_surface_area',
                  'enterocyte_density',
                  'sef')

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
#' 
#' 
#' ### Analyses for this data frame
#' 
#' `Model: Y ~ Taxon`
#' 
#' List of `Y`s:
#' 
#' - Intestinal length / body mass^0.4
#' - NSA / body mass^0.75
#' - Villus surface area / body mass^0.75
#' - Total number of enterocytes (body mass as covariate)
#'   * Calculated as such: `NSA * mean(enterocyte_density_for_all_positions)`
#' - Fractional absorption / (total intestinal surface / mass^0.75)
#'   * total intestinal surface = `NSA * SEF`
#' 
#' 
#' `Model: Y ~ Diet`
#' 
#' List of `Y`s:
#' 
#' - SEF
#' 
#+ sp_analyses

# (Since I don't yet have fractional absorption data, I'm skipping that for now)
spp_ys <- c("int_length_mass", "nsa_mass", "vill_area_mass", "log_total_enterocytes")

# # Takes ~6.5 min
# set.seed(88754829)
# spp_fits <- lapply(
#     spp_ys,
#     function(y) {
#         f <- paste(y, ' ~ taxon',
#                    ifelse(grepl('total_enterocytes', y), '+ log_mass', ''))
#         suppressWarnings(
#             do.call("phylolm", list(as.formula(f), data = as.name("spp_df"),
#                                     phy = as.name("tr"), model = 'lambda',
#                                     boot = 2000))
#         )
#     })
# names(spp_fits) <- spp_ys
# save(spp_fits, spp_df, file = './data/spp_models.rda')
load('./data/spp_models.rda')
lapply(spp_fits, summary)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ## Positions data frame
#' 
#+ make_pos_df
pos_measures <- c('mass',
                  'intestinal_diameter',
                  'villus_height', 
                  'villus_width',
                  'crypt_width',
                  'sef',
                  'enterocyte_diameter',
                  'enterocyte_density')

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
#' 
#' 
#' I need to do these analyses separately for each position because modelling 
#' within-species and within-individual variance due to position rather than measurement
#' error or process error would be difficult and not likely possible with this small 
#' dataset.
#' I'll now make three new `data.frame`s (rather than `tibble`s) 
#' of measurements for just one position, with species names as row names.
#' 
for (p in unique(pos_df$pos)) {
    out_df <- pos_df %>% 
        filter(pos == p) %>% 
        select(-pos) %>% 
        as.data.frame
    rownames(out_df) <- out_df$species
    assign(paste0(p, '_df'), out_df)
}; rm(p, out_df)

#' 
#' 
#' 
#' ### Analyses for this data frame
#' 
#' `Model: Y ~ Taxon * Segment`
#' 
#' (Segment = proximal, medial, or distal)
#' 
#' List of `Y`s:
#' 
#' - Intestinal diameter (body mass as covariate)
#' - Villus height
#' - Villus width
#' - Crypt width
#' - Surface enlargement factor (SEF)
#' - Enterocyte diameter
#' - Enterocytes per cm^2 NSA
#' 
#' 
#' 
#' 
#' 
#' 
pos_ys <- pos_measures[pos_measures != 'mass']
pos_ys[pos_ys == 'intestinal_diameter'] <- 'log_intestinal_diameter'
pos_ys[pos_ys == 'enterocyte_density'] <- 'log_enterocyte_density'






# # Took 21.7 minutes
# set.seed(25413535)
# pos_fits <- lapply(
#     unique(pos_df$pos),
#     function(pos) {
#         lapply(
#             pos_ys,
#             function(y) {
#                 f <- paste(y, ' ~ taxon',
#                            ifelse(grepl('intestinal_diameter|villus_height', y),
#                                   '+ log_mass', ''))
#                 # This model doesn't find the peak likelihood unless specifying a starting
#                 # value of 0.1.
#                 if (y == "log_enterocyte_density" & pos == "prox") {
#                     suppressWarnings(
#                         do.call("phylolm",
#                                 list(
#                                     as.formula(f),
#                                     data = as.name(paste0(pos, "_df")),
#                                     phy = as.name("tr"), model = 'lambda',
#                                     boot = 2000,
#                                     starting.value = 0.1)))
#                 } else {
#                     suppressWarnings(
#                         do.call("phylolm",
#                                 list(
#                                     as.formula(f),
#                                     data = as.name(paste0(pos, "_df")),
#                                     phy = as.name("tr"), model = 'lambda',
#                                     boot = 2000)))
#                 }
#             })
#     })
# names(pos_fits) <- unique(pos_df$pos)
# for (i in 1:length(pos_fits)) names(pos_fits[[i]]) <- pos_ys
# save(pos_fits, dist_df, med_df, prox_df, file = './data/pos_models.rda')

load('./data/pos_models.rda')
lapply(pos_fits$dist, summary)
lapply(pos_fits$med, summary)
lapply(pos_fits$prox, summary)






#' 
#' 
#' 
#' 
#' # Model: log(Y) ~ log(SEF)
#' 
#' - paracellular probe L-arabinose clearance
#'   ("we used reduced major axis regression (model II regression)... because both 
#'   variables [X and Y] were subject to error")
#' 
#' 
#' 
#' 
#' 
#' 
#+ make_sp_df



# # NOT SURE BELOW IS USEFUL, BUT IT *MIGHT* BE. SO I LEFT IT.
# sp_df <- prep_df(measures = unique(morph_df$measure))
# str(sp_df)
# 
# 
# phylolm
# 
# 
# 
# to_exam <- c('crypt_width', 'enterocyte_density', 'enterocyte_width', 
#              'intestinal_diameter', 'intestinal_length', 'nsa', 'sef',
#              'villa_surface_area', 'villus_height', 'villus_width')
# plot.new()
# par(mfrow = c(1, 2), mar=c(5.1, 4.1, 1, 1))
# for (p in to_exam) {
#     plot(sp_df[['log_mass']], sp_df[[paste0(p)]], 
#          ylab = paste0(p), xlab = 'log(mass)', main = NULL)
#     plot(sp_df[['log_mass']], sp_df[[paste0(p, '_log')]], 
#          ylab = paste0(p, '_log'), xlab = 'log(mass)', main = NULL)
# }; rm(p)
# 
# 
# # Keep logged: villa_surface_area, nsa, intestinal_length, intestinal_diameter

#' 
#' 
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
                      phylolm(nsa_log ~ log_mass + taxon, data = sp_df, phy = tr,
                              model = m, boot = 2000, 
                              upper.bound = ifelse(m == 'lambda', 1.2, Inf))})
names(nsa_fits) <- c('lambda', 'ou')
sef_fits <- lapply(c('lambda', 'OUfixedRoot'), 
                   function(m) {
                       phylolm(sef_log ~ log_mass + taxon, data = sp_df, phy = tr,
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
    mass = sapply(nsa_fits, function(m) mean(m$bootstrap[,'log_mass'] < 0) * 2)
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
    mass = sapply(sef_fits, function(m) mean(m$bootstrap[,'log_mass'] < 0) * 2)
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

