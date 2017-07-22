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
# Set the default ggplot theme
theme_set(theme_classic() %+replace% 
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 12),
                    legend.background = element_blank()))
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
    select(taxon, diet, species, id,
           int_length_mass, nsa_mass, vill_area_mass, 
           total_enterocytes, log_total_enterocytes, total_surface, sef, 
           mass, log_mass) %>% 
    # Taking mean by sample
    group_by(taxon, diet, species, id) %>% 
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

# Takes ~6.5 min
# model_fits <- lapply(
#     spp_ys,
#     function(x){
#         f <- paste(x, ' ~ taxon',
#                    ifelse(grepl('total_enterocytes', x), '+ log_mass', ''))
#         suppressWarnings(
#             do.call("phylolm", list(as.formula(f), data = as.name("spp_df"),
#                                     phy = as.name("tr"), model = 'lambda',
#                                     boot = 2000, upper.bound = 1.2))
#         )
#     })
# save(model_fits, file = './data/spp_models.rda')
load('./data/spp_models.rda')
lapply(model_fits, summary)




# Creates data frame containing 95% CI based on bootstrapping for one model
mod_ci <- function(.model, y_measure, mod_name = NULL){
    ci_matrix <- .model$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(head(x, -2) %*% t(.model$X))
              }) %>% 
        apply(1,
              function(x) as.numeric(quantile(x, probs = c(0.025, 0.975)))) %>% 
        t
    
    out_df <- cbind(as.numeric(.model$fitted.values), ci_matrix) %>% 
        as_data_frame %>% 
        rename_(.dots = setNames(colnames(.), c('predicted', 'low', 'high'))) %>% 
        mutate(measure = y_measure,
               log_mass = tryCatch(.model$X[,'log_mass'], 
                                   error = function(e) rep(NA, nrow(.model$X))),
               taxon = factor(.model$X[,'taxonBat'], levels = c(0,1), 
                              labels = c('Rodent', 'Bat'))) %>% 
        select(taxon, log_mass, measure, everything())
    if (!is.null(mod_name)) {
        out_df <- out_df %>% 
            mutate(model = mod_name) %>% 
            select(model, taxon, log_mass, measure, everything())
    }
    # Check if there is no need for log_mass or many of the rows
    if (ncol(.model$bootstrap) == 4) {
        out_df <- out_df %>% 
            select(-log_mass) %>% 
            distinct(taxon, measure, predicted, low, high)
    }
    return(out_df)
}





# Figure 1A
# fig1a <- 
mod_ci(model_fits[[1]], 'int_length_mass') %>% 
    ggplot(aes(taxon, predicted)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, size = 0.5) +
    geom_segment(aes(yend = predicted, 
                     x = as.numeric(taxon) - 0.05, 
                     xend = as.numeric(taxon) + 0.05)) +
    geom_point(data = spp_df, 
               aes(y = int_length_mass, shape = taxon),
               position = position_jitter(width = 0.2, height = 0),
               color = 'black', size = 2) +
    geom_text(data = NULL, label = 'A', x = 0.5, y = 11.7, size = 6, 
              vjust = 1, hjust = 0) +
    scale_shape_manual(values = c(19, 1)) +
    theme(legend.position = 'none', axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0))) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(expression(atop("Intestinal length / body" ~ mass^{0.4},
                                       "(" * cm / g^{0.4} * ")")),
                       limits = c(2.2, 11.7), expand = c(0, 0))

# Figure 1B
# fig1b <- 
mod_ci(model_fits[[2]], 'nsa_mass') %>% 
    ggplot(aes(taxon, predicted)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, size = 0.5) +
    geom_segment(aes(yend = predicted, 
                     x = as.numeric(taxon) - 0.05, 
                     xend = as.numeric(taxon) + 0.05)) +
    geom_point(data = spp_df, 
               aes(y = nsa_mass, shape = taxon),
               position = position_jitter(width = 0.2, height = 0),
               color = 'black', size = 2) +
    geom_text(data = NULL, label = 'B', x = 0.5, y = 2.2, size = 6, 
              vjust = 1, hjust = 0) +
    scale_shape_manual(values = c(19, 1)) +
    theme(legend.position = 'none', axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = -8, b = 0, l = 0))) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(expression(atop("NSA / body" ~ mass^{0.75},
                                       "(" * cm^2 / g^{0.75} * ")")),
                       limits = c(0.5, 2.2), expand = c(0, 0))
#

# Figure 4
# fig4 <- 
mod_ci(model_fits[[3]], 'vill_area_mass') %>% 
    ggplot(aes(taxon, predicted)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, size = 0.5) +
    geom_segment(aes(yend = predicted, 
                     x = as.numeric(taxon) - 0.05, 
                     xend = as.numeric(taxon) + 0.05)) +
    geom_point(data = spp_df, 
               aes(y = vill_area_mass, shape = taxon),
               position = position_jitter(width = 0.2, height = 0),
               color = 'black', size = 2) +
    scale_shape_manual(values = c(19, 1)) +
    theme(legend.position = 'none', axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = -8, b = 0, l = 0))) +
    scale_y_continuous(expression(atop("Villous surface area / body" ~ mass^{0.75},
                                       "(" * cm^2 / g^{0.75} * ")")))






# Figure 6
# fig6 <- 
mod_ci(model_fits[[4]], 'log_total_enterocytes') %>% 
    ggplot(aes(log_mass, predicted)) + 
    geom_ribbon(aes(group = taxon, ymin = low, ymax = high), 
                fill = 'gray70', alpha = 0.5) +
    geom_point(data = spp_df, aes(y = log_total_enterocytes, shape = taxon),
               color = 'black', size = 2) +
    geom_line(aes(linetype = taxon)) +
    theme(legend.position = c(0.15, 0.9), legend.title = element_blank(),
          legend.key.width = unit(0.05, "npc")) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_shape_manual(values = c(19, 1)) +
    scale_y_continuous(expression("Total enterocytes (" %*% 10^9 * ")"), 
                       breaks = log(seq(0, 2e9, 5e8)), labels = seq(0, 2, 0.5)) +
    scale_x_continuous("Body mass (g)", breaks = log(seq(0, 150, 50)), 
                       labels = seq(0, 150, 50)) +
    coord_trans(x = 'exp', y = 'exp')






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
    # Doing the transformation now, before taking any means
    mutate_(.dots = setNames(as.list(sprintf('%s(%s)', 'log', pos_measures)), 
                             paste0('log_', pos_measures))) %>% 
    # Grouping by, then taking mean of all measurement columns and transformed-
    # measurement columns
    group_by(taxon, diet, species, pos) %>% 
    summarize_at(.vars = c(pos_measures, paste0('log_', pos_measures)), mean) %>% 
    ungroup %>% 
    # Converting taxon and diet to factors
    mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
           diet = factor(diet, levels = c("Herbivorous", "Omnivorous", "Protein"))) %>% 
    arrange(taxon, diet, species, pos)
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
#' ### Model: log(Y) ~ log(SEF)
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
sp_df <- prep_df(measures = unique(morph_df$measure))
str(sp_df)




to_exam <- c('crypt_width', 'enterocyte_density', 'enterocyte_width', 
             'intestinal_diameter', 'intestinal_length', 'nsa', 'sef',
             'villa_surface_area', 'villus_height', 'villus_width')
plot.new()
par(mfrow = c(1, 2), mar=c(5.1, 4.1, 1, 1))
for (p in to_exam) {
    plot(sp_df[['log_mass']], sp_df[[paste0(p)]], 
         ylab = paste0(p), xlab = 'log(mass)', main = NULL)
    plot(sp_df[['log_mass']], sp_df[[paste0(p, '_log')]], 
         ylab = paste0(p, '_log'), xlab = 'log(mass)', main = NULL)
}; rm(p)


# Keep logged: villa_surface_area, nsa, intestinal_length, intestinal_diameter

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

