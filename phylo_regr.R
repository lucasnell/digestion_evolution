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
    library(phytools)
})
#' 
#' 
#' 
#' # Reading and cleaning csv files
#' 
#+ source_tidy_csv
source('tidy_csv.R')
#' 
#' 
#' ## Morphometric data
#' 
#' In 'tidy_csv.R` I first create two data sets from morphometric data in 
#' `data/clean_morph_data.csv`
#' with the different measurements as columns, plus separate columns for their 
#' log-transformed versions.
#' Each data set only needs certain measurements, so only those are included.
#' All data are present in `data.frame`s (rather than `tibble`s), with species names 
#' as row names.
#' 
#' The first data set consists of two `data.frame`s, each summarizing by species 
#' only (and takes means over positions).
#' The only difference between the two data frames is that one is designed for the
#' analysis of SEF in relation to diet. Thus it has species with no diet data removed.
#'  
#' The second data set consists of separate data frames for morphometric data for
#' each position (distal, medial, proximal).
#' I need to do my analyses separately for each position because modelling 
#' within-species and within-individual variance due to position rather than measurement
#' error or process error would be difficult and not likely possible with this small 
#' dataset.
#' 
#' 
#' ## Absorption and clearance data
#' 
#' `tidy_csv.R` also prepares absorption and clearance data 
#' (`data/clean_absorption_data.csv` and `data/clean_clearance_data.csv`, respectively).
#' Each of these is organized into a single `data.frame` with species names as row 
#' names.
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Phylogenetic trees
#' 
#' Reading phylogenetic tree, cleaning species names, and removing unnecessary species 
#' from it for each set of analyses that differs in the set of unique species.
#' These analyses are for `SEF ~ diet`, `clearance ~ SEF`, and 
#' `fractional absorption ~ taxon`, respectively.
#' The `tr` phylogeny is for all other analyses.
#' 
#+ make_trs
tr <- read.tree('data/tree.nwk')
tr$tip.label <- gsub('_', ' ', tr$tip.label)
diet_tr <- drop.tip(tr,  tip = tr$tip.label[!tr$tip.label %in% diet_df$species])
clear_tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% clear_df$species])
absorp_tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% absorp_df$species])
tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% spp_df$species])
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
#' 
#' ## Species data frame
#' 
#+ make_spp_df
spp_measures

# Function to get a p value from a bootstrapped phylolm model
pval <- function(model, parameter = 'taxonBat') {
    2 * min(mean(model$bootstrap[,parameter] < 0), 
            mean(model$bootstrap[,parameter] > 0))
}
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



set.seed(581120)
diet_mod <- phylolm(sef ~ diet, data = diet_df, phy = diet_tr, 
                    model = 'lambda', boot = 2000)
summary(diet_mod)
pval(diet_mod, 'dietOmnivorous'); pval(diet_mod, 'dietProtein')

#' 
#' 
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
# save(spp_fits, spp_df, diet_mod, file = 'data/spp_models.rda')
load('data/spp_models.rda')
sapply(spp_fits, pval)
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
pos_measures

#' 
#' 
#' 




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
#                 # This model doesn't find the peak likelihood unless specifying a 
#                 # starting value of 0.1.
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
# # (Saving the tree too in case you want to re-fit any models in other files)
# save(pos_fits, dist_df, med_df, prox_df, tr, file = 'data/pos_models.rda')

load('data/pos_models.rda')
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
#+ clear_sef


X <- log(clear_df$sef)
Y <- log(clear_df$clear)
names(X) <- names(Y) <- rownames(clear_df)
clear_rma <- phyl.RMA(X, Y, clear_tr, method = 'lambda')
clear_rma
plot(clear_rma)
#' 
#' 
#' 
#' 
#+ absorp_taxon

absorp_df



set.seed(454094511)
absorp_fit <- suppressWarnings(
    phylolm(fa_c ~ taxon, data = absorp_df, phy = absorp_tr, 
            model = 'lambda', boot = 2000)
)
summary(absorp_fit)


#' 
#' # Session info
#' 
#' This outlines the package versions I used for these analyses.
#' 

devtools::session_info()

