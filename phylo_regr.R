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
#' # Summary functions
#' 
#' Functions to get p values and 95 CIs, respectively, from a bootstrapped `phylolm`
#' model object.
#' 
pval <- function(model, parameter = 'taxonBat') {
    2 * min(mean(model$bootstrap[,parameter] < 0), 
            mean(model$bootstrap[,parameter] > 0))
}
ci <- function(model, parameter = 'taxonBat') model$bootconfint95[,parameter]
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
#' # `SEF ~ Diet`
#' 
#' 
set.seed(581120)
diet_mod <- phylolm(sef ~ diet, data = diet_df, phy = diet_tr, 
                    model = 'lambda', boot = 2000)
summary(diet_mod)
pval(diet_mod, 'dietOmnivorous'); pval(diet_mod, 'dietProtein')
#' 
#' 
#' 
#' # `Absorption ~ Taxon`
#' 
#' 
#+ absorp_taxon
set.seed(454094511)
absorp_fit <- suppressWarnings(
    phylolm(fa_c ~ taxon, data = absorp_df, phy = absorp_tr, 
            model = 'lambda', boot = 2000)
)
# summary(absorp_fit)
pval(absorp_fit, 'taxonBat')
ci(absorp_fit, 'taxonBat')
#' 
#' 
#' 
#' 
#' # `Morphometrics ~ Taxon`
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
#' These are the column names for the above parameters:
#+ sp_analyses_cols
spp_ys <- c("int_length_mass", "nsa_mass", "vill_area_mass", "log_total_enterocytes")
#' 
#' 
#' 
#' The actual analyses (takes ~6.5 min, which is why it's commented out):
#' 
#+ sp_analyses
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
#' # `Morphometrics ~ Taxon`, separately by segment
#' 
#' 
#' (Segment = proximal, medial, or distal)
#' 
#' List of `Y`s:
#' 
#' - Intestinal diameter (log-transformed, body mass as covariate)
#' - Villus height (body mass as covariate)
#' - Villus width
#' - Crypt width
#' - Surface enlargement factor (SEF)
#' - Enterocyte diameter
#' - Enterocytes per cm^2 NSA (log-transformed)
#' 
#' 
#' 
#' 
#' Constructing a vector of column names for these parameters using the previously 
#' constructed `pos_measures` vector.
#+ pos_measures
pos_ys <- pos_measures[pos_measures != 'mass']
pos_ys[pos_ys == 'intestinal_diameter'] <- 'log_intestinal_diameter'
pos_ys[pos_ys == 'enterocyte_density'] <- 'log_enterocyte_density'
#' 
#' 
#' 
#' The actual analyses (takes ~21.7 min, which is why it's commented out):
#' 
#+ pos_analyses
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
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # `log(Clearance) ~ log(SEF)`
#' 
#' Clearance = "paracellular probe L-arabinose clearance"
#' 
#' From the original manuscript: 
#' > ... we used reduced major axis regression (model II regression)... because both 
#' > variables [X and Y] were subject to error
#' 
#' Instead of an RMA regression, I'll be using `ape::corphylo` to estimate Pearson 
#' correlation coefficients.
#' 
#+ clear_sef


hist(clear_df$clear)

X <- clear_df$log_sef
Y <- clear_df$log_clear
names(X) <- rownames(clear_df)
names(Y) <- rownames(clear_df)
# clear_rma <- phyl.RMA(X, Y, clear_tr, method = 'lambda')
# clear_rma
# plot(clear_rma)

?corphylo
Xmat <- cbind(clear_df$log_sef, clear_df$log_clear)
rownames(Xmat) <- rownames(clear_df)

cp <- corphylo(Xmat, phy = clear_tr, method = "Nelder-Mead")
cp

seed = 2012700501


n <- nrow(Xmat)
phy <- clear_tr
d <- cp$d
p <- length(d)
R <- cp$R
B2 <- cp$B[2,1]
Vphy <- cp$Vphy
MM <- cp$MM
V <- cp$V
iD <- t(chol(V))



set.seed(seed)
star <- stree(n)
star$edge.length <- array(1, dim = c(n, 1))
star$tip.label <- phy$tip.label

cp$U

cp$XX; iD


# Perform Nrep simulations and collect the results
Nrep <- 100
cor.list <- matrix(0, nrow = Nrep, ncol = 1)
cor.noP.list <- matrix(0, nrow = Nrep, ncol = 1)
d.list <- matrix(0, nrow = Nrep, ncol = 2)
B.list <- matrix(0, nrow = Nrep, ncol = 3)
B.noP.list <- matrix(0, nrow = Nrep, ncol = 3)
for (rep in 1:Nrep) {
    XX <- iD
    X <- matrix(XX, nrow = n, ncol = 2)
    rownames(X) <- phy$tip.label
    
    # U <- list(NULL, matrix(rnorm(n, mean = 2, sd = 10), nrow = n, ncol = 1))
    # rownames(U[[2]]) <- phy$tip.label
    # colnames(U[[2]]) <- "V1"
    # X[,2] <- X[,2] + B2[1] * U[[2]][,1] - B2[1] * mean(U[[2]][,1])
    
    z <- corphylo(X = X, phy = phy, method = "Nelder-Mead")
    # z.noM <- corphylo(X = X, U = U, phy = phy, method = "Nelder-Mead")
    z.noP <- corphylo(X = X, phy = star, method = "Nelder-Mead")
    
    cor.list[rep] <- z$cor.matrix[1, 2]
    # cor.noM.list[rep] <- z.noM$cor.matrix[1, 2]
    cor.noP.list[rep] <- z.noP$cor.matrix[1, 2]
    # cor.noMP.list[rep] <- cor(cbind(lm(X[,1] ~ 1)$residuals, lm(X[,2] ~ U[[2]])$residuals))[1,2]
    
    d.list[rep, ] <- z$d
    d.noM.list[rep, ] <- z.noM$d
    
    B.list[rep, ] <- z$B
    B.noM.list[rep, ] <- z.noM$B
    B.noP.list[rep, ] <- z.noP$B
    
    show(c(rep, z$convcode, z$cor.matrix[1, 2], z$d))
}




#' 
#' 
#' 
#' # Session info
#' 
#' This outlines the package versions I used for these analyses.
#' 

devtools::session_info()

