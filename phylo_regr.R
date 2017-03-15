# Regression with phylogenetic covariance

source('tidy_csv.R')


library(magrittr)
library(phylolm)
library(ape)
library(nlme)
# library(Rphylopars)
# library(geiger)
library(ggplot2)
library(ggtree)



# package "Rphylopars" for version 3 of measurement error estimation
# "geiger" by Luke Harmon; function fitContinuous



tr <- read.tree('tree.nwk')
tr$tip.label <- gsub('_', ' ', tr$tip.label)
tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% (morph_df$species %>% unique)])
# plot(tr)
tr


#' 
#' All measures found in `morph_df`:
#' 
#' - villus height
#' - villus width
#' - crypt width
#' - enterocytes per surface area
#' - sef
#' - intestinal diameter
#' - intestinal length
#' - nsa
#' - body mass
#' - total villa surface area
#' - enterocite width
#' 
#' 




# Return 'wide' version of morph_df, given a character vector of measures
prep_df <- function(measures, df = morph_df, tree = tr, by_sp = TRUE, trans_fun = 'log'){
    
    meas_list <- lapply(measures, function(s) sprintf('`%s`', s))
    meas_clean <- gsub(' ', '_', measures)
    trans_meas <- paste(meas_clean, trans_fun, sep = '_')
    
    new_df <- df %>% 
        spread(measure, value) %>% 
        select_(.dots = append(list('diet', 'taxa', 'species', 'individual'), 
                               meas_list)) %>% 
        # Removing all rows with all NAs in measures columns
        filter(Reduce(`+`, lapply(.[,measures], is.na)) < length(measures)) %>% 
        # Replacing spaces in measures-column names with underscores
        rename_(.dots = setNames(as.list(sprintf('`%s`', measures)), meas_clean)) %>% 
        # Now taking mean by sample
        group_by(diet, taxa, species, individual) %>% 
        summarize_all(mean, na.rm = TRUE) %>% 
        ungroup %>% 
        # Doing the transformation now, before taking mean if aggregating by species
        mutate_(.dots = setNames(as.list(sprintf('%s(%s)', trans_fun, meas_clean)), 
                                 trans_meas))
    
    if (by_sp) {
        new_df <- new_df %>% 
            group_by(diet, taxa, species) %>% 
            summarize_at(.cols = c(trans_meas, meas_clean), mean) %>% 
            ungroup %>% 
            arrange(taxa, diet, species) %>% 
            # To change row names, it can't be a tibble
            as.data.frame
        rownames(new_df) <- new_df$species
        # Now sort by the phylogenetic tree's tip labels
        new_df <- new_df[match(tree$tip.label,rownames(new_df)),]
    } else {
        new_df <- new_df %>% select(species, everything()) %>% 
            as.data.frame
    }
    
    return(new_df)
}




sp_df <- prep_df(measures = c('nsa', 'sef', 'body mass'))
indiv_df <- prep_df(c('nsa', 'sef', 'body mass'), by_sp = FALSE) %>% 
    select(species, taxa, nsa_log, body_mass_log)


# Phylogenetic tree with body mass as species name color
gg_tr <- ggtree(tr)
gg_tr$data$x <- gg_tr$data$x - max(gg_tr$data$x)

tr_p <- gg_tr %<+% {sp_df %>% select(species, everything())} +
    theme_tree2(axis.title.x = element_text(size = 14)) +
    geom_tiplab(aes(color = body_mass_log), size = 2, fontface = 'bold.italic') +
    scale_x_continuous('Time (mya)', limits = c(-100, 10),
                       breaks = seq(-100, 0, 25), labels = seq(100, 0, -25)) +
    scale_color_gradient(low = 'yellow', high = 'black') +
    theme(legend.position = c(0.25, 0.75), 
          legend.background = element_rect(color = NA, fill = NA))
tr_p






# Takes ~X minutes
t0 <- Sys.time()
set.seed(352)
nsa_fits <- lapply(c('lambda', 'OUfixedRoot'), 
                  function(m) {
                      phylolm(nsa_log ~ body_mass_log + taxa, data = sp_df, phy = tr,
                              model = m, boot = 2000)})
sef_fits <- lapply(c('lambda', 'OUfixedRoot'), 
                   function(m) {
                       phylolm(sef_log ~ body_mass_log + taxa, data = sp_df, phy = tr,
                               model = m, boot = 2000)})
t1 <- Sys.time()
t1 - t0



summary(nsa_fits[[1]])

# P-value based on bootstrap replicates for coefficient of taxaRodent != 0
mean(nsa_fits[[1]]$bootstrap[,'taxaRodent'] < 0) * 2





# ou_fit2 <- phylolm(nsa_log ~ body_mass_log * taxa, data = sp_df, phy = tr,
#                    model = 'OUfixedRoot')
# summary(ou_fit2)






# ================================================================
# ================================================================

# Trying regression separately by taxon

# ================================================================
# ================================================================


# ====================
# With phylolm
# ====================

rod_mod <- phylolm(nsa_log ~ body_mass_log, 
                   data = sp_df[sp_df$taxa == 'Rodent',], 
                   phy = drop.tip(tr, tip = sp_df$species[sp_df$taxa != 'Rodent']),
                   model = 'lambda', boot = 100)
bat_mod <- phylolm(nsa_log ~ body_mass_log, 
                   data = sp_df[sp_df$taxa == 'Bat',], 
                   phy = drop.tip(tr, tip = sp_df$species[sp_df$taxa != 'Bat']),
                   model = 'lambda', boot = 100)
summary(rod_mod)
summary(bat_mod)


# Calculating confidence intervals
mm <- model.matrix(~ 1 + body_mass_log, data = sp_df[sp_df$taxa == 'Rodent',])
vc <- as.matrix(vcov(rod_mod))
ci_w <- apply(mm, 1, function(z) sqrt(t(z) %*% vc %*% z)) * 1.96
rod_ci <- data_frame(
    body_mass_log = sp_df$body_mass_log[sp_df$taxa == 'Rodent'],
    nsa_log = predict(rod_mod), 
    hi = predict(rod_mod) + ci_w, 
    lo = predict(rod_mod) - ci_w, 
    taxa = 'Rodent')

mm <- model.matrix(~ 1 + body_mass_log, data = sp_df[sp_df$taxa == 'Bat',])
vc <- as.matrix(vcov(bat_mod))
ci_w <- apply(mm, 1, function(z) sqrt(t(z) %*% vc %*% z)) * 1.96
bat_ci <- data_frame(
    body_mass_log = sp_df$body_mass_log[sp_df$taxa == 'Bat'],
    nsa_log = predict(bat_mod), 
    hi = predict(bat_mod) + ci_w, 
    lo = predict(bat_mod) - ci_w,
    taxa = 'Bat')

# Plotting model predictions with CI, plus raw data
bind_rows(rod_ci, bat_ci) %>% 
    ggplot(aes(body_mass_log, nsa_log, color = taxa)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = taxa), color = NA, 
                fill = 'gray20', alpha = 0.2) +
    geom_point(data = sp_df, shape = 1) +
    geom_line()




# ====================
# With gls
# ====================

rod_gls <- gls(nsa_log ~ body_mass_log, 
               data = sp_df[sp_df$taxa == 'Rodent',], 
               method = 'REML', 
               correlation = 
                   corPagel(1, 
                            drop.tip(tr, tip = sp_df$species[sp_df$taxa != 'Rodent']), 
                            fixed = FALSE))
bat_gls <- gls(nsa_log ~ body_mass_log, 
               data = sp_df[sp_df$taxa == 'Bat',], 
               method = 'REML', 
               correlation = 
                   corPagel(1, 
                            drop.tip(tr, tip = sp_df$species[sp_df$taxa != 'Bat']), 
                            fixed = FALSE))
summary(rod_gls)
summary(bat_gls)




# Calculating confidence intervals
mm <- model.matrix(~ 1 + body_mass_log, data = sp_df[sp_df$taxa == 'Rodent',])
vc <- as.matrix(vcov(rod_gls))
ci_w <- apply(mm, 1, function(z) sqrt(t(z) %*% vc %*% z)) * 1.96
rod_ci_gls <- data_frame(
    body_mass_log = sp_df$body_mass_log[sp_df$taxa == 'Rodent'],
    nsa_log = predict(rod_gls), 
    hi = predict(rod_gls) + ci_w, 
    lo = predict(rod_gls) - ci_w, 
    taxa = 'Rodent')

mm <- model.matrix(~ 1 + body_mass_log, data = sp_df[sp_df$taxa == 'Bat',])
vc <- as.matrix(vcov(bat_gls))
ci_w <- apply(mm, 1, function(z) sqrt(t(z) %*% vc %*% z)) * 1.96
bat_ci_gls <- data_frame(
    body_mass_log = sp_df$body_mass_log[sp_df$taxa == 'Bat'],
    nsa_log = predict(bat_gls), 
    hi = predict(bat_gls) + ci_w, 
    lo = predict(bat_gls) - ci_w,
    taxa = 'Bat')

# Plotting model predictions with CI, plus raw data
bind_rows(rod_ci_gls, bat_ci_gls) %>% 
    ggplot(aes(body_mass_log, nsa_log, color = taxa)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = taxa), color = NA, 
                fill = 'gray20', alpha = 0.2) +
    geom_point(data = sp_df, shape = 1) +
    geom_line()






# ====================
# With phylopars.lm
# ====================

rod_indiv <- phylopars.lm(nsa_log ~ body_mass_log, 
                          trait_data = indiv_df[indiv_df$taxa == 'Rodent',], 
                          tree = drop.tip(tr, tip = sp_df$species[sp_df$taxa != 'Rodent']),
                          model = 'lambda')

bat_indiv <- phylopars.lm(nsa_log ~ body_mass_log, 
                          trait_data = indiv_df[indiv_df$taxa == 'Bat',], 
                          tree = drop.tip(tr, tip = sp_df$species[sp_df$taxa != 'Bat']),
                          model = 'lambda')

summary(rod_indiv)
rod_indiv$PPE$model$lambda
summary(bat_indiv)
bat_indiv$PPE$model$lambda


