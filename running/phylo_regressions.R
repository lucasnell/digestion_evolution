## ----load_packages-------------------------------------------------------
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(phylolm)
    library(ape)
})

## ----source_R------------------------------------------------------------
invisible(sapply(list.files('R', '*.R', full.names = TRUE), source))

## ----diet_data-----------------------------------------------------------
spp_df <- get_df('spp')
tr <- get_tr('spp')

## ----sef_diet, eval = FALSE----------------------------------------------
## set.seed(581120)
## diet_fit <- phylolm(sef ~ diet, data = spp_df, phy = tr,
##                     model = 'lambda', boot = 2000)

## ----diet_fit_summ-------------------------------------------------------
summary(diet_fit)

## ----absorp_data---------------------------------------------------------
absorp_df <- get_df('absorp')
absorp_tr <- get_tr('absorp')

## ----absorp_taxon, eval = FALSE------------------------------------------
## set.seed(454094511)
## absorp_fit <- suppressWarnings(  # gives warning about lambda being very low
##     phylolm(absorp ~ taxon, data = absorp_df, phy = absorp_tr,
##             model = 'lambda', boot = 2000)
## )

## ----absorp_fit_summ-----------------------------------------------------
summary(absorp_fit)

## ----sp_analyses_cols----------------------------------------------------
spp_ys <- c("int_length_mass", "nsa_mass", "vill_area_mass", "log_total_enterocytes")

## ----sp_analyses, eval = FALSE-------------------------------------------
## set.seed(88754829)
## spp_fits <- lapply(
##     spp_ys,
##     function(y) {
##         f <- paste(y, ' ~ taxon',
##                    ifelse(grepl('total_enterocytes', y), '+ log_mass', ''))
##         suppressWarnings(
##             do.call("phylolm", list(as.formula(f), data = as.name("spp_df"),
##                                     phy = as.name("tr"), model = 'lambda',
##                                     boot = 2000))
##         )
##     })
## names(spp_fits) <- spp_ys
## readr::write_rds(spp_fits, 'output/models_spp.rds')

## ----pos_ys--------------------------------------------------------------
pos_ys <- c('log_intestinal_diameter', 'villus_height', 'villus_width', 
            'crypt_width', 'sef', 'enterocyte_diameter', 'log_enterocyte_density')
seg_types <- c('prox', 'med', 'dist')

## ----pos_analyses, eval = FALSE------------------------------------------
## set.seed(25413535)
## pos_fits <- lapply(
##     seg_types,
##     function(pos) {
##         # Assigning to obj named <pos>_df so that the call identifies the position
##         assign(paste0(pos, '_df'), get_df('pos', .pos = pos))
##         lapply(
##             pos_ys,
##             function(y) {
##                 f <- paste(y, ' ~ taxon',
##                            ifelse(grepl('intestinal_diameter|villus_height', y),
##                                   '+ log_mass', ''))
##                 arg_list <- list(
##                     as.formula(f),
##                     data = as.name(paste0(pos, "_df")),
##                     phy = as.name("tr"), model = 'lambda',
##                     boot = 2000)
##                 # These models don't find the peak likelihood unless specifying a
##                 # starting value of 0.1.
##                 if ((y == "log_enterocyte_density" & pos == "med") |
##                     (y == "crypt_width" & pos == "prox")) {
##                     arg_list <- c(arg_list, starting.value = 0.1)
##                 }
##                 # Now call phylolm
##                 suppressWarnings(do.call("phylolm", arg_list))
##             })
##     })
## names(pos_fits) <- seg_types
## for (i in 1:length(pos_fits)) names(pos_fits[[i]]) <- pos_ys; rm(i)
## readr::write_rds(pos_fits, 'output/models_pos.rds')

## ----clear_sef-----------------------------------------------------------

clear_df <- get_df('clear')
clear_se_df <- get_df('clear', .stat = 'se')  # <-- contains standard errors
clear_tr <- get_tr('clear')

Xmat <- cbind(clear_df$log_sef, clear_df$log_clear)
rownames(Xmat) <- rownames(clear_df)

MEmat <- cbind(clear_se_df$log_sef, clear_se_df$log_clear)
rownames(MEmat) <- clear_se_df$species

clear_cor <- corp(Xmat, phy = clear_tr, SeM = MEmat)

# Correlation with 95% CI
clear_cor['r',]

