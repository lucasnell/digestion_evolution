

# ==============================

# Script to do parametric bootstrapping on ape::corphylo

# ==============================


suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(phylolm)
    library(ape)
})
invisible(sapply(list.files('R', '*.R', full.names = TRUE)[list.files('R', '*.R', full.names = TRUE) != "R/bootstrap_corphylo.R"], source))
source(".Rprofile")

devtools::load_all('corphyloCpp')


clear_df <- get_df('clear')
clear_se_df <- get_df('clear', .stat = 'se')  # <-- contains standard errors
clear_tr <- get_tr('clear')



Xmat <- cbind(clear_df$log_enterocyte_density, clear_df$log_clear)
rownames(Xmat) <- rownames(clear_df)
Xmat <- Xmat[!is.na(rowSums(Xmat)),]

MEmat <- cbind(clear_se_df$log_enterocyte_density, clear_se_df$log_clear)
rownames(MEmat) <- clear_se_df$species
MEmat <- MEmat[!is.na(rowSums(MEmat)),]

Umat <- list(NULL, cbind(clear_df$log_mass[!is.na(clear_df$log_mass)]))
rownames(Umat[[2]]) <- rownames(Xmat)

clear_ed_tr <- ape::drop.tip(
    clear_tr, 
    tip = clear_tr$tip.label[!clear_tr$tip.label %in% rownames(Xmat)]
)




cp <- corphylo_cpp(Xmat, phy = clear_ed_tr, SeM = MEmat, U = Umat, method = "Nelder-Mead")


system.time(corrs <- boot_r(cp, 100L, seed = 1, n_cores = 4))

# p value
2 * min(mean(corrs > 0), mean(corrs < 0))

hist(corrs); abline(v = cp$cor.matrix[1,2], lty = 2)

mean(corrs); cp$cor.matrix[1,2]


