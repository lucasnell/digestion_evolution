

# ==============================

# Script to do parametric bootstrapping on ape::corphylo

# ==============================

# 
# Rcpp::sourceCpp(code = 
# '#include <RcppArmadillo.h>
# #include <numeric>
# #include <cmath>
# 
# using namespace Rcpp;
# 
# //[[Rcpp::depends(RcppArmadillo)]]
# //[[Rcpp::plugins(cpp11)]]
# 
# //[[Rcpp::export]]
# arma::mat test(arma::mat M) {
#     return arma::inv(M);
# }
# ')
# 


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
source("corphylo.R")

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

# X = Xmat
# phy = clear_ed_tr
# SeM = MEmat
# U = Umat
# REML = TRUE
# method = "Nelder-Mead"
# constrain.d = FALSE
# reltol = 10^-6
# maxit.NM = 1000
# maxit.SA = 1000
# temp.SA = 1
# tmax.SA = 1
# verbose = FALSE
# 
# Rcpp::sourceCpp('corphylo.cpp')
# source('corphylo.R')
# 
# phy <- reorder(phy, "postorder")
# n <- length(phy$tip.label)
# if (dim(X)[1] != n) 
#     stop("Number of rows of the data matrix does not match the length of the tree.")
# if (is.null(rownames(X))) {
#     warning("No tip labels on X; order assumed to be the same as in the tree.\n")
#     data.names = phy$tip.label
# } else data.names = rownames(X)
# order <- match(data.names, phy$tip.label)
# if (sum(is.na(order)) > 0) {
#     warning("Data names do not match with the tip labels.\n")
#     rownames(X) <- data.names
# } else {
#     temp <- X
#     rownames(X) <- phy$tip.label
#     X[order, ] <- temp[1:nrow(temp), ]
# }
# p <- dim(X)[2]
# if (!is.null(SeM)) {
#     if (dim(SeM)[1] != n) 
#         stop(
#             "Number of rows of the SeM matrix does not match the length of the tree.")
#     if (is.null(rownames(SeM))) {
#         warning(
#             "No tip labels on SeM; order assumed to be the same as in the tree.\n")
#         data.names = phy$tip.label
#     }
#     else data.names = rownames(SeM)
#     order <- match(data.names, phy$tip.label)
#     if (sum(is.na(order)) > 0) {
#         warning("SeM names do not match with the tip labels.\n")
#         rownames(SeM) <- data.names
#     }
#     else {
#         temp <- SeM
#         rownames(SeM) <- phy$tip.label
#         SeM[order, ] <- temp[1:nrow(temp), ]
#     }
# } else {
#     SeM <- matrix(0, nrow = n, ncol = p)
# }
# if (length(U) > 0) {
#     if (length(U) != p) 
#         stop(
#             "Number of elements of list U does not match the number of columns in X.")
#     for (i in 1:p) {
#         if (!is.null(U[[i]])) {
#             if (dim(U[[i]])[1] != n) 
#                 stop("Number of rows of an element of U does not match the tree.")
#             if (is.null(rownames(U[[i]]))) {
#                 warning(
#                     "No tip labels on U; order assumed to be the same as in the tree.\n")
#                 data.names = phy$tip.label
#             }
#             else data.names = rownames(U[[i]])
#             order <- match(data.names, phy$tip.label)
#             if (sum(is.na(order)) > 0) {
#                 warning("U names do not match with the tip labels.\n")
#                 rownames(U[[i]]) <- data.names
#             }
#             else {
#                 temp <- U[[i]]
#                 rownames(U[[i]]) <- phy$tip.label
#                 U[[i]][order, ] <- temp[1:nrow(temp), ]
#             }
#         }
#         else {
#             U[[i]] <- matrix(0, nrow = n, ncol = 1)
#             rownames(U[[i]]) <- phy$tip.label
#         }
#     }
# }
# Xs <- X
# for (i in 1:p) Xs[, i] <- (X[, i] - mean(X[, i]))/sd(X[, 
#                                                        i])
# if (!is.null(SeM)) {
#     SeMs <- SeM
#     for (i in 1:p) SeMs[, i] <- SeM[, i]/sd(X[, i])
# }
# if (length(U) > 0) {
#     Us <- U
#     for (i in 1:p) for (j in 1:ncol(U[[i]])) {
#         if (sd(U[[i]][, j]) > 0) {
#             Us[[i]][, j] <- (U[[i]][, j] - mean(U[[i]][, 
#                                                        j]))/sd(U[[i]][, j])
#         }
#         else {
#             Us[[i]][, j] <- U[[i]][, j] - mean(U[[i]][, j])
#         }
#     }
# }
# Vphy <- vcv(phy)
# Vphy <- Vphy/max(Vphy)
# Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
# XX <- matrix(as.matrix(Xs), ncol = 1)
# MM <- matrix(as.matrix(SeMs^2), ncol = 1)
# UU <- kronecker(diag(p), matrix(1, nrow = n, ncol = 1))
# if (length(U) > 0) {
#     zeros <- 0 * (1:p)
#     for (i in 1:p) {
#         dd <- zeros
#         dd[i] <- 1
#         u <- kronecker(dd, as.matrix(Us[[i]]))
#         for (j in 1:dim(u)[2]) if (sd(u[, j]) > 0) 
#             UU <- cbind(UU, u[, j])
#     }
# }
# if (length(U) > 0) {
#     eps <- matrix(nrow = n, ncol = p)
#     for (i in 1:p) {
#         if (ncol(U[[i]]) > 0) {
#             u <- as.matrix(Us[[i]])
#             z <- lm(Xs[, i] ~ u)
#             eps[, i] <- resid(z)
#         }
#         else {
#             eps[, i] <- Xs[, i] - mean(Xs[, i])
#         }
#     }
#     L <- t(chol(cov(eps)))
# } else {
#     L <- t(chol(cov(Xs)))
# }
# L.elements <- L[lower.tri(L, diag = T)]
# par <- c(L.elements, array(0.5, dim = c(1, p)))
# tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
# 
# 
# 
# corphylo_LL(par, XX, UU, MM, tau, Vphy, REML, constrain_d = constrain.d)
# corphylo.LL(par, XX, UU, MM, tau, Vphy, REML, constrain.d, verbose)
# 
# 
# opt <- optim(fn = corphylo.LL, par = par, XX = XX, UU = UU, 
#              MM = MM, tau = tau, Vphy = Vphy, REML = REML, verbose = verbose, 
#              constrain.d = constrain.d, method = "Nelder-Mead", 
#              control = list(maxit = maxit.NM, reltol = reltol))
# 
# opt2 <- optim(fn = corphylo_LL, par = par, XX = XX, UU = UU, 
#              MM = MM, tau = tau, Vphy = Vphy, REML = REML, 
#              constrain_d = constrain.d, method = "Nelder-Mead", 
#              control = list(maxit = maxit.NM, reltol = reltol))
# 
# rm(opt, opt2)
# 


cp <- corphylo(Xmat, phy = clear_ed_tr, SeM = MEmat, U = Umat, method = "Nelder-Mead")

# Bootstrap Pearson r from corphylo object
boot_corrs <- function(cp) {
    # Set up parameter values for simulating data
    p <- length(cp$d)
    n <- nrow(cp$XX) / p
    
    # 
    has_U <- ncol(cp$UU) > 2
    if (has_U) {
        U_inds <- strsplit(gsub('B', '', rownames(cp$B)), "\\.") %>% 
            lapply(as.numeric) %>% 
            do.call(what = rbind)
        U_inds <- U_inds[U_inds[,2] > 0,]
        if (is.null(attributes(U_inds)$dim)) attributes(U_inds)$dim <- c(1, 2)
    }
    
    ## Perform a Cholesky decomposition of Vphy. This is used to generate
    ## phylogenetic signal: a vector of independent normal random variables,
    ## when multiplied by the transpose of the Cholesky deposition of Vphy will
    ## have covariance matrix equal to Vphy.
    iD <- t(chol(cp$V))
    
    # Perform Nrep simulations and collect the results
    Nrep <- 50L
    corrs <- numeric(Nrep)
    
    for (rep in 1L:Nrep) {
        
        rnd <- cbind(rnorm(n * p))
        X <- matrix(iD %*% rnd, nrow = n, ncol = p)
        
        if (has_U) {
            UUcol <- 3
            for (i in U_inds[,1]) {
                mult <- cp$UU[,UUcol]
                mult <- mult[mult != 0]
                X[,i] <- X[,i] + as.numeric(cp$B[sprintf('B%i.1', i),1]) * mult
                UUcol <- UUcol + 1
            }
        }
        
        rownames(X) <- rownames(cp$Vphy)
        
        # z <- corphylo(X = X, SeM = SeM, U = U, phy = phy, method = "Nelder-Mead")
        # z <- corphylo_cpp(X = X, SeM = SeM, U = U, phy = phy, method = "Nelder-Mead")
        z <- cp
        
        corrs[rep] <- z$cor.matrix[1, 2]
        
    }
    return(corrs)
}


hist(corrs); abline(v = cp$cor.matrix[1,2], lty = 2)

mean(corrs); cp$cor.matrix[1,2]


