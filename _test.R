suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(phylolm)
    library(ape)
})
devtools::load_all('corphyloCpp')

invisible(sapply(list.files('R', '*.R', full.names = TRUE), source))



clear_df <- get_df('clear')
clear_se_df <- get_df('clear', .stat = 'se')  # <-- contains standard errors
clear_tr <- get_tr('clear')

Xmat <- cp_mat(clear_df, c('log_sef', 'log_clear', 'log_enterocyte_density'))
Xmat <- Xmat[!is.na(rowSums(Xmat)),]
MEmat <- cp_mat(clear_se_df, c('log_sef', 'log_clear', 'log_enterocyte_density'))
MEmat <- MEmat[!is.na(rowSums(MEmat)),]

X = Xmat; phy = filter_tr(clear_tr, rownames(Xmat)); SeM = MEmat
U = list(); REML = TRUE; constrain.d = FALSE; reltol = 10^-6; 
maxit.NM = 1000; maxit.SA = 1000; temp.SA = 1; tmax.SA = 1; 
verbose = FALSE



Rcpp::sourceCpp(code = 
'#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat corphylo_LL2(const arma::vec& par, const arma::mat& XX, const arma::mat& UU, 
                   const arma::mat& MM, const arma::mat& tau, const arma::mat& Vphy) {
    
    uint n = Vphy.n_rows;
    uint p = XX.n_rows / n;
    arma::vec L_elements = par(arma::span(0, (p + p * (p - 1)/2) - 1));
    arma::mat L(p, p, arma::fill::zeros);
    for (uint i = 0, j = 0, k = p - 1; i < p; i++) {
        L(arma::span(i, p-1), i) = L_elements(arma::span(j, k));
        j = k + 1;
        k += (p - i - 1);
    }
    arma::mat RR = L + L.t();
    RR.diag().ones();
    return RR;
}
')

x = matrix(1:5, 5, 1)
x
test(x, 1, 3)

corphylo_LL2(par, XX, UU, MM, tau, Vphy)


corphylo.LL <- function(par, XX, UU, MM, tau, Vphy, REML, 
                        constrain.d, verbose) {
    n <- nrow(X)
    p <- ncol(X)
    L.elements <- par[1:(p + p * (p - 1)/2)]
    L <- matrix(0, nrow = p, ncol = p)
    L[lower.tri(L, diag = T)] <- L.elements
    R <- t(L) %*% L
    if (constrain.d == TRUE) {
        logit.d <- par[(p + p * (p - 1)/2 + 1):length(par)]
        if (max(abs(logit.d)) > 10) 
            return(10^10)
        d <- 1/(1 + exp(-logit.d))
    } else {
        d <- par[(p + p * (p - 1)/2 + 1):length(par)]
        if (max(d) > 10) 
            return(10^10)
    }
    C <- matrix(0, nrow = p * n, ncol = p * n)
    for (i in 1:p) for (j in 1:p) {
        Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
        C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * 
                                                            n)] <- R[i, j] * Cd
    }
    V <- C + diag(as.numeric(MM))
    if (is.nan(rcond(V)) || rcond(V) < 10^-10) 
        return(10^10)
    iV <- solve(V)
    denom <- t(UU) %*% iV %*% UU
    if (is.nan(rcond(denom)) || rcond(denom) < 10^-10) 
        return(10^10)
    num <- t(UU) %*% iV %*% XX
    B <- solve(denom, num)
    B <- as.matrix(B)
    H <- XX - UU %*% B
    logdetV <- -determinant(iV)$modulus[1]
    if (is.infinite(logdetV)) 
        return(10^10)
    if (REML == TRUE) {
        LL <- 0.5 * (logdetV + determinant(t(UU) %*% iV %*% 
                                               UU)$modulus[1] + t(H) %*% iV %*% H)
    }
    else {
        LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
    }
    if (verbose == T) 
        show(c(as.numeric(LL), par))
    return(as.numeric(LL))
}
