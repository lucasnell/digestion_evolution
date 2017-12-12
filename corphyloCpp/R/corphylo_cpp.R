


# Assembling parameter names for SE output

get_par_names <- function(p, grepl_str = NULL) {
    L <- sum(1:p) + p
    par_names <- character(L)
    par_names[(L - p + 1):L] <- sprintf('d%i', 1:p)
    par_names[c(1, sapply(p:2, function(x) 1 + sum(p:x)))] <- sprintf('sig%i', 1:p)
    par_names[par_names == ""] <- unlist(lapply(1:(p-1), 
                                                function(i) {
                                                    sapply((i+1):p, 
                                                           function(j) {
                                                               sprintf('r%i%i', i, j)
                                                           })
                                                }))
    if (!is.null(grepl_str)) {
        par_names <- par_names[grepl(grepl_str, par_names)]
    }
    return(par_names)
}




# Below, `# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
# delimit areas that are edited from the original corphylo function


#' C++ version of ape::corphylo
#' 
#' For people in a hurry
#'
#' @inheritParams ape::corphylo
#' @param boot Number of bootstrap replicates. Defaults to 0.
#' @param boot_out Function to retrieve necessary info from corphylo object for each
#'     bootstrap replicate. If defining your own function for this, make sure that 
#'     the output is either a single number or a matrix row.
#'     Defaults to \code{NULL}, which retrieves the correlation(s) and the 
#'     \code{d}-values for the OU process (i.e., the measures of phylogenetic signal).
#' @param n_cores Number of cores to use. Defaults to 1.
#'
#' @return
#' 
#' @seealso \code{\link[ape]{corphylo}} \code{\link{boot_r}}
#' 
#' @useDynLib corphyloCpp
#' 
#' @export
#'
corphylo_cpp <- function(X, U = list(), SeM = NULL, phy = NULL, REML = TRUE, 
          method = c("Nelder-Mead", "SANN"), constrain.d = FALSE, reltol = 10^-6, 
          maxit.NM = 1000, maxit.SA = 1000, temp.SA = 1, tmax.SA = 1, 
          verbose = FALSE,
          boot = 0, boot_out = NULL, n_cores = 1, max_iter = 100) {
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Added this bc it returns warning if method if left as default
    method <- match.arg(method)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (!inherits(phy, "phylo")) stop("Object \"phy\" is not of class \"phylo\".")
    if (is.null(phy$edge.length)) stop("The tree has no branch lengths.")
    if (is.null(phy$tip.label)) stop("The tree has no tip labels.")
    
    phy <- reorder(phy, "postorder")
    n <- length(phy$tip.label)
    
    if (dim(X)[1] != n) {
        stop("Number of rows of the data matrix does not match the length of the tree.")
    }
    if (is.null(rownames(X))) {
        warning("No tip labels on X; order assumed to be the same as in the tree.\n")
        data.names = phy$tip.label
    } else data.names = rownames(X)
    order <- match(data.names, phy$tip.label)
    if (sum(is.na(order)) > 0) {
        warning("Data names do not match with the tip labels.\n")
        rownames(X) <- data.names
    } else {
        temp <- X
        rownames(X) <- phy$tip.label
        X[order, ] <- temp[1:nrow(temp), ]
    }
    p <- dim(X)[2]
    if (!is.null(SeM)) {
        if (dim(SeM)[1] != n) {
            stop(
                "Number of rows of the SeM matrix does not match the length of the tree.")
        }
        if (is.null(rownames(SeM))) {
            warning(
                "No tip labels on SeM; order assumed to be the same as in the tree.\n")
            data.names = phy$tip.label
        } else data.names = rownames(SeM)
        order <- match(data.names, phy$tip.label)
        if (sum(is.na(order)) > 0) {
            warning("SeM names do not match with the tip labels.\n")
            rownames(SeM) <- data.names
        } else {
            temp <- SeM
            rownames(SeM) <- phy$tip.label
            SeM[order, ] <- temp[1:nrow(temp), ]
        }
    } else {
        SeM <- matrix(0, nrow = n, ncol = p)
    }
    if (length(U) > 0) {
        if (length(U) != p) {
            stop(
                "Number of elements of list U does not match the number of columns in X.")
        }
        for (i in 1:p) {
            if (!is.null(U[[i]])) {
                if (dim(U[[i]])[1] != n) 
                    stop("Number of rows of an element of U does not match the tree.")
                if (is.null(rownames(U[[i]]))) {
                    warning(
                        "No tip labels on U; order assumed to be the same as in the tree.\n")
                    data.names = phy$tip.label
                } else data.names = rownames(U[[i]])
                order <- match(data.names, phy$tip.label)
                if (sum(is.na(order)) > 0) {
                    warning("U names do not match with the tip labels.\n")
                    rownames(U[[i]]) <- data.names
                } else {
                    temp <- U[[i]]
                    rownames(U[[i]]) <- phy$tip.label
                    U[[i]][order, ] <- temp[1:nrow(temp), ]
                }
            } else {
                U[[i]] <- matrix(0, nrow = n, ncol = 1)
                rownames(U[[i]]) <- phy$tip.label
            }
        }
    }
    Xs <- X
    for (i in 1:p) Xs[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
    if (!is.null(SeM)) {
        SeMs <- SeM
        for (i in 1:p) SeMs[, i] <- SeM[, i]/sd(X[, i])
    }
    if (length(U) > 0) {
        Us <- U
        for (i in 1:p) {
            for (j in 1:ncol(U[[i]])) {
                if (sd(U[[i]][, j]) > 0) {
                    Us[[i]][, j] <- (U[[i]][, j] - mean(U[[i]][, j])) / sd(U[[i]][, j])
                } else {
                    Us[[i]][, j] <- U[[i]][, j] - mean(U[[i]][, j])
                }
            }
        }
    }
    Vphy <- vcv(phy)
    Vphy <- Vphy/max(Vphy)
    Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
    XX <- matrix(as.matrix(Xs), ncol = 1)
    MM <- matrix(as.matrix(SeMs^2), ncol = 1)
    UU <- kronecker(diag(p), matrix(1, nrow = n, ncol = 1))
    if (length(U) > 0) {
        zeros <- 0 * (1:p)
        for (i in 1:p) {
            dd <- zeros
            dd[i] <- 1
            u <- kronecker(dd, as.matrix(Us[[i]]))
            for (j in 1:dim(u)[2]) if (sd(u[, j]) > 0) 
                UU <- cbind(UU, u[, j])
        }
    }
    if (length(U) > 0) {
        eps <- matrix(nrow = n, ncol = p)
        for (i in 1:p) {
            if (ncol(U[[i]]) > 0) {
                u <- as.matrix(Us[[i]])
                z <- lm(Xs[, i] ~ u)
                eps[, i] <- resid(z)
            } else {
                eps[, i] <- Xs[, i] - mean(Xs[, i])
            }
        }
        L <- t(chol(cov(eps)))
    } else {
        L <- t(chol(cov(Xs)))
    }
    L.elements <- L[lower.tri(L, diag = T)]
    par <- c(L.elements, rep(0.5, p))
    tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Below has corphylo.LL changed to corphylo_LL, verbose argument commented out, 
    # and constrain.d changed to constrain_d = constrain.d
    if (method == "Nelder-Mead")  {
        opt <- optim(fn = corphylo_LL, par = par, XX = XX, UU = UU, 
                     MM = MM, tau = tau, Vphy = Vphy, REML = REML, # verbose = verbose, 
                     constrain_d = constrain.d, method = "Nelder-Mead", 
                     control = list(maxit = maxit.NM, reltol = reltol))
    } else if (method == "SANN") {
        opt <- optim(fn = corphylo_LL, par = par, XX = XX, UU = UU, 
                     MM = MM, tau = tau, Vphy = Vphy, REML = REML, # verbose = verbose, 
                     constrain_d = constrain.d, method = "SANN", 
                     control = list(maxit = maxit.SA, temp = temp.SA, tmax = tmax.SA, 
                                    reltol = reltol))
        par <- opt$par
        opt <- optim(fn = corphylo_LL, par = par, XX = XX, UU = UU, 
                     MM = MM, tau = tau, Vphy = Vphy, REML = REML, # verbose = verbose, 
                     constrain_d = constrain.d, method = "Nelder-Mead", 
                     control = list(maxit = maxit.NM, reltol = reltol))
    }
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    par <- Re(opt$par)
    LL <- opt$value
    L.elements <- par[1:(p + p * (p - 1)/2)]
    L <- matrix(0, nrow = p, ncol = p)
    L[lower.tri(L, diag = T)] <- L.elements
    R <- t(L) %*% L
    Rd <- diag(diag(R)^-0.5)
    cor.matrix <- Rd %*% R %*% Rd
    if (constrain.d == TRUE) {
        logit.d <- par[(p + p * (p - 1)/2 + 1):length(par)]
        d <- 1/(1 + exp(-logit.d))
    } else {
        d <- par[(p + p * (p - 1)/2 + 1):length(par)]
    }
    C <- matrix(0, nrow = p * n, ncol = p * n)
    for (i in 1:p) for (j in 1:p) {
        Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
        C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
    }
    V <- C + diag(MM)
    iV <- solve(V)
    denom <- t(UU) %*% iV %*% UU
    num <- t(UU) %*% iV %*% XX
    B <- solve(denom, num)
    B <- as.matrix(B)
    B.cov <- solve(t(UU) %*% iV %*% UU)
    H <- XX - UU %*% B
    counter <- 0
    sd.list <- matrix(0, nrow = dim(UU)[2], ncol = 1)
    for (i in 1:p) {
        counter <- counter + 1
        B[counter] <- B[counter] + mean(X[, i])
        sd.list[counter] <- sd(X[, i])
        if (length(U) > 0) {
            for (j in 1:ncol(U[[i]])) {
                if (sd(U[[i]][, j]) > 0) {
                    counter <- counter + 1
                    B[counter] <- B[counter] * sd(X[, i])/sd(U[[i]][, j])
                    sd.list[counter] <- sd(X[, i])/sd(U[[i]][, j])
                }
            }
        }
    }
    B.cov <- diag(as.numeric(sd.list)) %*% B.cov %*% diag(as.numeric(sd.list))
    B.se <- as.matrix(diag(B.cov))^0.5
    B.zscore <- B/B.se
    B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
    if (length(U) > 0) {
        B.rownames <- NULL
        for (i in 1:p) {
            B.rownames <- c(B.rownames, paste("B", i, ".0", sep = ""))
            if (ncol(U[[i]]) > 0) 
                for (j in 1:ncol(U[[i]])) if (sd(U[[i]][, j]) > 
                                              0) {
                    if (is.null(colnames(U[[i]])[j])) {
                        B.rownames <- c(B.rownames, paste("B", i, 
                                                          ".", j, sep = ""))
                    } else {
                        B.rownames <- c(B.rownames, paste("B", i, 
                                                          ".", colnames(U[[i]])[j], 
                                                          sep = ""))
                    }
                }
        }
    } else {
        B.rownames <- NULL
        for (i in 1:p) {
            B.rownames <- c(B.rownames, paste("B", i, ".0", sep = ""))
        }
    }
    rownames(B) <- B.rownames
    rownames(B.cov) <- B.rownames
    colnames(B.cov) <- B.rownames
    rownames(B.se) <- B.rownames
    rownames(B.zscore) <- B.rownames
    rownames(B.pvalue) <- B.rownames
    if (REML == TRUE) {
        logLik <- -0.5 * ((n * p) - ncol(UU)) * log(2 * pi) + 
            0.5 * determinant(t(XX) %*% XX)$modulus[1] - LL
    } else {
        logLik <- -0.5 * (n * p) * log(2 * pi) - LL
    }
    k <- length(par) + ncol(UU)
    AIC <- -2 * logLik + 2 * k
    BIC <- -2 * logLik + k * (log(n) - log(pi))
    results <- list(cor.matrix = cor.matrix, d = d, B = B, B.se = B.se, 
                    B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue, 
                    logLik = logLik, AIC = AIC, BIC = BIC, REML = REML, 
                    constrain.d = constrain.d, 
                    XX = XX, UU = UU, MM = MM, Vphy = Vphy, R = R, V = V, 
                    C = C, convcode = opt$convergence, niter = opt$counts)
    class(results) <- "corphylo"
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Added bootstrapping to output
    if (boot > 0) {
        results$bootstrap <- boot_corphylo(results, boot, boot_out, n_cores, max_iter)
    } else {
        results$bootstrap <- matrix(NA, 0, 0)
    }
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    return(results)
}







# ================================================================================
# ================================================================================

#               Bootstrapping

# ================================================================================
# ================================================================================


# One bootstrap replicate for a single corphylo object, returning a new fit

one_boot_fit <- function(cp_obj, n, p, iD, U_add, SeM, U, phy) {
    rnd <- cbind(rnorm(n * p))
    X <- iD %*% rnd + U_add
    X <- matrix(X, n, p)
    rownames(X) <- rownames(cp_obj$Vphy)
    z <- corphylo_cpp(X = X, SeM = SeM, U = U, phy = phy, method = "Nelder-Mead")
    return(z)
}




# From r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2
# To make sure Cholesky decomposition works for V matrix

chol_fix <- function(V, max_iter = 100) {
    
    new_V <- V
    
    cholStatus <- try(u <- chol(V), silent = TRUE)
    cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    iter <- 0
    while (cholError & iter < max_iter) {
        
        # replace -ve eigen values with small +ve number
        new_Eig <- eigen(new_V)
        new_Eig2 <- ifelse(new_Eig$values < 0, 0, new_Eig$values)
        
        # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
        # eig vectors
        new_V <- new_Eig$vectors %*% diag(new_Eig2) %*% t(new_Eig$vectors)
        
        # normalize modified matrix eqn 6 from Brissette et al 2007
        new_V <- new_V / sqrt(diag(new_V) %*% t(diag(new_V)))
        
        # try chol again
        cholStatus <- try(u <- chol(new_V), silent = TRUE)
        cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
        
        iter <- iter + 1
    }
    if (cholError) stop("non positive-definite correlation matrix")
    
    return(new_V)
}


# Prepping info for input corphylo

prep_info <- function(cp_obj, max_iter = 100) {
    
    # Set up parameter values for simulating data
    p <- length(cp_obj$d)
    n <- nrow(cp_obj$XX) / p
    
    phy <- vcv2phylo(cp_obj$Vphy)
    
    if (ncol(cp_obj$UU) > 2) {
        U_inds <- strsplit(gsub('B', '', rownames(cp_obj$B)), "\\.") %>% 
            lapply(as.numeric) %>% 
            do.call(what = rbind)
        U_inds <- U_inds[U_inds[,2] > 0,]
        if (is.null(attributes(U_inds)$dim)) attributes(U_inds)$dim <- c(1, 2)
        # Setting up U itself
        U <- as.list(rep(list(NULL), p))
        UUcol <- 3
        for (i in U_inds[,1]) {
            U_i <- cp_obj$UU[,UUcol]
            U_i <- U_i[U_i != 0]
            U[[i]] = cbind(U_i)
            colnames(U[[i]]) <- NULL
            rownames(U[[i]]) <- rownames(cp_obj$Vphy)
            UUcol <- UUcol + 1
        }
        # For adding to X values
        U_add <- cp_obj$UU
        U_add[,1:2] <- 0 # <-- Keeps means at zero
        U_add <- t(t(cp_obj$B) %*% t(U_add))
    } else {
        U <- NULL
        U_inds <- NULL
        U_add <- matrix(0, n*p)
    }
    
    # Measurement error matrix
    SeM <- matrix(cp_obj$MM, n)
    rownames(SeM) <- rownames(cp_obj$Vphy)
    
    ## Perform a Cholesky decomposition of Vphy. This is used to generate
    ## phylogenetic signal: a vector of independent normal random variables,
    ## when multiplied by the transpose of the Cholesky deposition of Vphy will
    ## have covariance matrix equal to Vphy.
    V <- chol_fix(cp_obj$V, max_iter)
    
    iD <- t(chol(V))
    
    return(list(n = n, p = p, iD = iD, U_add = U_add, SeM = SeM, U = U, phy = phy))
    
}





# Default boot_out argument for boot_corphylo; returns correlation and d-values for OU
# process

boot_out_default <- function(cp_obj) {
    x <- cp_obj$cor.matrix
    out <- rbind(c(x[lower.tri(x, FALSE)], cp_obj$d))
    colnames(out) <- c(get_par_names(nrow(x), 'r'), get_par_names(nrow(x), 'd'))
    return(out)
}




# Internal function to do parametric bootstrapping from a corphylo object

boot_corphylo <- function(cp_obj, boot, boot_out = NULL, n_cores = 1, max_iter = 100) {
    
    # This allows output to be reproducible if someone uses set.seed outside this 
    # function, even when doing this in parallel
    seed <- sample.int(2^31-1, 1)
    
    if (is.null(boot_out)) {
        boot_out <- boot_out_default
    }
    one_boot <- function(i, cp_obj, n, p, iD, U_add, SeM, U, phy) {
        z <- one_boot_fit(cp_obj, n, p, iD, U_add, SeM, U, phy)
        return(boot_out(z))
    }
    
    if (!is(cp_obj, 'corphylo')) {
        stop("cp_obj argument must be a 'corphylo' object")
    }
    if (n_cores < 1 | n_cores %% 1 != 0) stop("n_cores must be an integer >= 1")

    # Adding objects from the corphylo object to a list for later using `do.call`
    arg_list <- prep_info(cp_obj, max_iter = 100)
    arg_list <- c(list(X = 1:boot, FUN = one_boot, cp_obj = cp_obj), arg_list)
    
    # Perform boot simulations and collect the results
    if (requireNamespace("parallel", quietly = TRUE) & .Platform$OS.type == 'unix' &
        n_cores > 1) {
        arg_list <- c(arg_list, list(mc.cores = n_cores))
        rng_orig <- RNGkind()
        RNGkind("L'Ecuyer-CMRG")
        set.seed(seed)
        corrs <- do.call(parallel::mclapply, arg_list)
        RNGkind(rng_orig[1])
    } else {
        set.seed(seed)
        corrs <- do.call(lapply, arg_list)
    }
    
    return(do.call(rbind, corrs))
}


# Redefine print.corphylo to include SE info if present
#' @export
print.corphylo <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("Call to corphylo\n\n")
    logLik = x$logLik
    AIC = x$AIC
    BIC = x$BIC
    names(logLik) = "logLik"
    names(AIC) = "AIC"
    names(BIC) = "BIC"
    print(c(logLik, AIC, BIC), digits = digits)
    cat("\ncorrelation matrix:\n")
    rownames(x$cor.matrix) <- 1:dim(x$cor.matrix)[1]
    colnames(x$cor.matrix) <- 1:dim(x$cor.matrix)[1]
    print(x$cor.matrix, digits = digits)
    cat("\nfrom OU process:\n")
    d <- data.frame(d = x$d)
    print(d, digits = digits)
    if (x$constrain.d == TRUE) {
        cat("\nvalues of d constrained to be in [0, 1]\n")
    }
    cat("\ncoefficients:\n")
    coef <- data.frame(Value = x$B, Std.Error = x$B.se, Zscore = x$B.zscore, 
                       Pvalue = x$B.pvalue)
    rownames(coef) <- rownames(x$B)
    printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
    if (x$convcode != 0) {
        cat("\nWarning: convergence in optim() not reached\n")
    }
    if (nrow(x$bootstrap) > 0) {
        cat("\nBootstrapped 95% CI:\n")
        if (!is.null(colnames(x$bootstrap))) {
            for (nn in colnames(x$bootstrap)) {
                cat(sprintf("  %-4s %9.3g [%9.3g %9.3g]\n", 
                            nn, 
                            mean(x$bootstrap[,nn]),
                            as.numeric(quantile(x$bootstrap[,nn], 0.025)),
                            as.numeric(quantile(x$bootstrap[,nn], 0.975))))
            }
        } else {
            for (nn in 1:ncol(x$bootstrap)) {
                cat(sprintf("  %-4s %9.3g [%9.3g %9.3g]\n", 
                            paste0('col', nn), 
                            mean(x$bootstrap[,nn]),
                            as.numeric(quantile(x$bootstrap[,nn], 0.025)),
                            as.numeric(quantile(x$bootstrap[,nn], 0.975))))
            }
        }
    }
    cat("\n")
}
