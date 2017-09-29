# This code defines a function 'corp' that is largely the same as Tony Ives' 
# ape::corphylo function, except that it is parameterized with the trait correlation 
# matrix rather than the trait covariance matrix.
# This made it easier to calculate confidence intervals on the estimated 
# correlation (using Fisher information).

# The output is the estimated trait correlation with confidence intervals; plus the 
# trait standard deviations


corp <-function(X, phy, SeM, U = NULL){

  #####Main program#####
  if (!inherits(phy, "phylo"))
    stop("Object \"phy\" is not of class \"phylo\".")
  if (is.null(phy$edge.length))
    stop("The tree has no branch lengths.")
  if (is.null(phy$tip.label))
    stop("The tree has no tip labels.")
  phy <- reorder(phy, "postorder")
  n <- length(phy$tip.label)



  ##Input X
  if (dim(X)[1] != n)
    stop("Number of rows of the data matrix does not match the length of the tree.")
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




  ##Input U
  if (length(U) > 0) {
    if (length(U) != p)
      stop("Number of elements of list U does not match the number of columns in X.")

    for (i in 1:p) {
      if (!is.null(U[[i]])){
        if (dim(U[[i]])[1] != n)
          stop("Number of rows of an element of U does not match the tree.")
        if (is.null(rownames(U[[i]]))) {
          warning("No tip labels on U; order assumed to be the same as in the tree.\n")
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
        U[[i]] <- matrix(0, nrow=n, ncol=1)
        rownames(U[[i]]) <- phy$tip.label
      }
    }
  }

  # SeM <-NULL
  # Input SeM
  if (!is.null(SeM)) {
    if (dim(SeM)[1] != n)
      stop("Number of rows of the SeM matrix does not match the length of the tree.")
    if (is.null(rownames(SeM))) {
      warning("No tip labels on SeM; order assumed to be the same as in the tree.\n")
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




  ##Standardize all variables
  Xs <- X
  for (i in 1:p) Xs[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])

  if (!is.null(SeM)) {
    SeMs <- SeM
    for (i in 1:p) SeMs[, i] <- SeM[, i]/sd(X[, i])
  }

  if (length(U) > 0) {
    Us <- U
    for (i in 1:p) for (j in 1:ncol(U[[i]])) {
      if (sd(U[[i]][, j]) > 0) {
        Us[[i]][, j] <- (U[[i]][, j] - mean(U[[i]][, j]))/sd(U[[i]][, j])
      } else {
        Us[[i]][, j] <- U[[i]][, j] - mean(U[[i]][, j])
      }
    }
  }



  ##Set up matrices
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



  # Compute initial estimates assuming no phylogeny if not provided
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
  par <- c(L.elements, array(0.5, dim = c(1, p)))
  names(par) <-c("sig1","r","sig2","d1","d2")
  tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy






#####All of the above is the same as the original 'corphylo'

  corphylo.LL <- function(par, XX, UU, MM, tau, Vphy, REML, constrain.d, verbose) {

    n <- nrow(X)
    p <- ncol(X)

    ##This is where my modifications begin
    RR <-matrix(0, nrow=p, ncol=p)
    RR[1,2] <-par["r"]
    RR[2,1] <-par["r"]
    diag(RR) <-1
    D <-diag(par[c("sig1","sig2")])
    R <-D %*% RR %*% D

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

    # OU transform
    C <- matrix(0, nrow = p * n, ncol = p * n)
    for (i in 1:p) for (j in 1:p) {
      Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
      C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
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
      # REML likelihood function
      LL <- 0.5 * (logdetV + determinant(t(UU) %*% iV %*% UU)$modulus[1] + 
                       t(H) %*% iV %*% H)
    } else {
      # ML likelihood function
      LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
    }

    if (verbose == T)
      show(c(as.numeric(LL), par))
    return(as.numeric(LL))
  }
  # End corphylo.LL

  # For SANN
  reltol = 10^-6
  maxit.NM = 1000
  maxit.SA = 1000
  temp.SA = 1
  tmax.SA = 1
  constrain.d = F
  # opt <- optim(fn = corphylo.LL, par = par, XX = XX, UU = UU, MM = MM, tau = tau, 
  #              Vphy = Vphy, verbose = F, constrain.d = constrain.d, REML = T, 
  #              method = "SANN", control = list(maxit = maxit.SA,temp = temp.SA, 
  #                                              tmax = tmax.SA, reltol = reltol), 
  #              hessian=T)
  # par <- opt$par
  opt <- optim(fn = corphylo.LL, par = par, XX = XX, UU = UU, MM = MM, tau = tau, 
               Vphy = Vphy, verbose = F, constrain.d = constrain.d, REML = T, 
               method = "Nelder-Mead", 
               control = list(maxit = maxit.NM, reltol = reltol), hessian=T)

###Extract parameters and calculate standard errors
  fisher_info<-solve(opt$hessian)
  prop_sigma<-sqrt(diag(fisher_info))
  upper<-opt$par + 1.96 * prop_sigma
  lower<-opt$par - 1.96 * prop_sigma
  df <-data.frame(value=opt$par, upper=upper, lower=lower)
  # df[4:5, ] <-1/(1 + exp(-df[4:5, ]))
  df
}
