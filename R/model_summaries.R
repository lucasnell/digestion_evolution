# 
# Functions for summarizing phylolm models
# 

pval <- function(.model, .parameters = 'taxonBat') {
    p_fun <- function(x) 2 * min(c(mean(x > 0), mean(x < 0)))
    if (is(.model, 'phylolm')) {
        out <- lapply(.parameters, function(p) p_fun(.model$bootstrap[, p])) %>% 
            c(recursive = TRUE)
    } else if (is(.model, 'corphylo')) {
        out <- as.numeric(apply(.model$bootstrap, 2, p_fun))
    } else {
        stop("only for phylolm or corphylo objects")
    }
    return(out)
}

ci <- function(.model, .parameters = 'taxonBat') {
    if (is(.model, 'phylolm')) {
        out <- matrix(as.numeric(.model$bootconfint95[,.parameters]), 
                      ncol = 2, nrow = length(.parameters), byrow = TRUE)
        colnames(out) <- c('lower', 'upper')
    } else if (is(.model, 'corphylo')) {
        out <- t(apply(.model$bootstrap, 2, quantile, probs = c(0.025, 0.975)))
        rownames(out) <- NULL
        colnames(out) <- c('lower', 'upper')
    } else {
        stop("only for phylolm or corphylo objects")
    }
    return(out)
}

# Row(s) of a data frame for a single model
summ_df <- function(.model, .pos = NA, .corr_pars = NA) {
    
    
    if (is(.model, 'phylolm')) {
        .parameters <- names(coef(.model))[names(coef(.model)) != '(Intercept)']
        estimates <- as.numeric(c(coef(.model)[.parameters], .model$optpar))
        .parameters <- c(.parameters, 'lambda')
        Y <- paste(.model$formula)[2]
    } else if (is(.model, 'corphylo')) {
        if (length(.corr_pars) != 2) stop("Correlation parameters must have length == 2")
        Y <- c(.corr_pars[1], .corr_pars)
        .parameters <- c(.corr_pars[2], rep('d', 2))
        estimates <- c(.model$cor.matrix[lower.tri(.model$cor.matrix)], .model$d)
    } else {
        stop("only for phylolm or corphylo objects")
    }
    
    
    .df <- as_data_frame(ci(.model, .parameters))
    .df <- .df %>% 
        mutate(X = .parameters,
               Y = Y, 
               pos = .pos,
               value = estimates,
               P = pval(.model, .parameters)) %>% 
        select(Y, X, pos, value, lower, upper, P)
    
    return(.df)
}
