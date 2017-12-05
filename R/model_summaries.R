# 
# Functions for summarizing phylolm models
# 

pval <- function(.model, .parameters = 'cladeBat') {
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

ci <- function(.model, .parameters = 'cladeBat') {
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





# Function to do jackknifing on phylolm object to find influential points (i.e., species)
jack_phylolm <- function(plm, phy) {
    
    return_fxn <- function(.plm) return(t(coef(.plm)))
    
    X_names <- gsub("\\s+", "", unlist(strsplit(paste(plm$formula)[3], '\\+|\\*')))
    
    df_ <- cbind(plm$X[,-1], as.data.frame(plm$y))
    colnames(df_)[length(colnames(df_))] <- paste(plm$formula)[2]
    # Dealing with factors
    for (xn in X_names) {
        ii <- which(grepl(xn, colnames(df_)))
        if (colnames(df_)[ii][1] == 'cladeBat') {
            # This one must be clade
            new_col <- factor(ifelse(df_$cladeBat == 1, 'Bat', 'Rodent'),
                              levels = c('Rodent', 'Bat'))
            df_[,xn] <- new_col
            df_ <- df_[,-ii]
        } else if (length(ii) == 2) {
            # This one must be diet
            new_col <- factor(ifelse(df_$dietOmnivorous == 1, 'Omnivorous', 
                                     ifelse(df_$dietProtein == 1, 
                                            'Protein', 'Herbivorous')),
                              levels = c('Herbivorous', 'Omnivorous', 'Protein'))
            df_[,xn] <- new_col
            df_ <- df_[,-ii]
        }
    }; rm(xn, ii, new_col)
    
    comp_ <- return_fxn(plm)
    
    jack_df <- lapply(1:nrow(df_), 
                      function(i) {
                          .df <- df_[-i,]
                          .tips_to_drop <- phy$tip.label[!phy$tip.label %in% rownames(.df)]
                          .tr <- ape::drop.tip(phy, .tips_to_drop)
                          suppressWarnings(
                              .plm <- update(plm, data = .df, phy = .tr, boot = 0)
                          )
                          return(return_fxn(.plm) - comp_)
                      }) %>% 
        do.call(what = rbind) %>% 
        as_tibble %>% 
        rename(intercept = `(Intercept)`) %>% 
        mutate(species = rownames(df_)) %>% 
        gather('estimate', 'influence', -species)
    
    return(jack_df)
}


