# 
# Functions for summarizing phylolm models
# 

pval <- function(model, parameter = 'taxonBat') {
    if (is(model, 'phylolm')) {
        2 * min(c(mean(model$bootstrap[,parameter] > 0), 
                  mean(model$bootstrap[,parameter] < 0)))
    } else if (is(model, 'corphylo')) {
        2 * min(c(mean(model$bootstrap[,parameter] > 0), 
                  mean(model$bootstrap[,parameter] < 0)))
    } else {
        stop("only for phylolm or corphylo objects")
    }
}

ci <- function(model, parameter = 'taxonBat') model$bootconfint95[,parameter]

ci_df <- function(.model, .pos = NA) {
    
    .parameters <- names(coef(.model))[names(coef(.model)) != '(Intercept)']
    
    .df <- as_data_frame(t(ci(.model, .parameters)))
    colnames(.df) <- c('lower', 'upper')
    .df <- .df %>% mutate(X = .parameters,
                          Y = paste(.model$formula)[2], 
                          pos = .pos,
                          value = coef(.model)[.parameters]) %>% 
        select(Y, X, pos, value, lower, upper)
    
    return(.df)
}

# Return AICc for a phylolm model
aicc <- function(m) {
    n = m$n
    k = m$p
    return(m$aic + { 2*k*(k+1) } / {n - k - 1})
}
