# 
# Functions for summarizing phylolm models
# 

pval <- function(model, parameter = 'taxonBat') {
    2 * min(c(mean(model$bootstrap[,parameter] > 0), 
              mean(model$bootstrap[,parameter] < 0)))
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