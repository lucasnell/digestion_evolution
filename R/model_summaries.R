# 
# Functions for summarizing phylolm and/or corphylo models
# 



#' Coefficient P-value(s) for a \code{phylolm} or \code{corphylo} object.
#'
#' @param .model A \code{phylolm} or \code{corphylo} object that has bootstrap 
#'     replicates within (i.e., was run with \code{boot > 0}).
#' @param .parameters A character vector of the parameter name(s) of interest. 
#'     Defaults to \code{'cladeBat'}.
#'
#' @return A numeric vector of p-values for each parameter not equalling zero.
#' 
#' @export
#'
pval <- function(.model, .parameters = 'cladeBat') {
    if (length(.model$bootstrap) == 0) {
        stop("The input model must include bootstrap replicates.")
    }
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

#' Coefficient confidence interval(s) for a \code{phylolm} or \code{corphylo} object.
#'
#' @inheritParams pval
#'
#' @return A two-column matrix of confidence intervals for each parameter, where the
#'     first colum is the lower limit and the second is the upper.
#' 
#' @export
#'
ci <- function(.model, .parameters = 'cladeBat') {
    if (is(.model, 'phylolm')) {
        out <- matrix(as.numeric(.model$bootconfint95[,.parameters]), 
                      ncol = 2, nrow = length(.parameters), byrow = TRUE)
    } else if (is(.model, 'corphylo')) {
        out <- t(apply(.model$bootstrap, 2, quantile, probs = c(0.025, 0.975)))
        rownames(out) <- NULL
    } else {
        stop("only for phylolm or corphylo objects")
    }
    colnames(out) <- c('lower', 'upper')
    return(out)
}



#' Extract phylogenetic signal string from a \code{phylolm} object.
#'
#' @param .model \code{phylolm} object.
#'
#' @return A string specifying the name of the coefficient for the phylogenetic signal.
#' @export
#'
phy_signal_str <- function(.model) {
    switch(.model$model,
           OUfixedRoot = "alpha",
           OUrandomRoot = "alpha",
           EB = "rate",
           .model$model)
}


#' Summarize a \code{phylolm} or \code{corphylo} object.
#' 
#'
#' @param .model A \code{phylolm} or \code{corphylo} object that has bootstrap 
#'     replicates within (i.e., was run with \code{boot > 0}).
#' @param .pos 
#' @param .corr_pars 
#'
#' @return
#' @export
#'
#' @examples
summ_df <- function(.model, .pos = NA, .corr_pars = NA) {
    
    if (is(.model, 'phylolm')) {
        error_model <- .model$model
        .parameters <- names(coef(.model))[names(coef(.model)) != '(Intercept)']
        estimates <- as.numeric(c(coef(.model)[.parameters], .model$optpar))
        .parameters <- c(.parameters, phy_signal_str(.model))
        Y <- paste(.model$formula)[2]
    } else if (is(.model, 'corphylo')) {
        if (length(.corr_pars) != 2) stop("Correlation parameters must have length == 2")
        error_model <- "OU"
        Y <- c(.corr_pars[1], .corr_pars)
        .parameters <- c(.corr_pars[2], rep('d', 2))
        estimates <- c(.model$cor.matrix[lower.tri(.model$cor.matrix)], .model$d)
    } else {
        stop("only for phylolm or corphylo objects")
    }
    
    .df <- as_data_frame(ci(.model, .parameters))
    .df <- .df %>% 
        mutate(Y = Y, 
               X = .parameters,
               pos = .pos,
               phy_model = error_model,
               value = estimates,
               P = pval(.model, .parameters)) %>% 
        select(Y, X, pos, phy_model, value, lower, upper, P)
    
    return(.df)
}






#' Jackknifing on a \code{phylolm} object to find influential points.
#' 
#' \emph{Note}: this function only works for my specific analyses.
#' The portion of this function dealing with factors is particularly non-generalizable,
#' plus it assumes an intercept is estimated.
#'
#' @param plm A \code{phylolm} object.
#' @param phy The \code{phylo} object, the same one that was used to fit the original
#'     \code{phylolm} model.
#'
#' @return An `n` by `p` data frame, where `n` is the sample size and `p` is the 
#'     number of coefficients estimated by the \code{phylolm} model. 
#'     A cell in row `i` and column `j` is the jackknife influence value on 
#'     coefficient `j` for removing row `i` from the analysis.
#' 
#' @export
#'
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
                          .to_drop <- phy$tip.label
                          .to_drop <- .to_drop[!.to_drop %in% rownames(.df)]
                          .tr <- ape::drop.tip(phy, .to_drop)
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







#' Create data frame of predictions and 95% CI from a \code{phylolm} object.
#'
#' @param .model \code{phylolm} object that has bootstrap replicates within 
#'     (i.e., was run with \code{boot > 0}).
#'
#' @return A dataframe with columns for clade, intestinal segment position, measurement 
#'     name, predicted estimate for measurement, lower limit of 95% CI for prediction, 
#'     and upper limit of 95% CI, respectively. 
#'     
#'     For models with log(mass) as a covariate, there should be 200 rows, 
#'     100 for the bat clade and 100 for the rodent clade.
#'     Within each clade in the output, there are 100 rows of log(mass) ranging between
#'     the minimum and maximum values of the input data for that clade.
#'     
#'     For models without log(mass) as a covariate, there should be two rows, one
#'     for bats and one for rodents.
#'     In these models, both rows have log(mass) set to zero.
#' 
#' @export
#'
predict_ci <- function(.model){
    
    stopifnot(is(.model, 'phylolm'))

    y_measure <- {paste(.model$formula) %>% purrr::discard(~ grepl('~', .x))}[1]
    
    # Name of position, if applicable
    pos_name <- paste(.model$call) %>% 
        purrr::keep(~ grepl('_df', .x)) %>% 
        gsub(pattern = '_df', replacement = '')
    if (pos_name %in% c('spp', 'absorp')) pos_name <- NA
    
    # Creating prediction data frames
    # If log_mass is a covariate in the model, I'm using, for each clade, a sequence 
    # ranging between the minimum and maximum values of the input data for that clade.
    # The new_data matrix is used for matrix multiplication for estimates later on.
    if ('log_mass' %in% colnames(.model$X)) {
        n_mass <- 100
        bat_log_mass <- as.numeric(.model$X[.model$X[,'cladeBat'] == 1,'log_mass'])
        rod_log_mass <- as.numeric(.model$X[.model$X[,'cladeBat'] == 0,'log_mass'])
        log_mass_ <- c(seq(min(bat_log_mass), max(bat_log_mass), length.out = n_mass),
                      seq(min(rod_log_mass), max(rod_log_mass), length.out = n_mass))
        new_data <- rbind(1.0, rep(c(1.0, 0.0), each = n_mass), log_mass_)
        rownames(new_data) <- NULL
        new_data_df <- data.frame(clade = rep(c(1.0, 0.0), each = n_mass), 
                                  log_mass = log_mass_)
    } else {
        log_mass_ <- 0
        n_mass <- 1
        new_data <- rbind(c(1.0, 1.0), c(1.0, 0.0))
        new_data_df <- data.frame(clade = c(1.0, 0.0))
    }
    
    # Column names coinciding with phylogenetic parameters from the model:
    phylo_cols <- which(colnames(.model$bootstrap) %in% 
                            c('sigma2', phy_signal_str(.model)))
    
    # Calculating confidence intervals for predicted y-value for bats, then rodents.
    # Bats are the first `n_mass` rows, rodents the second.
    # Columns are lower and upper limits, respectively.
    ci_matrix <- .model$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(x[-phylo_cols] %*% new_data)
              }) %>% 
        apply(1,
              function(x) as.numeric(quantile(x, probs = c(0.025, 0.975)))) %>% 
        t
    
    # Adding predicted outcomes, and converting to a data frame.
    out_df <- cbind(predict(.model, new_data_df), ci_matrix) %>% 
        as_data_frame %>% 
        rename_(.dots = setNames(colnames(.), c('estimate', 'low', 'high'))) %>% 
        mutate(measure = y_measure,
               clade = factor(rep(c(1.0, 0.0), each = n_mass), levels = c(0,1), 
                              labels = c('Rodent', 'Bat')),
               log_mass = log_mass_) %>% 
        select(clade, measure, everything()) %>% 
        mutate(pos = pos_name) %>% 
        select(clade, pos, measure, everything())
    
    return(out_df)
}



