# 
# Functions for summarizing phylolm and/or corphylo models
# 



#' Coefficient P-value(s) for a \code{phylolm} or \code{corphylo} object.
#'
#' @param mod A \code{phylolm} or \code{corphylo} object that has bootstrap 
#'     replicates within (i.e., was run with \code{boot > 0}).
#' @param params A character vector of the parameter name(s) of interest. 
#'     Defaults to \code{'cladeBat'}.
#'
#' @return A numeric vector of p-values for each parameter not equalling zero.
#' 
#' @export
#'
pval <- function(mod, params) UseMethod("pval")

pval.phylolm <- function(mod, params = 'cladeBat') {
    if (length(mod$bootstrap) == 0) {
        stop("The input model must include bootstrap replicates.")
    }
    p_fun <- function(x) 2 * min(c(mean(x > 0), mean(x < 0)))
    out <- lapply(params, function(p) p_fun(mod$bootstrap[, p]))
    out <- c(out, recursive = TRUE)
    return(out)
}

pval.cor_phylo <- function(mod, params) {
    lower.to.upper.tri.inds <- function (n) {
        n1 <- as.integer(n - 1)
        if (n1 < 1) {
            stop("'n' must be >= 2")
        } else if (n1 == 1) {
            1L
        } else {
            fxn <- function(k) cumsum(c(0L, (n - 2L):(n - k)))
            rep(seq_len(n1), seq_len(n1)) + c(0L, unlist(lapply(2:n1, fxn)))
        }
    }
    
    if (length(mod$bootstrap) == 0) {
        stop("The input model must include bootstrap replicates.")
    }
    p_fun <- function(x) 2 * min(c(mean(x > 0), mean(x < 0)))
    
    # Indices for failed convergences:
    if (phyr:::call_arg(mod$call,"method") %in% c("nelder-mead-r", "sann")) {
        f <- mod$bootstrap$inds[mod$bootstrap$codes != 0]
    } else {
        f <- mod$bootstrap$inds[mod$bootstrap$codes < 0]
    }
    if (length(f) == 0) f <- ncol(mod$bootstrap$d) + 1
    
    output <- rep(list(NA), length(params))
    for (i in 1:length(params)) {
        if (params[i] == "corrs") {
            
            lowers <- lapply({1:(dim(mod$bootstrap$corrs)[3])}[-f],
                             function(i) {
                                 logs <- lower.tri(mod$bootstrap$corrs[,,i])
                                 mod$bootstrap$corrs[,,i][logs]
                             })
            lowers <- do.call(rbind, lowers)
            
            out <- matrix(-1, nrow(mod$bootstrap$corrs), ncol(mod$bootstrap$corrs))
            out[lower.tri(out)] <- apply(lowers, 2, p_fun)
            out[upper.tri(out)] <- out[lower.tri(out)][lower.to.upper.tri.inds(ncol(out))]
            
            rownames(out) <- rownames(mod$corrs)
            colnames(out) <- colnames(mod$corrs)
        } else if (params[i] == "d") {
            out <- NA
        } else if (params[i] == "B0") {
            
            out <- cbind(apply(mod$bootstrap$B0[,-f,drop=FALSE], 1, p_fun))
            rownames(out) <- rownames(mod$B)
            colnames(out) <- "P-value"
            
        } else if (params[i] == "B_cov") {
            
            lowers <- lapply({1:(dim(mod$bootstrap$B_cov)[3])}[-f],
                             function(i) {
                                 logs <- lower.tri(mod$bootstrap$B_cov[,,i])
                                 mod$bootstrap$B_cov[,,i][logs]
                             })
            lowers <- do.call(rbind, lowers)
            
            out <- matrix(-1, nrow(mod$bootstrap$B_cov), ncol(mod$bootstrap$B_cov))
            out[lower.tri(out)] <- apply(lowers, 2, p_fun)
            out[upper.tri(out)] <- out[lower.tri(out)][lower.to.upper.tri.inds(ncol(out))]
            
            rownames(out) <- rownames(mod$B_cov)
            colnames(out) <- colnames(mod$B_cov)
        } else out <- NA
        output[[i]] <- out
    }
    if (length(output) == 1) output <- output[[1]]


    return(output)
}



#' Coefficient confidence interval(s) for a \code{phylolm} or \code{corphylo} object.
#'
#' @inheritParams pval
#'
#' @return A matrix of confidence intervals for each parameter.
#'     In most instances, it returns two-column matrices, where the first column 
#'     is the lower limit and the second is the upper.
#'     For `corrs` and `B_cov` fields in a `cor_phylo` object, it returns
#'     a square matrix, where values above the diagonal are upper limits and
#'     values below the diagonal are lower limits.
#' 
#' @export
#'
ci <- function(mod, params) UseMethod("ci")

ci.phylolm <- function(mod, params = "cladeBat") {
    
    out <- matrix(as.numeric(mod$bootconfint95[,params]), 
                  ncol = 2, nrow = length(params), byrow = TRUE)
    
    colnames(out) <- c("lower", "upper")
    return(out)
}

ci.cor_phylo <- function(mod, params) {
    if (any(!params %in% c("corrs", "d", "B0", "B_cov"))) {
        stop("params can only be corrs, d, B0, or B_cov for cor_phylo.")
    }
    out <- boot_ci(mod)
    output <- lapply(params, function(pp) out[[pp]])
    if (length(output) == 1) output <- output[[1]]
    return(output)
}


#' Extract phylogenetic signal string from a \code{phylolm} object.
#'
#' @param mod \code{phylolm} object.
#'
#' @return A string specifying the name of the coefficient for the phylogenetic signal.
#' @export
#'
phy_signal_str <- function(mod) {
    switch(mod$model,
           OUfixedRoot = "alpha",
           OUrandomRoot = "alpha",
           EB = "rate",
           mod$model)
}


#' Summarize a \code{phylolm} or \code{corphylo} object.
#' 
#'
#' @param mod A \code{phylolm} or \code{corphylo} object that has bootstrap 
#'     replicates within (i.e., was run with \code{boot > 0}).
#' @param .pos 
#'
#' @return
#' @export
#'
#' @examples
summ_df <- function(mod, .pos = NA) {
    
    if (inherits(mod, 'phylolm')) {
        error_model <- mod$model
        params <- names(coef(mod))[names(coef(mod)) != '(Intercept)']
        estimates <- as.numeric(c(coef(mod)[params], mod$optpar))
        params <- c(params, phy_signal_str(mod))
        Y <- paste(mod$formula)[2]
        .df <- as_data_frame(ci(mod, params))
        P <- pval(mod, params)
    } else if (inherits(mod, 'cor_phylo')) {
        .corr_pars <- rownames(mod$d)
        if (length(.corr_pars) > 2) stop("\nNot designed for >2 parameters.")
        error_model <- "OU"
        Y <- c(.corr_pars[1], .corr_pars)
        params <- c(.corr_pars[2], rep('d', 2))
        estimates <- c(mod$corrs[lower.tri(mod$corrs)], mod$d)
        
        corrs <- ci(mod, "corrs")
        corrs <- cbind(lower = corrs[lower.tri(corrs)],
                       upper = corrs[upper.tri(corrs)])
        ds <- ci(mod, "d")
        
        .df <- as_data_frame(rbind(corrs, ds))
        
        P <- pval(mod, c("corrs", rep("d", length(.corr_pars))))
        P <- c(lapply(P,
                      function(x) {
                          if (inherits(x, "matrix")) x[lower.tri(x)] else x
                      }),
               recursive = TRUE)
        
    } else {
        stop("only for phylolm or cor_phylo objects")
    }
    
    .df <- .df %>%
        mutate(Y = Y, 
               X = params,
               pos = .pos,
               phy_model = error_model,
               value = estimates,
               P = P) %>% 
        select(Y, X, pos, phy_model, value, lower, upper, P)
    
    return(.df)
}






#' Jackknifing on a \code{phylolm} object to find influential points.
#' 
#' \emph{Note}: this function only works for my specific analyses.
#' The portion of this function dealing with factors is particularly non-generalizable,
#' plus it assumes an intercept is estimated.
#'
#' @param mod A \code{phylolm} object.
#' @param phy A \code{phylo} object, the same one that was used to fit the original
#'     \code{phylolm} model.
#'
#' @return An `n` by `p` data frame, where `n` is the sample size and `p` is the 
#'     number of coefficients estimated by the \code{phylolm} model. 
#'     A cell in row `i` and column `j` is the jackknife influence value on 
#'     coefficient `j` for removing row `i` from the analysis.
#' 
#' @export
#'
jack <- function(mod, ...) UseMethod("jack")

jack.phylolm <- function(mod, phy) {
    
    return_fxn <- function(.mod) return(t(coef(.mod)))
    
    X_names <- gsub("\\s+", "", unlist(strsplit(paste(mod$formula)[3], '\\+|\\*')))
    
    df_ <- cbind(mod$X[,-1], as.data.frame(mod$y))
    colnames(df_)[length(colnames(df_))] <- paste(mod$formula)[2]
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
    
    comp_ <- return_fxn(mod)
    
    jack_df <- lapply(1:nrow(df_), 
                      function(i) {
                          .df <- df_[-i,]
                          .to_drop <- phy$tip.label
                          .to_drop <- .to_drop[!.to_drop %in% rownames(.df)]
                          .tr <- ape::drop.tip(phy, .to_drop)
                          suppressWarnings(
                              .mod <- update(mod, data = .df, phy = .tr, boot = 0)
                          )
                          return(return_fxn(.mod) - comp_)
                      }) %>% 
        do.call(what = rbind) %>% 
        as_tibble %>% 
        # Standardizing to SD of 1 and mean = 0
        mutate_if(is.numeric, function(x) (x - mean(x)) / sd(x)) %>% 
        rename(intercept = `(Intercept)`) %>% 
        mutate(species = rownames(df_)) %>% 
        gather('estimate', 'influence', -species)
    
    return(jack_df)
}


#' Jackknifing on a \code{cor_phylo} object to find influential points.
#' 
#' All arguments to this function should be the same as used for the original
#' call to \code{cor_phylo} or \code{corphylo_cpp}
#' 
#' \emph{Note}: this function only works for my specific analyses.
#' 
#'
#' @param mod A \code{cor_phylo} object.
#' @param mean_df A data frame of mean values by species.
#' @param se_df A data frame of standard errors by species.
#' @param phy A \code{phylo} object of all species.
#' @param par_names Names of columns for the parameters of interest.
#' 
#'
#' @return An `n` by 1 data frame, where `n` is the sample size. 
#'     A cell in row `i` is the jackknife influence value on 
#'     the correlation for removing row `i` from the analysis.
#' 
#' @export
#'
jack.cor_phylo <- function(mod) {
    
    comp_ <- mod$corrs[1,2]
    call_ <- mod$call
    call_$boot <- 0
    call_$data <- quote(clear_df_)
    call_$phy <- quote(clear_tr_)
    
    n <- nrow(eval(mod$call$data))
    
    out <- rep(list(NA), n)
    
    for (i in 1:n) {
        clear_df_ <- eval(mod$call$data)[-i,]
        to_drop_ <- eval(mod$call$phy)$tip.label
        to_drop_ <- to_drop_[!to_drop_ %in% clear_df_$species]
        clear_tr_ <- ape::drop.tip(eval(mod$call$phy), to_drop_)
        cp_ <- eval(call_)
        out[[i]] <- c(influence = cp_$corrs[1,2] - comp_)
    }
    
    out <- out %>%
        do.call(what = rbind) %>%
        as_tibble %>%
        # Standardizing to SD of 1 and mean = 0
        mutate_if(is.numeric, function(x) (x - mean(x)) / sd(x)) %>%
        mutate(species = eval(mod$call$data)$species)
    
    return(out)
}







#' Create data frame of predictions and 95% CI from a \code{phylolm} object.
#'
#' @param mod \code{phylolm} object that has bootstrap replicates within 
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
predict_ci <- function(mod){
    
    stopifnot(is(mod, 'phylolm'))

    y_measure <- {paste(mod$formula) %>% purrr::discard(~ grepl('~', .x))}[1]
    
    # Name of position, if applicable
    pos_name <- paste(mod$call) %>% 
        purrr::keep(~ grepl('_df', .x)) %>% 
        gsub(pattern = '_df', replacement = '')
    if (pos_name %in% c('spp', 'absorp')) pos_name <- NA
    
    # Creating prediction data frames
    # If log_mass is a covariate in the model, I'm using, for each clade, a sequence 
    # ranging between the minimum and maximum values of the input data for that clade.
    # The new_data matrix is used for matrix multiplication for estimates later on.
    if ('log_mass' %in% colnames(mod$X)) {
        n_mass <- 100
        bat_log_mass <- as.numeric(mod$X[mod$X[,'cladeBat'] == 1,'log_mass'])
        rod_log_mass <- as.numeric(mod$X[mod$X[,'cladeBat'] == 0,'log_mass'])
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
    phylo_cols <- which(colnames(mod$bootstrap) %in% 
                            c('sigma2', phy_signal_str(mod)))
    
    # Calculating confidence intervals for predicted y-value for bats, then rodents.
    # Bats are the first `n_mass` rows, rodents the second.
    # Columns are lower and upper limits, respectively.
    ci_matrix <- mod$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(x[-phylo_cols] %*% new_data)
              }) %>% 
        apply(1,
              function(x) as.numeric(quantile(x, probs = c(0.025, 0.975)))) %>% 
        t
    
    # Adding predicted outcomes, and converting to a data frame.
    out_df <- cbind(predict(mod, new_data_df), ci_matrix) %>% 
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




#' Add title to the top left corner of a plot.
#'
#' @param .p \code{ggplot} object to add the title to.
#' @param .title String to add as the plot title.
#' @param .mult List specifying multipliers for offset for x- and y-axes. 
#'     Offsets are the amount of space from the upper-left corner of the plot 
#'     the upper-left corner of the title text. Defaults to 1 for both.
#' @param .data Data frame used for the call to `geom_text`. This is often useful
#'     for faceted plots to make sure the title ends up on only one facet.
#'     Defaults to \code{NULL}.
#' @param font_size Numeric font size for the title (in pts). Defaults to \code{14}.
#'
#' @return
#' @export
#'
#' @examples
add_title <- function(.p, .title, .mult = list(x = 1, y = 1), .data = NULL, 
                      font_size = 14) {
    
    if (missing(.title)) stop("Why is .title missing?")
    if (is.null(.title)) return(.p)
    
    x_range <- ggplot_build(.p)$layout$panel_ranges[[1]]$x.range
    y_range <- ggplot_build(.p)$layout$panel_ranges[[1]]$y.range
    if (!is.null(.p$coordinates$trans$x)) {
        x_range <- .p$coordinates$trans$x$inverse(x_range)
    }
    if (!is.null(.p$coordinates$trans$y)) {
        y_range <- .p$coordinates$trans$y$inverse(y_range)
    }
    min_x <- min(x_range) + .mult$x * 0.02 * diff(x_range)
    max_y <- max(y_range) - .mult$y * 0.02 * diff(y_range)
    .p <- .p +
        geom_text(data = .data, 
                  label = .title,
                  x = min_x, y = max_y,
                  color = 'black',
                  hjust = 0, vjust = 1, 
                  size = font_size * (25.4/72), # <-- 25.4/72 is to convert from mm to pt
                  fontface = 'plain')
    return(.p)
}