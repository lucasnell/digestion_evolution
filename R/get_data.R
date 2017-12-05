# 
# Functions for retrieving data of mean, standard errors, or phylogenies among species
# 


get_df <- function(.df, .pos = NA, .stat = c('mean', 'se')) {
    
    .df <- match.arg(.df, c('spp', 'clear', 'absorp', 'pos'))
    if (!is.na(.pos)) .pos  <- match.arg(.pos, c('prox', 'med', 'dist'))
    .stat <- match.arg(.stat)
    
    if (.df == 'pos' & is.na(.pos)) stop('provide a position when .df == "pos"')
    if (.df != 'pos' & !is.na(.pos)) warning('.df != "pos" so position is ignored')
    if (!.df %in% c('clear', 'absorp') & .stat == 'se') {
        message("You shouldn't need SE for this analysis set. Returning anyway...")
    }
    
    fn <- sprintf('output/tidy_%s.csv', .df)
    
    out_df <- suppressMessages(readr::read_csv(fn)) %>% 
        dplyr::filter(stat == .stat) %>% 
        dplyr::select(-stat)
    
    if (!is.na(.pos) & .df == 'pos') {
        out_df <- dplyr::filter(out_df, pos == .pos) %>% 
            dplyr::select(-pos)
    }
    
    # Making sure they're all sorted the same
    out_df <- dplyr::arrange(out_df, clade, species)
    
    # Converting clade to factor
    out_df$clade <-  factor(out_df$clade, levels = c('Rodent', 'Bat'))
    
    # Converting diet to factor, if present
    if (.df %in% c('spp', 'pos')) {
        out_df$diet <-  factor(out_df$diet, levels = c('Herbivorous', 'Omnivorous', 
                                                       'Protein'))
    } else if (.df == 'clear') {
        out_df$diet <-  factor(out_df$diet, levels = c('Carb', 'Protein'))
    }
    
    # phylolm requires that the rownames match the species names
    # (tibbles cannot have row names)
    out_df <- as.data.frame(out_df)
    rownames(out_df) <- out_df$species
    
    return(out_df)
}


get_tr <- function(.df) {
    
    # Just choosing parameters for get_df. These won't affect anything bc this is just
    # to retrieve species names.
    .pos <- ifelse(.df == 'pos', 'med', NA)

    focal_df <- get_df(.df, .pos = .pos)

    out_tr <- ape::read.tree('data/tree.nwk')
    out_tr$tip.label <- gsub('_', ' ', out_tr$tip.label)
    tips_to_rm <- out_tr$tip.label[!out_tr$tip.label %in% focal_df$species]
    out_tr <- ape::drop.tip(out_tr, tip = tips_to_rm)
    
    return(out_tr)
}



cp_mat <- function(.df, .pars) {
    out <- as.matrix(.df[,.pars])
    colnames(out) <- NULL
    rownames(out) <- .df$species
    return(out)
}

filter_tr <- function(tr, keep_spp) {
    cleared_tr <- ape::drop.tip(
        tr, 
        tip = tr$tip.label[!tr$tip.label %in% keep_spp]
    )
    return(cleared_tr)
    
}
