# 
# Functions for retrieving data of mean, standard errors, or phylogenies among species
# 


get_df <- function(.df, .pos = NA) {
    
    .df <- match.arg(.df, c('spp', 'clear', 'absorp', 'pos'))
    if (!is.na(.pos)) .pos  <- match.arg(.pos, c('prox', 'med', 'dist'))
    
    if (.df == 'pos' & is.na(.pos)) stop('provide a position when .df == "pos"')
    if (.df != 'pos' & !is.na(.pos)) warning('.df != "pos" so position is ignored')
    
    fn <- paste0('output/tidy_', .df, '.csv')
    
    out_df <- suppressMessages(readr::read_csv(fn))
    
    
    if (!is.na(.pos) & .df == 'pos') {
        out_df <- dplyr::filter(out_df, pos == .pos)
        out_df <- dplyr::select(out_df, -pos)
    }
    
    # Making sure they're all sorted the same
    out_df <- dplyr::arrange(out_df, taxon, species)
    
    # Converting taxon to factor
    out_df$taxon <-  factor(out_df$taxon, levels = c('Rodent', 'Bat'))
    
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
    
    # Just choosing parameters for get_df. These won't affect anything.
    .pos <- ifelse(.df == 'pos', 'med', NA)

    focal_df <- get_df(.df, .pos = .pos)

    out_tr <- ape::read.tree('data/tree.nwk')
    out_tr$tip.label <- gsub('_', ' ', out_tr$tip.label)
    tips_to_rm <- out_tr$tip.label[!out_tr$tip.label %in% focal_df$species]
    out_tr <- ape::drop.tip(out_tr, tip = tips_to_rm)
    
    return(out_tr)
}

