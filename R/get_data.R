# 
# Function for retrieving data of mean or standard errors among species
# 


get_df <- function(.df, .stat = c('mean', 'se'), 
                   .pos = c('none', 'prox', 'med', 'dist')) {
    
    .df <- match.arg(.df, c('spp', 'diet', 'clear', 'absorp', 'pos'))
    .stat  <- match.arg(.stat)
    .pos  <- match.arg(.pos)
    
    if (.df == 'pos' & .pos == 'none') stop('provide non-"none" .pos when .df == "pos"')
    if (.df != 'pos' & .pos != 'none') warning('.df != "pos" so .pos argument ignored')
    
    fn <- paste0('output/', ifelse(.df == 'diet', 'spp', .df), '_df.csv')
    
    out_df <- suppressMessages(readr::read_csv(fn))
    out_df <- as.data.frame(out_df)
    
    
    if (.pos != 'none' & .df == 'pos') {
        out_df <- out_df[out_df$pos == .pos,]
    }
    if (.df == 'diet') {
        out_df <- out_df[!is.na(out_df$diet),]
    }
    
    out_df <- out_df[out_df$stat == .stat,]
    out_df <- out_df[,!(colnames(out_df) %in% c('stat', 'pos'))]
    # Making sure they're all sorted the same
    out_df <- dplyr::arrange(out_df, taxon, species)
    
    # Converting taxon to factor
    out_df$taxon <-  factor(out_df$taxon, levels = c('Rodent', 'Bat'))
    
    # Converting diet to factor, if present
    if (.df %in% c('spp', 'diet', 'pos')) {
        out_df$diet <-  factor(out_df$diet, levels = c('Herbivorous', 'Omnivorous', 
                                                       'Protein'))
    } else if (.df == 'clear') {
        out_df$diet <-  factor(out_df$diet, levels = c('Carb', 'Protein'))
    }
    
    # phylolm requires that the rownames match the species names
    rownames(out_df) <- out_df$species
    
    return(out_df)
}


get_tr <- function(.df) {
    
    # Just choosing parameters for get_df. These won't affect anything.
    .pos <- ifelse(.df == 'pos', 'med', 'none')

    focal_df <- get_df(.df, .pos = .pos)

    out_tr <- ape::read.tree('data/tree.nwk')
    out_tr$tip.label <- gsub('_', ' ', out_tr$tip.label)
    tips_to_rm <- out_tr$tip.label[!out_tr$tip.label %in% focal_df$species]
    out_tr <- ape::drop.tip(out_tr, tip = tips_to_rm)
    
    return(out_tr)
}
