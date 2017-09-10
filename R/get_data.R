# 
# Function for retrieving data of mean or standard errors among species
# 

get_data <- function(.df, .stat, 
                   .pos = c(NA, 'prox', 'med', 'dist')) {
    
    .stat  <- match.arg(.stat, c('mean', 'se'))
    .pos  <- match.arg(.pos)
    
    stopifnot(.df %in% c('spp', 'diet', 'clear', 'absorp', 'pos'))
    
    if (.df == 'pos' & is.na(.pos)) stop('provide non-NA .pos when .df == "pos"')
    if (.df != 'pos' & !is.na(.pos)) warning('.df != "pos" so .pos argument ignored')
    
    
    fn <- paste0('output/', ifelse(.df == 'diet', 'spp', .df), '_df.csv')
    
    out_df <- suppressMessages(read_csv(fn))
    out_df <- as.data.frame(out_df)
    
    
    if (!is.na(.pos) & .df == 'pos') {
        out_df <- out_df[out_df$pos == .pos,]
    }
    if (.df == 'diet') {
        out_df <- out_df[!is.na(out_df$diet),]
    }
    
    out_df <- out_df[out_df$stat == .stat,]
    out_df <- out_df[,!(colnames(out_df) %in% c('stat', 'pos'))]
    
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
}