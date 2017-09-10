get_df <- function(.df, .stat, 
                   .pos = c(NA, "prox", "med", "dist")) {
    
    .stat  <- match.arg(.stat, c("mean", "se"))
    .pos  <- match.arg(.pos)
    
    stopifnot(.df %in% c("spp", "diet", "clear", "absorp", "pos"))
    
    if (.df == "pos" & is.na(.pos)) stop("provide non-NA .pos when .df == 'pos'")
    if (.df != "pos" & !is.na(.pos)) warning(".df != 'pos' so .pos argument ignored")
    
    
    fn <- paste0('output/', ifelse(.df == "diet", "spp", .df), '_df.csv')
    
    out_df <- suppressMessages(read_csv(fn))
    out_df <- as.data.frame(out_df)
    
    
    # spp_df
    # %>%
    #     # Converting taxon and diet to factors
    #     mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
    #            diet = factor(diet, levels = c("Herbivorous", "Omnivorous", "Protein")))
    
    
    # pos_df
    # %>%
    #     # Converting taxon and diet to factors
    #     mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
    #            diet = factor(diet, levels = c("Herbivorous", "Omnivorous", "Protein")))
    
    # clear_df
    # %>% 
    #     mutate(taxon = factor(taxon, levels = c('Rodent', 'Bat')),
    #            diet = factor(diet, levels = c('Carb', 'Protein')))
    
    if (.df =="diet") {
        out_df <- out_df[!is.na(out_df$diet),]
    }
    if (!is.na(.pos) & .df == "pos") {
        out_df <- out_df[out_df$pos == .pos,]
    }
    
    out_df <- out_df[out_df$stat == .stat,]
    out_df <- out_df[,!(colnames(out_df) %in% c("stat", "pos"))]
    
    out_df$taxon <-  factor(out_df$taxon, levels = c('Rodent', 'Bat'))
    
    # phylolm requires that the rownames match the species names
    rownames(out_df) <- out_df$species
}