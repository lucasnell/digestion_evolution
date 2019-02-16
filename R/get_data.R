# 
# Functions to retrieve and prepare data for phylolm and corphylo
# 



#' Get data frame for an analysis set
#'
#' @param .df refers to the analysis set and must be one of the following:
#' \describe{
#'     \item{ \code{'spp'} }{ morphometric analyses by species }
#'     \item{ \code{'clear'} }{ analyses for clearance }
#'     \item{ \code{'absorp'} }{ analyses for absorption }
#'     \item{ \code{'pos'} }{ morphometric analyses by species and position }
#'     }
#' @param .pos Position of intestine to return data from. Takes values of \code{'prox'}, 
#'     \code{'mid'}, or {'dist'}. This argument is ignored if \code{.df != 'pos'}.
#'     Defaults to \code{NA}, which will return an error if \code{.df == 'pos'}.
#' @param .stat Statistic to return, which can be \code{'mean'} or \code{'se'} for 
#'     mean and standard error, respectively. Defaults to \code{'mean'}.
#'
#' @return A data frame with factors consistently formatted and species names as 
#'     row names
#' 
#' @export
#'
#'
get_df <- function(.df, .pos = NA, .stat = c('mean', 'se')) {
    
    .df <- match.arg(.df, c('spp', 'clear', 'absorp', 'pos'))
    if (!is.na(.pos)) .pos  <- match.arg(.pos, c('prox', 'mid', 'dist'))
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



#' Get a phylogenetic tree for an analysis set
#'
#' @param .df See \code{\link{get_df}}.
#'
#' @return A phylogenetic tree for a given analysis set
#' 
#' @export
#'
#'
get_tr <- function(.df) {
    
    .df <- match.arg(.df, c('spp', 'clear', 'absorp', 'pos'))
    
    # Just choosing parameters for get_df. These won't affect anything bc this is just
    # to retrieve species names.
    .pos <- ifelse(.df == 'pos', 'mid', NA)

    focal_df <- get_df(.df, .pos = .pos)

    out_tr <- ape::read.tree('data/tree.nwk')
    out_tr$tip.label <- gsub('_', ' ', out_tr$tip.label)
    tips_to_rm <- out_tr$tip.label[!out_tr$tip.label %in% focal_df$species]
    out_tr <- ape::drop.tip(out_tr, tip = tips_to_rm)
    
    return(out_tr)
}


#' Filter a phylogenetic tree
#'
#' @param tr A \code{phylo} object.
#' @param keep_spp A character vector of species names to keep.
#'
#' @return A \code{phylo} object with only the species specified.
#' 
#' @export
#'
#'
filter_tr <- function(tr, keep_spp) {
    cleared_tr <- ape::drop.tip(
        tr, 
        tip = tr$tip.label[!tr$tip.label %in% keep_spp]
    )
    return(cleared_tr)
}




#' Create a matrix for a call to corphylo
#'
#' @param .df A data frame. It should contain the species names and parameters of 
#'     interest.
#' @param .pars A character vector. Column names for the parameters of interest.
#'
#' @return A matrix with proper rownames for corphylo to use.
#' 
#' @export
#'
cp_mat <- function(.df, .pars) {
    out <- as.matrix(.df[,.pars])
    colnames(out) <- NULL
    rownames(out) <- .df$species
    return(out)
}

