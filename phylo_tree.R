#' ---
#' title: "Visualizing phylogenetic tree"
#' author: "Lucas Nell"
#' date: "`r format(Sys.Date(), '%d %b %Y')`"
#' output: github_document
#' ---
#' 
#+ setup, include = FALSE, cache = FALSE
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages({
    library(ape)
    library(ggplot2)
    library(ggtree)
})
source('tidy_csv.R')
sp_df <- prep_df(measures = c('nsa', 'sef', 'mass'))
tr <- read.tree('./data/tree.nwk')
tr$tip.label <- gsub('_', ' ', tr$tip.label)
tr <- drop.tip(tr, tip = tr$tip.label[!tr$tip.label %in% (morph_df$species %>% unique)])
#' 
#' 
#' 
#' 
#' # Visualizing tree
#' 
#' Here is the phylogenetic tree with `log(NSA)` as tip color and `log(SEF)` as tip size.
#' 
#' 
#+ phylo_plot, fig.width=8, fig.height=6, echo = FALSE
x_end <- 32
gg_tr <- ggtree(tr)
gg_tr$data$x <- gg_tr$data$x - max(gg_tr$data$x)
gg_tr %<+% {sp_df %>% select(species, everything())} +
    theme_tree2(axis.title.x = element_text(size = 14),
                legend.position = c(0.25, 0.75), legend.box = 'horizontal',
                legend.background = element_rect(color = NA, fill = NA)) +
    geom_tiplab(aes(x = x + 3), size = 3, fontface = 'bold.italic') +
    geom_tippoint(aes(x = x + 1, size = sef_log, color = nsa_log)) +
    geom_text(data = data_frame(x = rep(x_end + 1, 2), y = c(14, 5), 
                                label = c('Rodents', 'Bats')), 
              aes(label = label), angle = -90, vjust = 0, size = 6) +
    geom_segment(data = data_frame(x = rep(x_end, 2), xend = rep(x_end, 2), 
                                   y = c(1, 10), yend = c(9, 18)), 
              aes(xend = xend, yend = yend)) +
    scale_x_continuous('Time (mya)', limits = c(-100, x_end + 2),
                       breaks = seq(-100, 0, 25), labels = seq(100, 0, -25)) +
    scale_color_gradient('log(NSA)', low = 'darkblue', high = 'cadetblue1',
                         guide = guide_colorbar(direction = "horizontal",
                                                title.position = 'top', 
                                                title.hjust = 0.5)) +
    scale_size_continuous('log(SEF)', 
                         guide = guide_legend(direction = "horizontal",
                                              title.position = 'top', 
                                              title.hjust = 0.5))
