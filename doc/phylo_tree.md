Visualizing phylogenetic tree
================
Lucas Nell
10 Sep 2017

Visualizing tree
================

Here is how I created the phylogenetic tree.

``` r
# Commented areas are for when I was basing tip points on `sp_df` data frame values
x_end <- 32
gg_tr <- ggtree(tr)
gg_tr$data$x <- gg_tr$data$x - max(gg_tr$data$x)
phylo <- gg_tr + #%<+% {sp_df %>% select(species, everything())} +
    theme_tree2(axis.title.x = element_text(size = 14),
                legend.position = c(0.25, 0.75), legend.box = 'horizontal',
                legend.background = element_rect(color = NA, fill = NA)) +
    geom_tiplab(aes(x = x + 3), size = 2, fontface = 'bold.italic') +
    geom_tippoint(aes(x = x)) +  # , size = log_mass, color = diet)) +
    geom_text(data = data_frame(x = rep(x_end + 1, 2), y = c(14, 5), 
                                label = c('Rodents', 'Bats')), 
              aes(label = label), angle = -90, vjust = 0, size = 6) +
    geom_segment(data = data_frame(x = rep(x_end, 2), xend = rep(x_end, 2), 
                                   y = c(1, 10), yend = c(9, 18)), 
              aes(xend = xend, yend = yend)) +
    scale_x_continuous('Time (mya)', limits = c(-100, x_end + 2),
                       breaks = seq(-100, 0, 25), labels = seq(100, 0, -25))# +
    # scale_color_manual('Diet', values = c('#e66101','#fdb863','#b2abd2','#5e3c99'),
    #                      guide = guide_legend(direction = "vertical",
    #                                             title.position = 'top', 
    #                                             title.hjust = 0.5)) +
    # scale_size_continuous('log(mass)', 
    #                      guide = guide_legend(direction = "horizontal",
    #                                           title.position = 'top', 
    #                                           title.hjust = 0.5))
```

Saving this figure.

``` r
ggsave('figs/phylo.pdf', phylo, width = 6, height = 4)
```
