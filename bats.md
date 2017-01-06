Bat phylogenetic tree
================
Lucas Nell
2017-01-06

Background
==========

The paper by [Agnarsson et al. (2011)](http://dx.doi.org/10.1371/currents.RRN1212) used only one mitochrondrial DNA gene and is the reason the Open Tree of Life (OTOL) version doesn't have a single resolved branch for *Myotis lucifugus*. A more recent paper by [Shi and Rabosky (2015)](http://onlinelibrary.wiley.com/doi/10.1111/evo.12681/abstract) uses multiple mitochondrial and nuclear genes and contains *M. lucifugus*. This paper will be used for this tree.

Preamble
========

I first need to `source` the `'tree_preamble.R'` file to load necessary R packages and make the data frame of species names.

``` r
source('tree_preamble.R')
```

Retrieving tree
===============

You can download the tree by [Shi and Rabosky (2015)](http://onlinelibrary.wiley.com/doi/10.1111/evo.12681/abstract) (OTOL link [here](https://tree.opentreeoflife.org/curator/study/view/ot_254)), adjust its tip labels to remove underscores, and drop outgroups as such:

``` r
bat_tr <- get_study(study_id = 'ot_254', object_format = 'phylo')
bat_tr$tip.label <- gsub('_', ' ', bat_tr$tip.label)
bat_tr <- drop.tip(bat_tr, tip = c('Mus musculus', 'Sorex araneus', 'Canis lupus'))
bat_tr
```


    Phylogenetic tree with 812 tips and 811 internal nodes.

    Tip labels:
        Neoromicia roseveari, Neoromicia brunneus, Neoromicia tenuipinnis, Neoromicia rendalli, Neoromicia nanus, Laephotis wintoni, ...

    Rooted; includes branch lengths.

Checking species names
======================

All bat species are present in this tree.

``` r
sp_df$species[sp_df$type == 'bat' & ! sp_df$species %in% bat_tr$tip.label]
```

    character(0)

Drawing tree
============

If you want to make a pdf of this tree that's (somewhat) readable:

``` r
pdf('bat_tree.pdf', width = 6, height = 24)
    plot(bat_tr, cex = 0.1, edge.width = 0.25, no.margin = TRUE, label.offset = 0.5)
    nodelabels(frame = 'none', cex = 0.1, adj = c(1, 0))
    tiplabels(frame = 'none', cex = 0.1, adj = c(0, 0.5))
dev.off()
```
