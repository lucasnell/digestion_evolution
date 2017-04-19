Digestion Evolution
========

Phylogenetic analyses related to the mechanistic basis of higher paracellular absorption in flying mammals
-------

Lucas Nell

The `phylo_regr.R` and `phylo_regr.md` files contains the analyses; the former has
the raw R code, and the latter has a compiled description of the steps.

`tidy_csv.R` and `xl_to_csv.R` contain code to prepare the initial Excel file for the 
analyses.

`model_fits.RData` contains the phylogenetic regression model fits.
You can load this file yourself by using the following:
```{r}
library(phylolm)
load('model_fits.RData')
```

`tree.nwk` is the phylogenetic tree from [timetree.org](http://timetree.org/).

