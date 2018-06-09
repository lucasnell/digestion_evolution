Digestion Evolution
========

Phylogenetic analyses related to the mechanistic basis of higher paracellular absorption in flying mammals
-------

Lucas A. Nell



## Description of contents

*(Alphabetical order)*

* `data`: raw data files
    - `raw_data.xlsx`: the Excel sheet exactly as it was sent to me
    - `tree.nwk`: the phylogenetic tree as downloaded from
      [timetree.org](http://timetree.org/)

* `doc`: Files to create the electronic supplementary material document.
    - `phylo_tree`: create plot of phylogeny.
    - `results`: table of coefficient estimates.


* `figs`: PDFs of all figures. These are not present in this repo, but can be
  recreated using R code in [`reports/05-plots.Rmd`](reports/05-plots.Rmd).

* `output`: output and intermediate datasets. Below, files in parentheses are where
  a given data file was created.
    - `include_mass_pos.csv`: whether to include log(mass) as a covariate for `phylolm`
      models by species and intestinal position
      ([`doc/03-include_mass`](doc/03-include_mass.md)).
    - `models_summaries.csv`: summaries for all `phylolm` and `corphylo` models
      ([`doc/04-phylo_fits`](doc/04-phylo_fits.md)).
    - `tidy_absorp.csv`: tidy version of absorption dataset
      ([`doc/02-aggregate`](doc/02-aggregate.md))
    - `tidy_clear.csv`: tidy version of clearance dataset 
      ([`doc/02-aggregate`](doc/02-aggregate.md)).
    - `tidy_pos.csv`: tidy version of morphometric dataset by species and position
      ([`doc/02-aggregate`](doc/02-aggregate.md)).
    - `tidy_spp.csv`: tidy version of morphometric dataset by species
      ([`doc/02-aggregate`](doc/02-aggregate.md)).

* `R`: R scripts with useful functions.
    - `get_data.R`: retrieve and prepare data for `phylolm` and `corphylo`.
    - `model_summaries.R`: summarize `phylolm` and `corphylo` models.

* `reports`: Documents with R code to clean, analyze, and plot data.
  Each step has an `R` or `Rmd` file paired with an `md` file.
  The former has the raw code with comments for explanation.
  The latter should be viewed for documentation.
    - `01-xl_to_dfs`: convert raw Excel file sheets to readable, "clean" data frames.
    - `02-aggregate`: convert cleaned data frames by individual to new ones by species
      and/or intestinal segment, then write to CSV files.
    - `03-include_mass`: test whether to include log(mass) as a covariate for all 
      analyses.
    - `04-phylo_fits`: `phylolm` and `corphylo` analyses.
    - `05-plots`: plots used in the main portion of the text.

