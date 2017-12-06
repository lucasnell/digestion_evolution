Digestion Evolution
========

Phylogenetic analyses related to the mechanistic basis of higher paracellular absorption in flying mammals
-------

Lucas A. Nell



## Description of contents

*(Alphabetical order)*

* `corphyloCpp`  __NOT UPDATED__:

* `data`: contains the raw data
    - `raw_data.xlsx`: the Excel sheet exactly as it was sent to me
    - `tree.nwk`: the phylogenetic tree as downloaded from
      [timetree.org](http://timetree.org/)

* `doc` __NOT UPDATED__: 
    - `plots`: plots from analyses.
    - `phylo_regressions`: `phylolm` and `corphylo` analyses
    - `phylo_tree`: create pdf of phylogeny in `figs/`
    - `results`: table of coefficient estimates with confidence intervals.
    - `tidy_csvs`: convert cleaned CSVs to tidy ones that are directly used in
      analyses and plots
    - `xl_to_csvs`: convert raw Excel file sheets to readable, "clean" CSV files

* `figs`: PDFs of all figures. These are not present in this repo, but can be
  recreated using files from `doc`.

* `output`: output and intermediate datasets. Below, files in parentheses are where
  a given data file was created.
    - `include_mass_pos.csv`: whether to include log(mass) as a covariate for `phylolm`
      models by species and intestinal position
      ([`doc/03-include_mass`](doc/03-include_mass.md))
    - `models_summaries.csv`: summaries for all `phylolm` and `corphylo` models
      ([`doc/04-phylo_regressions`](doc/04-phylo_regressions.md))
    - `tidy_absorp.csv`: tidy version of absorption dataset
      ([`doc/02-aggregate`](doc/02-aggregate.md))
    - `tidy_clear.csv`: tidy version of clearance dataset 
      ([`doc/02-aggregate`](doc/02-aggregate.md))
    - `tidy_pos.csv`: tidy version of morphometric dataset by species and position
      ([`doc/02-aggregate`](doc/02-aggregate.md))
    - `tidy_spp.csv`: tidy version of morphometric dataset by species
      ([`doc/02-aggregate`](doc/02-aggregate.md))

* `R`: R scripts with useful functions
    - `get_data.R`: retrieve and prepare data for `phylolm` and `corphylo`
    - `model_summaries.R`: summarize `phylolm` and `corphylo` models
