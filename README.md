Digestion Evolution
========

Phylogenetic analyses related to the mechanistic basis of higher paracellular absorption in flying mammals
-------

Lucas Nell



## Description of contents

*(Alphabetical order)*

* `data`: contains the raw data
    - `raw_data.xlsx`: the Excel sheet exactly as it was sent to me
    - `tree.nwk`: the phylogenetic tree as downloaded from
      [timetree.org](http://timetree.org/)

* `doc`: 
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
    - `clean_absorption.csv`: cleaned absorption data
      ([`doc/xl_to_csvs`](doc/xl_to_csvs.md))
    - `clean_clearance.csv`: cleaned clearance data
      ([`doc/xl_to_csvs`](doc/xl_to_csvs.md))
    - `clean_morph.csv`: cleaned morphometric data
      ([`doc/xl_to_csvs`](doc/xl_to_csvs.md))
    - `models_absorp.rds`: `phylolm` object for regression of fractional
      absorption ~ taxon ([`doc/phylo_regressions`](doc/phylo_regressions.md))
    - `models_pos.rds`: `phylolm` objects related to by-species-and-position model
      fits ([`doc/phylo_regressions`](doc/phylo_regressions.md))
    - `models_spp.rds`: `phylolm` objects related to by-species model fits 
      ([`doc/phylo_regressions`](doc/phylo_regressions.md))
    - `models_summaries.csv`: summaries for all `phylolm` models
      ([`doc/phylo_regressions`](doc/phylo_regressions.md))
    - `tidy_absorp.csv`: tidy version of absorption dataset
      ([`doc/tidy_csvs`](doc/tidy_csvs.md))
    - `tidy_clear.csv`: tidy version of clearance dataset 
      ([`doc/tidy_csvs`](doc/tidy_csvs.md))
    - `tidy_pos.csv`: tidy version of morphometric dataset by species and position
      ([`doc/tidy_csvs`](doc/tidy_csvs.md))
    - `tidy_spp.csv`: tidy version of morphometric dataset by species
      ([`doc/tidy_csvs`](doc/tidy_csvs.md))

* `R`: R scripts with useful functions
    - `get_data.R`: for a given [analysis set](doc/tidy_csvs.md#summary-of-output),
      retrieve a phylogeny or data frames of means or standard errors among species
    - `model_summaries.R`: summarize `phylolm` models
    - `Modified_corphylo.R`: modified version of `ape::corphylo` that outputs 
      confidence intervals using Fisher information

* `running`: files to fully run script in `phylo_regressions.Rmd`