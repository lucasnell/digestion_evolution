Digestion Evolution
========

Phylogenetic analyses related to the mechanistic basis of higher paracellular absorption in flying mammals
-------

Lucas Nell

The `phylo_regr.R` and `phylo_regr.md` files contains the analyses; the former has
the raw R code, and the latter has a compiled description of the steps.

`tidy_csv.R` and `xl_to_csv.R` contain code to prepare the initial Excel file for the 
analyses.

The `data` folder contains the data:

- `clean_data.csv` is the cleaned data file of morphometrics measurements. This is the 
  file that should be used for analysis.
- `model_fits.RData` contains the phylogenetic regression model fits.
- `raw_data.xlsx` is the Excel file before any processing. It is not the file that should
  be used for analysis.
- `tree.nwk` is the phylogenetic tree from [timetree.org](http://timetree.org/).
