
# This is to update the `running` directory
# You can ignore this file if you've only been sent this directory
purl('doc/phylo_regressions.Rmd', 'running/phylo_regressions.R')
file.copy('R', 'running', recursive = TRUE)
file.copy('data/tree.nwk', 'running/data', recursive = TRUE)
for (f in list.files('output', 'tidy_*', full.names = TRUE)) {
    file.copy(f, to = 'running/output', overwrite = TRUE)
}; rm(f)
