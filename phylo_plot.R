suppressPackageStartupMessages({
    library(phylolm)
    library(ggplot2)
    library(dplyr)
    library(grid)
})

# Set the default ggplot theme
theme_set(theme_classic() %+replace% 
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 12),
                    legend.background = element_blank()))

# Load model fits
load('./data/model_fits.RData')


# Function to return 'wide' version of an input data frame of morphometric measurements 
# (defaults to morph_df), given a character vector of measures
prep_df <- function(measures, input_df = morph_df, by_sp = TRUE, 
                    trans_fun = 'log'){
    
    # List of measures wrapped in ticks (``) so they'll play nice with dplyr and tidyr
    # as column names even if they have spaces
    meas_list <- lapply(measures, function(s) sprintf('`%s`', s))
    # Replacing all spaces with underscores for better naming in the output data frame
    meas_clean <- gsub(' ', '_', measures)
    # Names of transformed columns, also with underscores rather than spaces
    trans_meas <- paste(meas_clean, trans_fun, sep = '_')
    
    new_df <- input_df %>% 
        # Changing from tall to wide format
        spread(measure, value) %>% 
        # Selecting measurement columns, plus the identifying columns
        select_(.dots = append(list('diet', 'taxon', 'species', 'id'), 
                               meas_list)) %>% 
        # Removing all rows with all NAs in measures columns
        filter(Reduce(`+`, lapply(.[,measures], is.na)) < length(measures)) %>% 
        # Replacing spaces in measures-column names with underscores
        rename_(.dots = setNames(as.list(sprintf('`%s`', measures)), meas_clean)) %>% 
        # Doing the transformation now, before taking any means
        mutate_(.dots = setNames(as.list(sprintf('%s(%s)', trans_fun, meas_clean)), 
                                 trans_meas)) %>% 
        # Taking mean by sample
        group_by(diet, taxon, species, id) %>% 
        summarize_all(mean, na.rm = TRUE) %>% 
        ungroup
    
    if (by_sp) {
        new_df <- new_df %>% 
            # Grouping by, then taking mean of all measurement columns and transformed-
            # measurement columns
            group_by(diet, taxon, species) %>% 
            summarize_at(.vars = c(trans_meas, meas_clean), mean) %>% 
            ungroup %>% 
            arrange(taxon, diet, species) %>% 
            # To change row names, it can't be a tibble, so I'm reverting back to normal
            # data frame
            as.data.frame
        # phylolm requires that the rownames match the species names
        rownames(new_df) <- new_df$species
    } else {
        new_df <- new_df %>% select(species, everything()) %>% 
            as.data.frame
    }
    
    return(new_df)
}

# Get raw data, remove unnecessary columns and data frame
source('tidy_csv.R')
sp_df <- prep_df(measures = c('nsa', 'sef', 'mass')) %>% 
    as_data_frame %>% 
    select(-nsa, -sef, -mass) %>% 
    gather('measure', 'value', sef_log, nsa_log) %>% 
    mutate(measure = gsub('_log', '', measure))
rm(morph_df)



# Creates data frame containing 95% CI based on bootstrapping for one model
mod_ci <- function(.model, mod_name, y_measure){
    ci_matrix <- .model$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(x[1:3] %*% t(.model$X))
              }) %>% 
        apply(1,
              function(x) as.numeric(quantile(x, probs = c(0.025, 0.975)))) %>% 
        t
    
    out_df <- cbind(as.numeric(.model$fitted.values), ci_matrix) %>% 
        as_data_frame %>% 
        rename_(.dots = setNames(colnames(.), c('predicted', 'low', 'high'))) %>% 
        mutate(model = mod_name, 
               measure = y_measure,
               log_mass = .model$X[,'log_mass'],
               taxon = factor(.model$X[,'taxonRodent'], levels = c(0,1), 
                              labels = c('Bat', 'Rodent'))) %>% 
        select(model, taxon, log_mass, measure, everything())
    return(out_df)
}

# Apply above function to both Pagel's lambda model fits and add column of observed data
ci_df <- mapply(mod_ci, list(nsa_fits[['lambda']], sef_fits[['lambda']]),
       rep('lambda', 2), c('nsa', 'sef'),
       SIMPLIFY = FALSE, USE.NAMES = FALSE) %>% 
    bind_rows %>% 
    select(-model) %>% 
    group_by(measure, log_mass) %>% 
    mutate(observed = sp_df$value[sp_df$measure == measure & 
                                      sp_df$log_mass == log_mass]) %>% 
    ungroup



nsa_plot <- ci_df %>%
    filter(measure == 'nsa') %>% 
    mutate_if(is.numeric, funs(exp)) %>% 
    ggplot(aes(log_mass, predicted)) +
    geom_ribbon(aes(group = taxon, ymin = low, ymax = high), 
                fill = 'gray70', alpha = 0.5) +
    geom_point(aes(y = observed, color = taxon, shape = taxon), size = 2) +
    geom_line(aes(color = taxon), size = 0.75) +
    scale_color_brewer(NULL, palette = 'Dark2') +
    scale_shape_manual(values = c(1, 2), guide = FALSE) +
    theme(legend.position = c(0.15, 0.8),
          axis.title.x = element_blank()) +
    scale_y_continuous(expression('NSA (' * cm^2 * ')'), trans = 'log',
                       breaks = c(5, 10, 20, 40)) +
    scale_x_continuous(trans = 'log', breaks = c(15, 30, 60, 120)) +
    guides(color = guide_legend(override.aes = list(shape = c(1, 2))))



sef_plot <- ci_df %>%
    filter(measure == 'sef') %>% 
    mutate_if(is.numeric, funs(exp)) %>% 
    ggplot(aes(log_mass, predicted)) +
    geom_ribbon(aes(group = taxon, ymin = low, ymax = high), 
                fill = 'gray70', alpha = 0.5) +
    geom_point(aes(y = observed, color = taxon, shape = taxon), size = 2) +
    geom_line(aes(color = taxon), size = 0.75, linetype = 2) +
    scale_color_brewer(NULL, palette = 'Dark2') +
    scale_shape_manual(values = c(1, 2), guide = FALSE) +
    theme(legend.position = 'none', strip.text = element_blank(),
          axis.title.x = element_text(size = 12, margin = margin(t = 12))) +
    scale_y_continuous('SEF', trans = 'log',
                       breaks = c(8, 12, 18)) +
    scale_x_continuous('Mass (g)', trans = 'log', breaks = c(15, 30, 60, 120)) +
    guides(color = guide_legend(override.aes = list(shape = c(1, 2))))


pdf('phylo_plot.pdf', 4, 6, title = 'Phylogenetic Regression', useDingbats = FALSE)
grid.newpage()
grid.draw(rbind(ggplotGrob(nsa_plot), 
                ggplotGrob(sef_plot), 
                size = "last"))
dev.off()
