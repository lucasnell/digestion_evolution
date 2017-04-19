suppressPackageStartupMessages({
    library(phylolm)
    library(ggplot2)
    library(dplyr)
    library(grid)
})

# Set the default ggplot theme
theme_set(theme_classic() %+replace% 
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 12)))

# Load model fits
load('./data/model_fits.RData')

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
               mass_log = .model$X[,'mass_log'],
               taxon = factor(.model$X[,'taxonRodent'], levels = c(0,1), 
                              labels = c('Bat', 'Rodent'))) %>% 
        select(model, taxon, mass_log, measure, everything())
    return(out_df)
}

# Apply above function to all model fits, change the model column so it shows up
# more nicely in plots, and add column of observed data
ci_df <- mapply(mod_ci, append(nsa_fits, sef_fits), 
       rep(c('lambda', 'ou'), 2),
       rep(c('nsa', 'sef'), each = 2),
       SIMPLIFY = FALSE, USE.NAMES = FALSE) %>% 
    bind_rows %>% 
    mutate(model = recode_factor(model, lambda = "\"Pagel's\" ~ lambda", 
                                 ou = "Ornstein-Uhlenbeck")) %>% 
    group_by(model, measure, mass_log) %>% 
    mutate(observed = sp_df$value[sp_df$measure == measure & 
                                      sp_df$mass_log == mass_log]) %>% 
    ungroup



nsa_plot <- ci_df %>%
    filter(measure == 'nsa') %>% 
    mutate_if(is.numeric, funs(exp)) %>% 
    ggplot(aes(mass_log, predicted)) +
    geom_ribbon(aes(group = taxon, ymin = low, ymax = high), 
                fill = 'gray70', alpha = 0.5) +
    geom_point(aes(y = observed, color = taxon, shape = taxon), size = 2) +
    geom_line(aes(color = taxon), size = 0.75) +
    scale_color_brewer(NULL, palette = 'Dark2') +
    scale_shape_manual(values = c(1, 2), guide = FALSE) +
    theme(legend.position = c(0.1, 0.8),
          axis.title.x = element_blank()) +
    scale_y_continuous(expression('NSA (' * cm^2 * ')'), trans = 'log',
                       breaks = c(5, 10, 20, 40)) +
    scale_x_continuous(trans = 'log', breaks = c(15, 30, 60, 120)) +
    facet_wrap(~ model, labeller = label_parsed) +
    guides(color = guide_legend(override.aes = list(shape = c(1, 2))))



sef_plot <- ci_df %>%
    filter(measure == 'sef') %>% 
    mutate_if(is.numeric, funs(exp)) %>% 
    ggplot(aes(mass_log, predicted)) +
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
    facet_wrap(~ model) +
    guides(color = guide_legend(override.aes = list(shape = c(1, 2))))


pdf('2model_plot.pdf', 6, 6, title = 'Phylogenetic Regression')
grid.newpage()
grid.draw(rbind(ggplotGrob(nsa_plot), 
                ggplotGrob(sef_plot), 
                size = "last"))
dev.off()
