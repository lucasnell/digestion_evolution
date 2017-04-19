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

# Data frame containing 95% CI based on bootstrapping for one model
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

# Apply above function to all model fits and change the model column so it shows up
# more nicely in plots
ci_df <- mapply(mod_ci, append(nsa_fits, sef_fits), 
       rep(c('lambda', 'ou'), 2),
       rep(c('nsa', 'sef'), each = 2),
       SIMPLIFY = FALSE, USE.NAMES = FALSE) %>% 
    bind_rows %>% 
    mutate(model = recode_factor(model, lambda = "Pagel's lambda", 
                                 ou = "Ornstein-Uhlenbeck"))




nsa_plot <- ci_df %>% 
    filter(measure == 'nsa') %>% 
    mutate_if(is.numeric, funs(exp)) %>% 
    ggplot(aes(mass_log, predicted)) +
    geom_ribbon(aes(group = taxon, ymin = low, ymax = high), 
                fill = 'gray70', alpha = 0.5) +
    geom_line(aes(color = taxon), size = 0.75) +
    scale_color_brewer(NULL, palette = 'Dark2') +
    theme(legend.position = c(0.1, 0.8),
          axis.title.x = element_blank()) +
    ylab(expression('NSA (' * cm^2 * ')')) +
    facet_wrap(~ model)


sef_plot <- ci_df %>% 
    filter(measure == 'sef') %>% 
    mutate_if(is.numeric, funs(exp)) %>% 
    ggplot(aes(mass_log, predicted)) +
    geom_ribbon(aes(group = taxon, ymin = low, ymax = high), 
                fill = 'gray70', alpha = 0.5) +
    geom_line(aes(color = taxon), size = 0.75) +
    scale_color_brewer(NULL, palette = 'Dark2') +
    theme(legend.position = 'none', strip.text = element_blank(),
          axis.title.x = element_text(size = 12, margin = margin(t = 12))) +
    xlab('Mass (g)') +
    ylab('SEF') +
    facet_wrap(~ model)



pdf('phylo_plot.pdf', 6, 6, title = 'Phylogenetic Regression')
grid.newpage()
grid.draw(rbind(ggplotGrob(nsa_plot), 
                ggplotGrob(sef_plot), 
                size = "last"))
dev.off()
