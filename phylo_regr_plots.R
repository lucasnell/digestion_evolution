suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(magrittr)
    library(phylolm)
    library(ape)
    library(ggplot2)
    library(grid)
    library(gridExtra)
})
# Set the default ggplot theme
theme_set(theme_classic() %+replace% 
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 10),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))


load('data/spp_models.rda')






# Creates data frame containing 95% CI based on bootstrapping for one species model
spp_mod_ci <- function(.model, y_measure){
    
    if ('log_mass' %in% colnames(.model$X)) {
        new_data <- rbind(c(1.0, 1.0), c(1.0, 0.0), 
                          rep(mean(.model$X[,'log_mass']), 2))
        new_data_df <- data.frame(taxon = c(1.0, 0.0), 
                                  log_mass = rep(mean(.model$X[,'log_mass']), 2))
    } else {
        new_data <- rbind(c(1.0, 1.0), c(1.0, 0.0))
        new_data_df <- data.frame(taxon = c(1.0, 0.0))
    }
    
    # Column names coinciding with phylogenetic parameters from models lambda, BM, and
    # OUfixed
    phylo_cols <- which(colnames(.model$bootstrap) %in% c('lambda', 'sigma2', 'alpha'))
    ci_matrix <- .model$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(x[-phylo_cols] %*% new_data)
              }) %>% 
        apply(1,
              function(x) as.numeric(quantile(x, probs = c(0.025, 0.975)))) %>% 
        t
    
    out_df <- cbind(predict(.model, new_data_df), ci_matrix) %>% 
        as_data_frame %>% 
        rename_(.dots = setNames(colnames(.), c('predicted', 'low', 'high'))) %>% 
        mutate(measure = y_measure,
               taxon = factor(c(1.0, 0.0), levels = c(0,1), 
                              labels = c('Rodent', 'Bat'))) %>% 
        select(taxon, measure, everything())
    return(out_df)
}





# Same for species by position models
pos_mod_ci <- function(.model, y_measure, pos_name = NULL){
    
    
    if ('log_mass' %in% colnames(.model$X)) {
        new_data <- rbind(c(1.0, 1.0), c(1.0, 0.0), 
                          rep(mean(.model$X[,'log_mass']), 2))
        new_data_df <- data.frame(taxon = c(1.0, 0.0), 
                                  log_mass = rep(mean(.model$X[,'log_mass']), 2))
    } else {
        new_data <- rbind(c(1.0, 1.0), c(1.0, 0.0))
        new_data_df <- data.frame(taxon = c(1.0, 0.0))
    }
    
    # Column names coinciding with phylogenetic parameters from models lambda, BM, and
    # OUfixed:
    phylo_cols <- which(colnames(.model$bootstrap) %in% c('lambda', 'sigma2', 'alpha'))
    
    ci_matrix <- .model$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(x[-phylo_cols] %*% new_data)
              }) %>% 
        apply(1,
              function(x) as.numeric(quantile(x, probs = c(0.025, 0.975)))) %>% 
        t
    
    out_df <- cbind(predict(.model, new_data_df), ci_matrix) %>% 
        as_data_frame %>% 
        rename_(.dots = setNames(colnames(.), c('predicted', 'low', 'high'))) %>% 
        mutate(measure = y_measure,
               taxon = factor(c(1.0, 0.0), levels = c(0,1), 
                              labels = c('Rodent', 'Bat'))) %>% 
        select(taxon, measure, everything())
    if (!is.null(pos_name)) {
        out_df <- out_df %>% 
            mutate(pos = pos_name) %>% 
            select(taxon, pos, measure, everything())
    }
    return(out_df)
}






# =======================================================================================
# =======================================================================================
# =======================================================================================

#       By taxon only

# =======================================================================================
# =======================================================================================
# =======================================================================================


taxon_only_no_mass <- function(.model, y_name, y_axis_title, title = NULL) {
    .p <- spp_mod_ci(.model, y_name) %>%
        ggplot(aes(taxon, predicted)) +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, size = 0.5) +
        geom_segment(aes(yend = predicted, 
                         x = as.numeric(taxon) - 0.1, 
                         xend = as.numeric(taxon) + 0.1)) +
        geom_point(data = spp_df, 
                   aes_string(y = y_name, shape = 'taxon'),
                   position = position_jitter(width = 0.2, height = 0),
                   color = 'black', size = 2, fill = 'gray60') +
        scale_shape_manual(values = c(21, 1)) +
        theme(legend.position = 'none', axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(color = 'black', size = 10),
              axis.title.y = element_text(margin = margin(t = 0, r = -8, 
                                                          b = 0, l = 0))) +
        ylab(y_axis_title)
    if (!is.null(title)) .p <- .p + ggtitle(title)
    return(.p)
}



# Figure 1A
fig1a <- taxon_only_no_mass(spp_fits$int_length_mass, 'int_length_mass', 
                            expression(atop("Intestinal length / body" ~ mass^{0.4},
                                            "(" * cm / g^{0.4} * ")")),
                            '(a)')

# Figure 1B
fig1b <- taxon_only_no_mass(spp_fits$nsa_mass, 'nsa_mass', 
                            expression(atop("NSA / body" ~ mass^{0.75},
                                            "(" * cm^2 / g^{0.75} * ")")),
                            '(b)')

# Figure 4
fig4 <- taxon_only_no_mass(spp_fits$vill_area_mass, 'vill_area_mass', 
                           expression(atop("Villous surface area / body" ~ mass^{0.75},
                                           "(" * cm^2 / g^{0.75} * ")")))


# Figure 6
fig6 <- taxon_only_no_mass(spp_fits$log_total_enterocytes, 'log_total_enterocytes',
                   expression("Total enterocytes (" %*% 10^9 * ")")) +
    theme(axis.title.y = element_text(margin = margin(0, 5.5, 0, 0))) +
    scale_y_continuous(breaks = log(c(5e8, 1e9, 1.5e9)), labels = seq(0.5, 1.5, 0.5),
                       limits = log(c(1, 1.55e9))) +
    coord_trans(y = 'exp')
# Mention that bars represent model predictions for mean body mass among all species








# =======================================================================================
# =======================================================================================
# =======================================================================================

#       By taxon and position

# =======================================================================================
# =======================================================================================
# =======================================================================================


load('./data/pos_models.rda')

# Nest pos_fits by parameter, not position, bc the former is how they'll be plotted
pos_fits <- lapply(names(pos_fits$dist), 
                   function(n) {
                       list(prox = pos_fits$prox[[n]], 
                            med = pos_fits$med[[n]], 
                            dist = pos_fits$dist[[n]])
                   })
names(pos_fits) <- c('log_intestinal_diameter',
                     'villus_height', 
                     'villus_width',
                     'crypt_width',
                     'sef',
                     'enterocyte_diameter',
                     'log_enterocyte_density')

# Making data frame of confidence intervals
pos_ci <- lapply(names(pos_fits), 
                 function(n) {
                     bind_rows(list(pos_mod_ci(pos_fits[[n]]$prox, n, 'prox'),
                                    pos_mod_ci(pos_fits[[n]]$med, n, 'med'),
                                    pos_mod_ci(pos_fits[[n]]$dist, n, 'dist')))
                 }) %>% 
    bind_rows %>% 
    mutate(pos = factor(pos, levels = c('prox','med', 'dist'), 
                        labels = c('Proximal', 'Medial', 'Distal')))

# Combining the three data frame separated by position into one
pos_df <- lapply(c('prox','med', 'dist'),
                 function(p) {
                     df <- eval(parse(text = paste0(p, '_df')))
                     as_data_frame(df) %>% 
                         mutate(pos = p) %>% 
                         select(taxon, diet, species, pos, everything())
                 }) %>% 
    bind_rows %>% 
    mutate(pos = factor(pos, levels = c('prox','med', 'dist'), 
                        labels = c('Proximal', 'Medial', 'Distal')))

# Table of y-axis names for each parameter
plot_names <- read_csv('og,new
log_intestinal_diameter,"Intestinal ~ diameter ~ \'(cm)\'"
villus_height,"Villus ~ height ~ \'(mm)\'"
villus_width,"Villus ~ width ~ \'(mm)\'"
crypt_width,"Crypt ~ width ~ \'(mm)\'"
sef,"Surface ~ enlargement ~ factor ~ \'(SEF)\'"
enterocyte_diameter,"Enterocyte ~ diameter ~ \'(µm)\'"
log_enterocyte_density,"Enterocyte ~ density ~ \'(\' %*% 10^6 * \')\'"
')


# Plots for each parameter. I'm avoiding facets bc they make `ggplotGrob`s 
# annoying to combine
pos_plots <- lapply(names(pos_fits), 
       function(n) {
           plot_n <- pos_ci %>%
               filter(measure == n) %>% 
               mutate(pos = as.numeric(pos)) %>% 
               group_by(taxon) %>% 
               mutate(pos = pos + ifelse(taxon == 'Bat', 0.2, -0.2)) %>% 
               ungroup %>% 
               ggplot(aes(pos, group = taxon)) + 
               geom_point(data = pos_df %>% 
                              mutate(pos = as.numeric(pos)) %>% 
                              group_by(taxon) %>% 
                              mutate(pos = pos + ifelse(taxon == 'Bat', 0.2, -0.2)) %>% 
                              ungroup, 
                          aes_string(y = n, shape = 'taxon'),
                          color = 'black', size = 2, fill = 'gray60',
                          position = position_jitter(0.075, 0)) +
               geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, size = 0.5) +
               geom_segment(aes(y = predicted, yend = predicted,
                                x = pos - 0.1, xend = pos + 0.1), 
                            size = 0.5) +
               theme(legend.position = 'none', legend.margin = margin(0,0,0,0),
                     axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                     axis.text.x = element_blank(), legend.title = element_blank()) +
               scale_shape_manual(values = c(21, 1)) +
               ylab(eval(parse(text = plot_names[plot_names$og == n,]$new))) +
               scale_x_continuous(breaks = 1:3, 
                                  labels = c('Proximal', 'Medial', 'Distal'))
           return(plot_n)
       })
names(pos_plots) <- names(pos_fits)


# Figure 1c
fig1c <- pos_plots$log_intestinal_diameter +
    theme(legend.position = 'bottom', 
          axis.text.x = element_text(color = 'black', size = 10)) +
    ggtitle("(c)") +
    scale_y_continuous(breaks = log(seq(0.6, 1.2, 0.2)), labels = seq(0.6, 1.2, 0.2)) +
    coord_trans(y = 'exp')


# Figure 2a
fig2a <- pos_plots$villus_height +
    ggtitle('(a)')


# Figure 2b
fig2b <- pos_plots$villus_width +
    ggtitle('(b)')

# Figure 2c
fig2c <- pos_plots$crypt_width +
    theme(legend.position = 'bottom', 
          axis.text.x = element_text(color = 'black', size = 10)) +
    ggtitle("(c)")



combine_fig <- function(fig_num) {
    grob_list <- c(lapply(ls(envir = .GlobalEnv)[grepl(paste0('^fig', fig_num), 
                                                       ls(envir = .GlobalEnv))], 
                        function(n) ggplotGrob(eval(parse(text = n)))),
                   size = "first")
    grid.newpage()
    grid.draw(do.call(rbind, grob_list))
}


combine_fig(1)
# 3.875" wide, 9.4375 " tall
combine_fig(2)



# figure 3
fig3 <- pos_plots$sef +
    theme(legend.position = 'top', 
          axis.text.x = element_text(color = 'black', size = 10))

# Figure 5
fig5a <- pos_plots$enterocyte_diameter +
    scale_y_continuous(breaks = seq(2e-3, 10e-3, 2e-3), labels = seq(2, 10, 2)) +
    ggtitle('(a)')

fig5b <- pos_plots$log_enterocyte_density +
    ggtitle('(b)') +
    theme(axis.text.x = element_text(color = 'black', size = 10),
          legend.position = 'bottom')

combine_fig(5)
