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


load('./data/spp_models.rda')

load('./data/pos_models.rda')






# Creates data frame containing 95% CI based on bootstrapping for one species model
spp_mod_ci <- function(.model, y_measure, mod_name = NULL){
    ci_matrix <- .model$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(head(x, -2) %*% t(.model$X))
              }) %>% 
        apply(1,
              function(x) as.numeric(quantile(x, probs = c(0.025, 0.975)))) %>% 
        t
    
    out_df <- cbind(as.numeric(.model$fitted.values), ci_matrix) %>% 
        as_data_frame %>% 
        rename_(.dots = setNames(colnames(.), c('predicted', 'low', 'high'))) %>% 
        mutate(measure = y_measure,
               log_mass = tryCatch(.model$X[,'log_mass'], 
                                   error = function(e) rep(NA, nrow(.model$X))),
               taxon = factor(.model$X[,'taxonBat'], levels = c(0,1), 
                              labels = c('Rodent', 'Bat'))) %>% 
        select(taxon, log_mass, measure, everything())
    if (!is.null(mod_name)) {
        out_df <- out_df %>% 
            mutate(model = mod_name) %>% 
            select(model, taxon, log_mass, measure, everything())
    }
    # Check if there is no need for log_mass or many of the rows
    if (ncol(.model$bootstrap) == 4) {
        out_df <- out_df %>% 
            select(-log_mass) %>% 
            distinct(taxon, measure, predicted, low, high)
    }
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
    ci_matrix <- .model$bootstrap %>% 
        apply(1, 
              function(x) {
                  matrix(head(x, -2) %*% new_data)
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



# Figure 1A
fig1a <- spp_mod_ci(spp_fits[[1]], 'int_length_mass') %>%
    ggplot(aes(taxon, predicted)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, size = 0.5) +
    geom_segment(aes(yend = predicted, 
                     x = as.numeric(taxon) - 0.05, 
                     xend = as.numeric(taxon) + 0.05)) +
    geom_point(data = spp_df, 
               aes(y = int_length_mass, shape = taxon),
               position = position_jitter(width = 0.2, height = 0),
               color = 'black', size = 2, fill = 'gray60') +
    # geom_text(data = NULL, label = 'A', x = 0.5, y = 11.7, size = 6, 
    #           vjust = 1, hjust = 0) +
    scale_shape_manual(values = c(21, 1)) +
    ggtitle('(a)') +
    theme(legend.position = 'none', axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0))) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(expression(atop("Intestinal length / body" ~ mass^{0.4},
                                       "(" * cm / g^{0.4} * ")")))

# Figure 1B
fig1b <- spp_mod_ci(spp_fits[[2]], 'nsa_mass') %>%
    ggplot(aes(taxon, predicted)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, size = 0.5) +
    geom_segment(aes(yend = predicted, 
                     x = as.numeric(taxon) - 0.05, 
                     xend = as.numeric(taxon) + 0.05)) +
    geom_point(data = spp_df, 
               aes(y = nsa_mass, shape = taxon),
               position = position_jitter(width = 0.2, height = 0),
               color = 'black', size = 2, fill = 'gray60') +
    # geom_text(data = NULL, label = 'B', x = 0.5, y = 2.2, size = 6, 
    #           vjust = 1, hjust = 0) +
    ggtitle('(b)') +
    scale_shape_manual(values = c(21, 1)) +
    theme(legend.position = 'none', axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = -8, b = 0, l = 0))) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(expression(atop("NSA / body" ~ mass^{0.75},
                                       "(" * cm^2 / g^{0.75} * ")")))


# Figure 4
fig4 <- spp_mod_ci(spp_fits[[3]], 'vill_area_mass') %>%
    ggplot(aes(taxon, predicted)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, size = 0.5) +
    geom_segment(aes(yend = predicted, 
                     x = as.numeric(taxon) - 0.05, 
                     xend = as.numeric(taxon) + 0.05)) +
    geom_point(data = spp_df, 
               aes(y = vill_area_mass, shape = taxon),
               position = position_jitter(width = 0.2, height = 0),
               color = 'black', size = 2, fill = 'gray60') +
    scale_shape_manual(values = c(21, 1)) +
    theme(legend.position = 'none', axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = -8, b = 0, l = 0))) +
    scale_y_continuous(expression(atop("Villous surface area / body" ~ mass^{0.75},
                                       "(" * cm^2 / g^{0.75} * ")")))



# Figure 6
fig6 <- spp_mod_ci(spp_fits[[4]], 'log_total_enterocytes') %>%
    ggplot(aes(log_mass, predicted)) + 
    geom_ribbon(aes(group = taxon, ymin = low, ymax = high), 
                fill = 'gray70', alpha = 0.5) +
    geom_point(data = spp_df, aes(y = log_total_enterocytes, shape = taxon),
               color = 'black', size = 2, fill = 'gray60') +
    geom_line(aes(linetype = taxon)) +
    theme(legend.position = c(0.15, 0.9), legend.title = element_blank(),
          legend.key.width = unit(0.05, "npc")) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_shape_manual(values = c(21, 1)) +
    scale_y_continuous(expression("Total enterocytes (" %*% 10^9 * ")"), 
                       breaks = log(seq(0, 2e9, 5e8)), labels = seq(0, 2, 0.5)) +
    scale_x_continuous("Body mass (g)", breaks = log(seq(0, 150, 50)), 
                       labels = seq(0, 150, 50)) +
    coord_trans(x = 'exp', y = 'exp')










# =======================================================================================
# =======================================================================================
# =======================================================================================

#       By taxon and position

# =======================================================================================
# =======================================================================================
# =======================================================================================


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

pos_ci <- lapply(names(pos_fits), 
                 function(n) {
                     bind_rows(list(pos_mod_ci(pos_fits[[n]]$prox, n, 'prox'),
                                    pos_mod_ci(pos_fits[[n]]$med, n, 'med'),
                                    pos_mod_ci(pos_fits[[n]]$dist, n, 'dist')))
                 }) %>% 
    bind_rows %>% 
    mutate(pos = factor(pos, levels = c('prox','med', 'dist'), 
                        labels = c('Proximal', 'Medial', 'Distal')))

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

plot_names <- read_csv('og,new
log_intestinal_diameter,"Intestinal ~ diameter ~ \'(cm)\'"
villus_height,"Villus ~ height ~ \'(mm)\'"
villus_width,"Villus ~ width ~ \'(mm)\'"
crypt_width,"Crypt ~ width ~ \'(mm)\'"
sef,"Surface ~ enlargement ~ factor ~ \'(SEF)\'"
enterocyte_diameter,"Enterocyte ~ diameter ~ \'(µm)\'"
log_enterocyte_density,"Enterocyte ~ density ~ \'(\' %*% 10^6 * \')\'"
')


pos_plots <- lapply(names(pos_fits), 
       function(n) {
           plot_n <- pos_ci %>% 
               filter(measure == n) %>% 
               ggplot(aes(taxon)) + 
               geom_point(data = pos_df, aes_string(y = n, shape = 'taxon'),
                          color = 'black', size = 2, fill = 'gray60',
                          position = position_jitter(width = 0.2, height = 0)) +
               geom_errorbar(aes(ymin = low, ymax = high), width = 0.3) +
               geom_segment(aes(y = predicted, yend = predicted, 
                                x = as.numeric(taxon) - 0.2, 
                                xend = as.numeric(taxon) + 0.2)) +
               theme(legend.position = 'none', axis.title.x = element_blank()) +
               scale_shape_manual(values = c(21, 1)) +
               facet_wrap(~ pos, nrow = 1) +
               ylab(eval(parse(text = plot_names[plot_names$og == n,]$new)))
           return(plot_n)
       })
names(pos_plots) <- names(pos_fits)



n = names(pos_fits)[1]
fig1c <- pos_ci %>%
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
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.3) +
    geom_segment(aes(y = predicted, yend = predicted,
                     x = pos - 0.2,
                     xend = pos + 0.2)) +
    theme(legend.position = 'bottom', legend.margin = margin(0,0,0,0),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_text(color = 'black', size = 10)) +
    scale_shape_manual(values = c(21, 1)) +
    ylab(eval(parse(text = plot_names[plot_names$og == n,]$new))) +
    ggtitle("(c)") +
    # geom_text(data = data_frame(pos = 1:3, taxon = 0.5, y = 0.37 + 0.05), aes(y=y),
    #           label = c('Proximal', 'Medial', 'Distal'), 
    #           size = 4, vjust = 1, hjust = 0.5) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(limits = c(-0.96, 0.37 + 0.1), expand = c(0, 0),
                       breaks = log(seq(0.6, 1.2, 0.2)), labels = seq(0.6, 1.2, 0.2)) +
    scale_x_continuous(# breaks = c(1 + c(-0.2, 0.2), 2 + c(-0.2, 0.2), 
                       #            3 + c(-0.2, 0.2)), 
                       # labels = rep(c('Rodent', 'Bat'), 3)
        breaks = 1:3, labels = c('Proximal', 'Medial', 'Distal')) + 
    guides(shape = guide_legend(title = NULL))

#




# Figure 1c
fig1c <- pos_plots$log_intestinal_diameter +
    geom_text(data = data_frame(pos = 'Distal', taxon = 0.5, y = 0.37), aes(y=y),
              label = 'C', size = 6, vjust = 1, hjust = 0) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(limits = c(-0.96, 0.37), expand = c(0, 0),
                       breaks = log(seq(0.6, 1.2, 0.2)), labels = seq(0.6, 1.2, 0.2)) +
    coord_trans(y = 'exp') +
    coord_cartesian(xlim = c(0.4, 2.6), expand = FALSE)


# Figure 2a
fig2a <- pos_plots$villus_height +
    theme(axis.text.x = element_blank()) +
    geom_text(data = data_frame(pos = 'Distal', taxon = 0.5, y = 0.922), aes(y=y),
              label = 'A', size = 6, vjust = 1, hjust = 0) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(limits = c(0.14, 0.922), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.4, 2.6), expand = FALSE)


# Figure 2b
fig2b <- pos_plots$villus_width +
    theme(axis.text.x = element_blank(), strip.text = element_blank()) +
    geom_text(data = data_frame(pos = 'Distal', taxon = 0.5, y = 0.13), aes(y=y),
              label = 'B', size = 6, vjust = 1, hjust = 0) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(limits = c(0.036, 0.13), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.4, 2.6), expand = FALSE)


# Figure 2c
fig2c <- pos_plots$crypt_width +
    theme(strip.text = element_blank()) +
    geom_text(data = data_frame(pos = 'Distal', taxon = 0.5, y = 0.059), aes(y=y),
              label = 'C', size = 6, vjust = 1, hjust = 0) +
    # Used ggplot_build(<obj>)$layout$panel_ranges[[1]]$y.range to get y range
    scale_y_continuous(limits = c(0.014, 0.059), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.4, 2.6), expand = FALSE)


# fig1 <- function() {
    g1 <- ggplotGrob(fig1a)
    g2 <- ggplotGrob(fig1b)
    g3 <- ggplotGrob(fig1c)
    # colnames(g1) <- paste0(seq_len(ncol(g1)))
    # colnames(g2) <- paste0(seq_len(ncol(g2)))
    # colnames(g3) <- paste0(seq_len(ncol(g3)))
    # grid.draw(combine(g1, g2, g3, along=2))
    # 
    # newWidth = unit.pmax(g1$widths[2:3], g2$widths[2:3], g3$widths[2:3])
    # 
    # g1$widths[2:3] = as.list(newWidth)
    # g2$widths[2:3] = as.list(newWidth)
    # g3$widths[2:3] = as.list(newWidth)
    
    grid.arrange(g1, g2, g3, ncol=1)
    
    grid.newpage()
    grid.draw(rbind(ggplotGrob(fig1a),
                          ggplotGrob(fig1b),
                    ggplotGrob(fig1c),
                    size = "first"))
    
    
    
# }

fig1()


fig2 <- function() {
    g1 <- ggplotGrob(fig2a)
    g2 <- ggplotGrob(fig2b)
    g3 <- ggplotGrob(fig2c)
    colnames(g1) <- paste0(seq_len(ncol(g1)))
    colnames(g2) <- paste0(seq_len(ncol(g2)))
    colnames(g3) <- paste0(seq_len(ncol(g3)))
    grid.draw(combine(g1, g2, g3, along=2))
}

fig2()
# 3.875" wide, 9.4375 " tall

