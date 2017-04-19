library(phylolm)
library(ggplot2)
library(dplyr)


# Load model fits
load('./data/model_fits.RData')


model <- nsa_fits[['lambda']]
# # Approximate 95% CI based on estimated standard errors from the model
vc <- vcov(model)
mm <- model.matrix(~1 + mass_log + taxon, data = sp_df)
p.se <- apply(mm, 1, function(x) sqrt(t(x) %*% vc %*% x))
p.ci <- p.se * 1.96


# Now doing it based on bootstrapping
model$bootstrap %>% 
    apply(1, 
          function(x) {
              matrix(x[1:3] %*% t(model$X))
          }) %>% 
    apply(1,
          function(x) quantile(x, probs = c(0.025, 0.975)))
