
rm(list = ls())
gc()

set.seed(1954)

# adjust to your setting
.libPaths("~/Rlib/")
setwd("~/Desktop/Teaching/Stan-tutorial-Lueven/scripts")

library(rjson)
library(bayesplot)
library(posterior)
library(cmdstanr)
library(ggplot2)
source("tools.r")

##########################################################################
## If running on windows, this can be used for debugging purposes.
if (FALSE) {
  set_tbb(cmdstan_path, "model")
}

##########################################################################

data <- fromJSON(file = "data/linear.data.json")

# define starting distribution
init <- function() {
  list(sigma = rgamma(1, 1),
       beta = rnorm(1, mean = 1, sd = 1))
}

# transpile (translate Stan to C++ and then compile)
mod <- cmdstan_model("model/linear2.stan")

# run sampler
n_chains <- 4
fit <- mod$sample(data = data, chains = n_chains,
                  init = init,
                  save_warmup = TRUE,
                  parallel_chains = 4)

# Examine Stan's default summaries
fit$summary()

# Construct diagnostic plots
pars <- c("beta", "sigma")
bayesplot::mcmc_trace(fit$draws(inc_warmup = TRUE),
                      n_warmup = 1000, pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

# Extract posterior predictive checks
yrep <- as.matrix(
  as_draws_df(fit$draws(variables = c("y_pred"))))
head(yrep)

# We don't need the chain, iteration and draw ID, so let's remove them.
yrep <- yrep[, -(11:13)]

# Plot the posterior predictions and compare it to the real data.
bayesplot::ppc_ribbon(y = data$y, yrep = yrep, x = data$x,
                      y_draw = "point") + 
  theme_bw() +
  ylab("y")
