
rm(list = ls())
gc()

set.seed(1954)

# adjust to your setting
.libPaths("~/Rlib/")
setwd("~/Desktop/Teaching/Stan-tutorial-Lueven/scripts")

library(ggplot2)
library(bayesplot)
library(posterior)
library(cmdstanr)
library(rjson)

cmdstan_make_local(dir = cmdstan_path(),
                   cpp_options = list(
  "CXXFLAGS += -Wno-deprecated-declarations -Wno-ignored-attributes"
                   ), 
                   append = TRUE)

mc.cores = 4
iter_warmup <- 1000
iter_sampling <- 1000

#####################################################################
## 8 schools

n_schools <- 8
y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)

stan_data <- list(n_schools = n_schools,
                  y = y,
                  sigma = sigma)

mod <- cmdstan_model("model/8schools.stan")
# other models: 8schools_nc, 8schools_marginal

fit <- mod$sample(data = stan_data,
                  chains = 4, parallel_chains = 4,
                  iter_warmup = iter_warmup,
                  iter_sampling = iter_warmup,
                  # adapt_delta = 0.8,
                  seed = 1234)

fit$summary()

mcmc_draws <- fit$draws()

# plot log tau against theta to identify divergent transitions
np <- nuts_params(fit)
mcmc_scatter(mcmc_draws, pars = c("theta[1]", "tau"),
             transform = list(tau = "log"), np = np)

mcmc_pairs(mcmc_draws, pars = "tau", regex_pars = "theta\\[[1,8]\\]",
           transform = list(tau = "log"), np = np)


#####################################################################
## Disease map of Finland

data <- fromJSON(file = "data/disease_100.json")

mod <- cmdstan_model("model/disease_map.stan")

fit <- mod$sample(data = data,
                  chains = 4, parallel_chains = 4,
                  iter_warmup = iter_warmup,
                  iter_sampling = iter_warmup,
                  adapt_delta = 0.95,
                  seed = 1234)

pars <- c("alpha", "rho", "theta[1]")
fit$summary(variables = pars)
