
rm(list = ls())
gc()

set.seed(1954)

# adjust to your setting
.libPaths("~/Rlib/")
setwd("~/Desktop/Teaching/Stan-tutorial-Lueven/scripts")


library(outbreaks)
library(ggplot2)

library(rjson)
library(bayesplot)
library(posterior)
library(cmdstanr)
library(loo)
library(tidyverse)

source("tools_is.r")

mc.cores = 4

#####################################################################
## Prep data

# plot the data, obtained from outbreaks
theme_set(theme_bw())
ggplot(data = influenza_england_1978_school) + 
  geom_point(mapping = aes(x = date, y = in_bed)) + 
  labs(y = "Number of students in bed")

# create a data list to be passed to Stan
cases <- influenza_england_1978_school$in_bed
N <- 763;
n_days <- length(cases) 
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions
i0 <- 1
s0 <- N - i0
r0 <- 0
y0 = c(S = s0, I = i0, R = r0)

data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t,
                 N = N, cases = cases)

#####################################################################
## Run and fit Stan

# define starting distribution
init <- function() {
  list(beta = abs(rnorm(1, mean = 2, sd = 1)),
       gamma = abs(rnorm(1, mean = 0.4, sd = 0.5)),
       phi_inv = rexp(1, rate = 5))
}

# transpile (translate Stan to C++ and then compile)
mod <- cmdstan_model("model/sir.stan")

n_chains <- 4
fit <- mod$sample(data = data_sir,
                  chains = n_chains,
                  parallel_chains = n_chains,
                  init = init,
                  save_warmup = TRUE)


#####################################################################
## Check the inference
pars <- c("gamma", "beta", "phi", "R0")
fit$summary(variables = pars)

bayesplot::mcmc_trace(fit$draws(inc_warmup = TRUE),
                      n_warmup = 1000, pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

# Extract posterior predictive checks
pred_cases <- as.matrix(
  as_draws_df(fit$draws(variables = c("pred_cases"))))[, -(15:17)]

bayesplot::ppc_ribbon(y = data_sir$cases, yrep = pred_cases, 
                      x = data_sir$ts, y_draw = "point") + 
  theme_bw() +
  ylab("cases") + xlab("days")

#####################################################################
## Run same model with a Poisson likelihood

mod <- cmdstan_model("model/sir_poisson.stan")

fit_poisson <- mod$sample(data = data_sir,
                          chains = n_chains,
                          parallel_chains = n_chains,
                          init = init,
                          save_warmup = TRUE)

fit_poisson$summary(variables = pars)

pred_cases_poisson <- as.matrix(
  as_draws_df(fit_poisson$draws(variables = c("pred_cases"))))[, -(15:17)]

bayesplot::ppc_ribbon(y = data_sir$cases, yrep = pred_cases_poisson, 
                      x = data_sir$ts, y_draw = "point") + 
  theme_bw() +
  ylab("cases") + xlab("days")

#####################################################################
# compute PSIS-loo estimate

log_lik_draws <- fit$draws("log_lik")
loo_estimate <- loo(log_lik_draws, r_eff = relative_eff(log_lik_draws))


log_lik_draws_poisson <- fit_poisson$draws("log_lik")
loo_estimate_poisson <-
  loo(log_lik_draws_poisson, r_eff = relative_eff(log_lik_draws_poisson))

print(loo_estimate_poisson)
print(loo_estimate)


#####################################################################
## run SIR model with custom tolerance for the ODE solver

mod <- cmdstan_model("model/sir_tol.stan")

tol <- 5e-1
data_sir$tol <- tol
fit_tol <- mod$sample(data = data_sir,
                      chains = n_chains,
                      parallel_chains = n_chains,
                      init = init,
                      save_warmup = TRUE)

fit_tol$time()

log_lik <- fit_tol$draws("log_lik")

fit_tol$summary(variables = pars)
loo_estimate <- loo(log_lik, r_eff = relative_eff(log_lik))
print(loo_estimate)

log_ratios <- fit_tol$draws("log_ratios")

psis_fit <- psis(log_ratios, r_eff =  relative_eff(log_ratios))
psis_fit

# Correct Monte Carlo samplers, using importance weights.
# Only works if the log ratios don't go to 0!
is_summary(fit_tol, pars, psis_fit, log_ratios)


