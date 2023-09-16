
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
library(loo)

source("tools_is.r")

cmdstan_path <- "~/Desktop/Teaching/Torsten/cmdstan/"
set_cmdstan_path(cmdstan_path)
bayesplot::color_scheme_set("mix-blue-green")

data <- fromJSON(file = "data/twoCpt.data.json")

# Plot drug concentration profile
p <- ggplot(data = data.frame(time = data$time[-1],
                              cObs = data$cObs),
            aes(x = time, y = cObs)) +
  geom_point() +
  theme_bw()
p

# Draw initial conditions from the prior
init <- function() {
  list(CL = exp(rnorm(1, log(10), 0.25)),
       Q = exp(rnorm(1, log(15), 0.5)),
       VC = exp(rnorm(1, log(35), 0.25)),
       VP = exp(rnorm(1, log(105), 0.5)),
       ka = exp(rnorm(1, log(2), 1)),
       sigma = abs(rnorm(1, 0, 1)))
}

model_name <- "twoCpt_ode"
mod <- cmdstan_model(paste0("model/", model_name, ".stan"))

n_chains <- 4
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = n_chains,
                  iter_warmup = 500, iter_sampling = 500)

# Save fit (useful for expensive models!)
fit$save_object(paste0("deliv/", model_name, ".fit.RDS"))


fit$time()

# check inference
pars = c("CL", "Q", "VC", "VP", "ka", "sigma")
fit$summary(variables = pars)

bayesplot::mcmc_trace(fit$draws(), pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

# posterior predictive checks
yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))
  ))[, -(52:54)]

yobs <- data$cObs
time <- data$time[-1]

p <- bayesplot::ppc_ribbon(y = yobs, yrep = yrep, 
                           x = time, y_draw = "point") +
  xlab("time(h)") + ylab("drug plasma concentration (mg/L)") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18))
p


# PSIS check for ODE tolerance
log_ratios <- fit$draws("log_ratios")

psis_fit <- psis(log_ratios, r_eff =  relative_eff(log_ratios))
psis_fit


#####################################################################
## Use built-in analytical solution

# two compartment model
mod <- cmdstan_model("model/twoCpt.stan")
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = n_chains,
                  iter_warmup = 500, iter_sampling = 500)

fit$time()
fit$summary()

yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))
  ))[, -(52:54)]

yobs <- data$cObs
time <- data$time[-1]

p <- bayesplot::ppc_ribbon(y = yobs, yrep = yrep, 
                           x = time, y_draw = "point") +
  xlab("time(h)") + ylab("drug plasma concentration (mg/L)") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18))
p

# one compartment model
mod <- cmdstan_model("model/oneCpt.stan")
fit_one <- mod$sample(data = data, chains = n_chains, init = init,
                      parallel_chains = n_chains,
                      iter_warmup = 500, iter_sampling = 500)

fit_one$time()
fit_one$summary()

yrep <- as.matrix(
  as_draws_df(
    fit_one$draws(variables = c("concentrationObsPred"))
  ))[, -(52:54)]

yobs <- data$cObs
time <- data$time[-1]

p <- bayesplot::ppc_ribbon(y = yobs, yrep = yrep, 
                           x = time, y_draw = "point") +
  xlab("time(h)") + ylab("drug plasma concentration (mg/L)") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18))
p


# estimate ELPD_loo for both models
log_lik_draws <- fit$draws("log_lik")
loo_estimate <- loo(log_lik_draws, r_eff = relative_eff(log_lik_draws))


log_lik_draws_one <- fit_one$draws("log_lik")
loo_estimate_one <-
  loo(log_lik_draws_one, r_eff = relative_eff(log_lik_draws_one))

print(loo_estimate_one)
print(loo_estimate)

#####################################################################
## Population model

data <- fromJSON(file = "data/twoCptPop.data.json")

# plot data
patientID <- rep(NA, data$nSubjects)
for (i in 1:data$nSubjects) {
  patientID[data$start[i]:data$end[i]] <- i 
}

p <- ggplot(data = data.frame(cObs = data$cObs, time = data$time[data$iObs],
                              ID = patientID[data$iObs]), aes(x = time, y = cObs)) +
  theme_bw() + geom_point() + facet_wrap(~ID)
p

# Draw initial points from the prior
init <- function () {
  n_subjects <- data$nSubjects
  pop_var <- c(0.2, 0.2, 0.2, 0.2, 0.2)
  
  CL_pop <- exp(rnorm(1, log(10), pop_var[1]))
  Q_pop <- exp(rnorm(1, log(15), pop_var[2]))
  VC_pop <- exp(rnorm(1, log(35), pop_var[3]))
  VP_pop <- exp(rnorm(1, log(105), pop_var[4]))
  ka_pop <- exp(rnorm(1, log(2.5), pop_var[5]))
  omega <- abs(rnorm(5, 0, pop_var))
  
  theta_pop <- c(CL_pop, Q_pop, VC_pop, VP_pop, ka_pop)
  theta <- matrix(NA, n_subjects, length(theta_pop))
  for (j in 1:n_subjects) {
    theta[j, ] <- exp(rnorm(length(theta_pop), log(theta_pop), omega))
  }
  
  list(CL_pop = CL_pop, Q_pop = Q_pop, VC_pop = VC_pop, VP_pop = VP_pop,
       ka_pop = ka_pop, omega = omega, theta = theta,
       sigma = abs(rnorm(1, 0, 1)))   
}


model_name <- "twoCptPop"
mod <- cmdstan_model(paste0("model/", model_name, ".stan"))

n_chains <- 4
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = n_chains,
                  iter_warmup = 500, iter_sampling = 500,
                  seed = 123, adapt_delta = 0.8)

fit$save_object(paste0("deliv/", model_name, ".fit.RDS"))

fit$time()

pars = c("lp__", "CL_pop", "Q_pop", "VC_pop", "VP_pop", "ka_pop", "sigma")
fit$summary(variables = pars)
bayesplot::mcmc_trace(fit$draws(), pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

# posterior predictive checks
yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))))
yrep <- yrep[, -((ncol(yrep) - 2):ncol(yrep))]

yobs <- data$cObs
time <- data$time[data$iObs]
patientID <- with(data, rep(1:nSubjects, each = nObs / nSubjects))

# within patient predictions
bayesplot::ppc_ribbon_grouped(y = yobs, yrep = yrep, x = time, patientID,
                              y_draw = "point")

# predictions for new patient
yrepNew <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("cObsNewPred"))))
yrepNew <- yrepNew[, -((ncol(yrepNew) - 2):ncol(yrepNew))]

bayesplot::ppc_ribbon_grouped(y = yobs, yrep = yrepNew, x = time, patientID,
                              y_draw = "point")


###############################################################################
## Fit model using within-chain parallelization

model_name <- "twoCptPop_rs"
mod <- cmdstan_model(paste0("model/", model_name, ".stan"),
                     cpp_options = list(stan_threads = TRUE))

n_chains <- 1
fit_rs <- mod$sample(data = data, chains = n_chains, init = init,
                     parallel_chains = n_chains,
                     threads_per_chain = 8,
                     iter_warmup = 500, iter_sampling = 500,
                     seed = 123, adapt_delta = 0.8)

fit_rs$time()
fit_rs$summary(variables = pars)
