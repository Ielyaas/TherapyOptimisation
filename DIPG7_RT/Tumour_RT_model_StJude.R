# Error using macOS Catalina 10.15 --> installl c++ toolchain from this website 
# https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-Source#mac

Sys.setenv('STAN_NUM_THREADS' = 4)

library(Rcpp)
library(rstan)
library(bayesplot)
library(deSolve)
library(reshape2)
library(dplyr)
library(tidyr)
library(loo)


# Run different chains in parallel on different cores
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

theme_set(theme_bw())
setwd('~/Dropbox/Projects/GBM/Codes/DIPG/DIPG7_RT')
set.seed(2000)

# CHANGE BASED ON MODEL
stan_file <- 'Logistic_RT5.stan'     # listed in "model keys" spreadsheet
outputfile <-                     # pdf with diagnostics and predicted fits
  'Logistic/RT5_StJudeL7_output.pdf'
rdsfile <-                        # file to save model fit
  'Logistic/RT5_StJudeL7_modelfit.pdf'

################################################################################################################

# Experimental data

data <- read.csv('StJude_L7.csv')

n_mice = length(unique(data$mouse))
n_control_mice <- 21
first_obs_time <- 0

# Determine number of time points before the first and following the last treatment administration
first_tx_times <- array(0, dim = n_mice)
final_tx_times <- array(0, dim = n_mice)
n_pre_tx_times <- array(0, dim = n_mice)
n_post_tx_times <- array(0, dim = n_mice)

for (i in unique(data[data$observation != 1, 'mouse'])) {
  times <- data[data$mouse == i, 'time']
  tx_times <- data[data$mouse == i & data$observation != 1, 'time']
  first_tx_index <- which(times == tx_times[1])
  final_tx_index <- which(times == tx_times[length(tx_times)])
  first_tx_times[i] <- tx_times[1]
  final_tx_times[i] <- tx_times[length(tx_times)]
  n_pre_tx_times[i] <- first_tx_index - 1
  n_post_tx_times[i] <- length(times) - final_tx_index
}

schedule_start <- array(0, dim = n_mice)
schedule_end <- array(0, dim = n_mice)

for (i in unique(data[data$observation != 1, 'mouse'])) {
  schedule_start[i] <-
    which(data[data$observation != 1 & data$mouse == i, 'time'] == first_tx_times[i]) +
    nrow(data[data$observation != 1 & data$mouse < i,])
  schedule_end[i] <-
    c((which(data[data$observation != 1 & data$mouse == i, 'time'] == first_tx_times[i]) - 1)[-1],
      length(data[data$observation != 1 & data$mouse == i, 'time'])) +
    nrow(data[data$observation != 1 & data$mouse < i,])
}

################################################################################################################

# Fit Stan model

# For the stan model we need the following variables:
stan_data <- list(n_times = nrow(data[data$time > first_obs_time,]),
                  n_obs = nrow(data[data$observation == 1 & data$time > first_obs_time,]),
                  n_admin = length(data[data$observation != 1, 'time']),
                  n_mice = length(unique(data$mouse)),
                  n_control_mice = n_control_mice,
                  tx_times = data[data$observation != 1, 'time'],
                  tx_doses = data[data$observation != 1, 'dose'],
                  schedule_start = schedule_start,
                  schedule_end = schedule_end,
                  start = c(1, which(diff(data[data$time > first_obs_time, 'time']) < 0) + 1),
                  end = c(c(1, which(diff(data[data$time > first_obs_time, 'time']) < 0))[-1],
                          nrow(data[data$time > first_obs_time,])),
                  n_pre_tx_times = n_pre_tx_times,
                  n_post_tx_times = n_post_tx_times,
                  t0 = data[data$time == first_obs_time, 'time'],
                  V_t0 = data[data$time == first_obs_time, 'volume'],
                  obs = which(data[data$time > first_obs_time,]$observation == 1),
                  t_obs = data[data$observation == 1 & data$time > first_obs_time, 'time'],
                  V_obs = data[data$observation == 1 & data$time > first_obs_time, 'volume'],
                  times = data[data$time > first_obs_time, 'time'],
                  run_estimation = 1)

# Model parameters to plot
params_plot <- c('loggamma_hat', 'logK_hat', 'logdelta_hat', 'logalpha1_hat',
                 'logalpha2_hat', 'logsigma')


# Additional parameters to monitor
additional_params <- c('V_pred', 'log_likelihood')

parameters <- c(params_plot, additional_params)

# Generate initial parameter estimates
# This is important as default initial parameter estimates by Stan are often poor
# Select random value from distribution so different chains have different initial parameter estimates
# This helps for ensuring convergence

initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.5, 1),
    logK_hat = rnorm(1, 2.75, 1),
    logdelta_hat = rnorm(1, -1, 1),
    logalpha1_hat = rnorm(1, -2, 1),
    logalpha2_hat = rnorm(1, -2, 1),
    omega = exp(rnorm(5, log(0.25), 0.5)),
    rho = diag(5),
    logsigma = rnorm(1, -1.5, 1),
    logtheta = matrix(rep(c(-2.5, 2.75, -1, -2, -2), ea = n_mice), nrow = n_mice)
)


# Test/debug the model
test <- stan(stan_file,
             data = stan_data,
             chains = 1,
             iter = 10,
             control = list(adapt_delta = 0.8,
                            stepsize = 0.1,
                            max_treedepth = 10),
             init = initialParameterValues)

# Fit and sample from the posterior
fit <- stan(fit = test,
            data = stan_data,
            pars = parameters,
            cores = 4,
            chains = 4,
            warmup = 500,
            iter = 1500,
            control = list(adapt_delta = 0.999, #0.8
                           stepsize = 0.05,
                           max_treedepth = 20),
            init = initialParameterValues)

################################################################################################################
# Diagnostics on Test
# mcmc_rhat(rhat(test, pars = params_plot)) + yaxis_text(hjust = 1)
# posterior <- as.array(test, pars = params_plot)
# mcmc_trace(posterior)


# Diagnostics
# pdf(outputfile)
# That should be less than 1.05
color_scheme_set('brightblue')
mcmc_rhat(rhat(fit, pars = params_plot)) + yaxis_text(hjust = 1)

# MCMC diagnostics
posterior <- as.array(fit, pars = params_plot)
# Chains should be consistent
mcmc_trace(posterior)
# Red points - divergent transitions after warmup
# Yellow points - maximum treedepth exceeded
# Divergent transitions above the diagonal of the pairs() plot (meaning that the amount of
# numerical error was above the median over the iterations) can often be eliminated simply by
# increasing the value of the adapt_delta parameter
# see if can label divergent chains - or just use pairs
mcmc_pairs(posterior, np = nuts_params(fit), max_treedepth = 10)
# 
# ################################################################################################################

# Posterior predictive distributions
mcmc_dens(posterior)

# Prediction of future observations in a new study, i.e. posterior predictive distributions
observations <- data.frame(volume = data[data$observation == 1, 'volume'],
                           mouse = data[data$observation == 1, 'mouse'],
                           time = data[data$observation == 1, 'time'])

predictions <- as.data.frame(fit, pars = 'V_pred') %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lowerBound = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            upperBound = quantile(value, probs = 0.95)) %>%
  bind_cols(observations)
p1 <- ggplot(predictions, aes(x = time, y = volume))
p1 <- p1 + geom_point() + scale_y_log10() +
  labs(x = 'Time (h)', y = 'Volume (log) (p/s)') +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = 'none', strip.text = element_text(size = 8))
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lowerBound, ymax = upperBound), alpha = 0.25) +
  facet_wrap(~factor(mouse))
dev.off()

saveRDS(fit, rdsfile) # use readRDS(filename) to re-load stored model fit


# ################################################################################################################
# 
# # Evaluate model performance
# http://mc-stan.org/loo/articles/loo2-with-rstan.html
# # Extract pointwise log-likelihood and compute LOO
log_lik <- extract_log_lik(fit, 'log_likelihood', merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik))
loo <- loo(log_lik, r_eff = r_eff, cores = 2)
print(loo)
