# R code to Stan for Bayesian modelling and parameter inference

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
setwd('~/Dropbox/Projects/GBM/Codes/DIPG/DIPG7_ctrl')
set.seed(2000)

# CHANGE BASED ON MODEL
dataset <- 'StJude_ctrl'               # Choose dataset
model_id <- 'Gompertz'           # see "model keys" spreadsheet for reference
n_params <- 3                     # total number of local parameters in model (length of vector 'theta')
n_tumor_volume_params <- 2        # number of parameters in system of ODEs describing tumor volume dynamics between radiation fractions
n_conc_params <- 3                # number of parameters in system of ODEs describing the contribution a single administration of cct to 
# the concentration of the drug in the plasma
inc_cct <- FALSE                   # indicates whether a model incorporates the effect of cct
stan_file <- 'Gompertz_ctrl.stan' # listed in "model keys" spreadsheet
outputfile <-                     # pdf with diagnostics and predicted fits
  'Gompertz/ctrl_StJude_output.pdf'
rdsfile <-                        # file to save model fit
  'Gompertz/ctrl_StJude_modelfit.rds'


################################################################################################################
# convert data to long format, with one row per observation
convert_to_long <- function(data, control_cols) {
  control <- na.omit(gather(data[c(control_cols, "days")], mouse, volume, control_cols, factor_key=TRUE))
  long <- rbind(control)
  print('here')
  rownames(long) <- NULL
  long$days = long$days*24
  colnames(long)[1] <- "time"
  long$observation = 1
  n_mice <- length(unique(long$mouse))
  mouseids <- 1:n_mice
  names(mouseids) <- c(control_cols)
  long$mouse <- recode(long$mouse, !!!mouseids)
  long$time <- as.numeric(long$time)
  long$volume <- as.numeric(long$volume)
  long <-long[order(long$mouse, long$time),]
}

# Experimental data

# Vardon Experiment 1
if(dataset == 'Vardon_ctrl1') {
  control_cols_Vardon1 <- c('ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5')
  Vardon_ctrl1 <- read.csv('BLI_data_Vardon_Exp1.csv')
  data <- convert_to_long(Vardon_ctrl1, control_cols_Vardon1)
  n_control_mice <- 5
  first_obs_time <- 14

# Vardon Experiment 2 
} else if(dataset == 'Vardon_ctrl2') {
  control_cols_Vardon2 <- c('ctrl1', 'ctrl2', 'ctrl3')
  Vardon_ctrl2 <- read.csv('BLI_data_Vardon_Exp2.csv')
  data <- convert_to_long(Vardon_ctrl2, control_cols_Vardon2)
  n_control_mice <- 3
  first_obs_time <- 14
  
# St Jude control
} else if(dataset == 'StJude_ctrl') {
  control_cols_StJude <- c('ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5', 'ctrl6', 'ctrl7',
                           'ctrl8', 'ctrl9', 'ctrl10', 'ctrl11', 'ctrl12', 'ctrl13', 'ctrl14',
                           'ctrl15', 'ctrl16', 'ctrl17', 'ctrl18', 'ctrl19', 'ctrl20', 'ctrl21')
  StJude <- read.csv('StJude_DIPG7_ctrl.csv')
  data <- convert_to_long(StJude, control_cols_StJude)
  n_control_mice <- 21
  first_obs_time <- 0
}

n_mice = length(unique(data$mouse))
################################################################################################################

# Fit Stan model
# For the rt model (1 arm) we need the following variables:
stan_data <- list(n_times = nrow(data),
                  n_obs = nrow(data[data$observation == 1,]),
                  n_mice = length(unique(data$mouse)),
                  n_control_mice = n_control_mice,
                  start = c(1, which(diff(data[, 'time']) < 0) + 1),
                  end = c(c(1, which(diff(data[, 'time']) < 0))[-1],
                          nrow(data)),
                  t0 = data[data$time == first_obs_time, 'time'],
                  V_t0 = data[data$time == first_obs_time, 'volume'],
                  obs = which(data$observation == 1),
                  t_obs = data[data$observation == 1, 'time'],
                  V_obs = data[data$observation == 1, 'volume'],
                  times = data[, 'time'],
                  obs_per_mouse = matrix(sapply(1:n_mice, function(i) nrow(data[data$mouse == i & data$observation == 1,])), ncol = 1),
                  max_obs = max(sapply(1:n_mice, function(i) nrow(data[data$mouse == i & data$observation == 1,]))),
                  n_params = n_params)

# Model parameters to plot
if (model_id == 'Exponential') {
  params_plot <- c('loggamma_hat', 'logsigma', 'logV_t0_hat')
} else if (model_id == 'Logistic') {
  params_plot <- c('loggamma_hat', 'logK_hat', 'logsigma', 'logV_t0_hat')
} else if (model_id == 'Gompertz') {
  params_plot <- c('loggamma_hat', 'logbeta_hat' , 'logsigma', 'logV_t0_hat')
} 

# Additional parameters to monitor
additional_params <- c('V_pred', 'log_likelihood')

parameters <- c(params_plot, additional_params)

# Generate initial parameter estimates
# This is important as default initial parameter estimates by Stan are often poor
# Select random value from distribution so different chains have different initial parameter estimates
# This helps for ensuring convergence

if (model_id == 'Exponential') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.9, 1),
    logV_t0_hat = rnorm(1, 5.5, 0.5),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2.9, 5.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Logistic') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.9, 1),
    logK_hat = rnorm(1, 3, 1),
    logV_t0_hat = rnorm(1, 5.5, 0.5),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2.9, 3, 5.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Gompertz') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.9, 1),
    logbeta_hat = rnorm(1, 3, 1),
    logV_t0_hat = rnorm(1, 5.5, 0.5),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2.9, 3, 5.5), ea = n_mice), nrow = n_mice)
  )
} 

# Test/debug the model
test <- stan(stan_file,
             data = stan_data,
             chains = 2,
             iter = 15,
             control = list(adapt_delta = 0.999,
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
                           max_treedepth = 15),
            init = initialParameterValues)

################################################################################################################
# Diagnostics on Test
# mcmc_rhat(rhat(test, pars = params_plot)) + yaxis_text(hjust = 1)
# posterior <- as.array(test, pars = params_plot)
# mcmc_trace(posterior)


# Diagnostics
pdf(outputfile)
# That should be less than 1.05
color_scheme_set('viridis')
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
mcmc_pairs(posterior, np = nuts_params(fit), max_treedepth = 12)
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
  labs(x = 'Time (h)', y = 'Volume (p/s)') +
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
