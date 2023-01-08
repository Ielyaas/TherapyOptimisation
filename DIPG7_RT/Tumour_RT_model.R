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
dataset <- 'L7S11'           # 'Vardon_RT1', 'Vardon_RT2', 'L7S1', 'L7S2', 'L7S3'
model_id <- 'Logistic4'           # 'Exponential', 'Exponential1', 'Exponential2', 'Exponential3', 'Exponential4' , 
                                  # 'Logistic', 'Logistic2', 'Logistic3', 'Logistic4',
                                  # 'Gompertz1', 'Gompertz2', 'Gompertz3', 'Gompertz4'
n_params <- 6                     # total number of local parameters in model (length of vector 'theta')
n_tumor_volume_params <- 3        # number of parameters in system of ODEs describing tumor volume dynamics between radiation fractions
stan_file <- 'Logistic_RT4.stan'     # listed in "model keys" spreadsheet
outputfile <-                     # pdf with diagnostics and predicted fits
  'Logistic/RT4_L7S11_output.pdf'
rdsfile <-                        # file to save model fit
  'Logistic/RT4_L7S11_modelfit.pdf'

################################################################################################################
# convert data to long format, with one row per observation
convert_to_long <- function(data, control_cols, rt_cols, rt_times, rt_doses) {
  control <- na.omit(gather(data[c(control_cols, "days")], mouse, volume, control_cols, factor_key=TRUE))
  rt <- na.omit(gather(data[c(rt_cols, "days")], mouse, volume, rt_cols, factor_key=TRUE))
  long <- rbind(control, rt)
  print('here')
  rownames(long) <- NULL
  long$days = long$days*24
  colnames(long)[1] <- "time"
  long$observation = 1
  long$dose = NA
  long$rt_fraction = 0
  for (rtmouse in rt_cols) {
    for (i in 1:length(rt_times)) {
      long = rbind(long, c(rt_times[i], rtmouse, NA, 0, rt_doses[i], 1, 0))
    }
  }
  n_mice <- length(unique(long$mouse))
  mouseids <- 1:n_mice
  names(mouseids) <- c(control_cols, rt_cols)
  long$mouse <- recode(long$mouse, !!!mouseids)
  long$time <- as.numeric(long$time)
  long$dose <- as.numeric(long$dose)
  long$volume <- as.numeric(long$volume)
  long <-long[order(long$mouse, long$time),]
}

# Experimental data

# Vardon Experiment 1
if(dataset == 'Vardon_RT1') {
  control_cols_Vardon1 <- c('ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5')
  rt_cols_Vardon1 <- c("rt1", "rt2", "rt3", "rt4", "rt5")
  rt_times_6 <- c(840, 888, 936, 1008, 1056, 1104)
  rt_doses_2Gyx6 <- c(2, 2, 2, 2, 2, 2)
  Vardon_RT1 <- read.csv('BLI_data_Vardon_Exp1.csv')
  Vardon_RT1
  data <- convert_to_long(Vardon_RT1, control_cols_Vardon1, rt_cols_Vardon1, rt_times_6, rt_doses_2Gyx6)
  n_control_mice <- 5
  n_rt_mice <- 5
  n_untreated_mice <- n_control_mice
  first_rt_time <- min(rt_times_6)

# Vardon Experiment 2  
} else if(dataset == 'Vardon_RT2') {
  control_cols_Vardon2 <- c('ctrl1', 'ctrl2', 'ctrl3')
  rt_cols_Vardon2 <- c("rt_8Gy1", "rt_8Gy2", "rt_8Gy3", "rt_12Gy1", "rt_12Gy2", "rt_12Gy3")
  rt_times_4 <- c(840, 936, 1008, 1104)
  rt_doses_2Gyx4 <- c(2, 2, 2, 2)
  rt_times_6 <- c(840, 888, 936, 1008, 1056, 1104)
  rt_doses_2Gyx6 <- c(2, 2, 2, 2, 2, 2)
  Vardon_RT2 <- read.csv('BLI_data_Vardon_Exp2.csv')
  Vardon_RT2
  data <- convert_to_long2(Vardon_RT2, control_cols_Vardon2, rt_cols_Vardon2, 
                          rt_times_4, rt_doses_2Gyx4, rt_times_6, rt_doses_2Gyx6)
  n_control_mice <- 3
  n_rt_mice <- 6
  n_untreated_mice <- n_control_mice
  first_rt_time <- min(rt_times_4)

# St Jude L7S11
} else if(dataset == 'L7S11') {
  control_cols_L7S11 <- c('ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5', 'ctrl6')
  rt_cols_L7S11 <- c("rt1", "rt2", "rt3", "rt4", "rt5", "rt6")
  rt_times_20Gy <- c(648, 672, 696, 720, 744, 816, 840, 864, 888, 912)
  rt_doses_20Gy <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
  L7S11 <- read.csv('StJude_L7S11.csv')
  L7S11
  data <- convert_to_long(L7S11, control_cols_L7S11, rt_cols_L7S11, rt_times_20Gy,
                           rt_doses_20Gy)
  n_control_mice <- 6
  n_rt_mice <- 6
  n_untreated_mice <- n_control_mice
  first_rt_time <- min(rt_times_20Gy)

# St Jude L7S12
} else if(dataset == 'L7S12') {
  control_cols_L7S12 <- c('ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5', 'ctrl6')
  rt_cols_L7S12 <- c("rt1", "rt2", "rt3", "rt4", "rt5", "rt6")
  rt_times_54Gy <- c(648, 672, 696, 720, 744, 816, 840, 864, 888, 912,
                     984, 1008, 1032, 1056, 1080, 1152, 1176, 1200, 1224, 1248,
                     1320, 1344, 1368, 1392, 1416, 1488, 1512)
  rt_doses_54Gy <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                     2, 2, 2, 2, 2, 2, 2)
  L7S12 <- read.csv('StJude_L7S12.csv')
  L7S12
  data <- convert_to_long(L7S12, control_cols_L7S12, rt_cols_L7S12, rt_times_54Gy, rt_doses_54Gy)
  n_control_mice <- 6
  n_rt_mice <- 6
  n_untreated_mice <- n_control_mice
  first_rt_time <- min(rt_times_54Gy)

# St Jude L7S2
} else if(dataset == 'L7S2') {
  control_cols_L7S2 <- c('ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5', 'ctrl6', 'ctrl7')
  rt_cols_L7S2 <- c("rt1", "rt2", "rt3", "rt4", "rt5", "rt6", "rt7")
  rt_times_20Gy <- c(648, 672, 696, 720, 744, 816, 840, 864, 888, 912,
                     984, 1008, 1032, 1056, 1080, 1152, 1176, 1200, 1224, 1248,
                     1320, 1344, 1368, 1392, 1416, 1488, 1512)
  rt_doses_20Gy <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
  L7S2 <- read.csv('StJude_L7S2.csv')
  L7S2
  data <- convert_to_long(L7S2, control_cols_L7S2, rt_cols_L7S2, rt_times_20Gy, rt_doses_20Gy)
  n_control_mice <- 7
  n_rt_mice <- 7
  n_untreated_mice <- n_control_mice
  first_rt_time <- min(rt_times_20Gy)
}

n_mice = length(unique(data$mouse))
################################################################################################################

# Fit Stan model
stan_data <- list(n_times = nrow(data),
                  n_obs = nrow(data[data$observation == 1,]),
                  n_fractions = nrow(data[data$rt_fraction == 1,]),
                  n_rt_fractions = nrow(data[data$rt_fraction == 1,]),
                  n_mice = length(unique(data$mouse)),
                  n_control_mice = n_control_mice,
                  n_untreated_mice = n_untreated_mice,
                  n_rt_mice = n_rt_mice,
                  rt_times = data[data$rt_fraction == 1, 'time'],
                  rt_doses = data[data$rt_fraction == 1, 'dose'],
                  schedule_start = which(data[data$rt_fraction == 1,]$time == first_rt_time),
                  schedule_end = c((which(data[data$rt_fraction == 1,]$time == first_rt_time) - 1)[-1],
                                   length(data[data$rt_fraction == 1, 'time'])),
                  start = c(1, which(diff(data[, 'time']) < 0) + 1),
                  end = c(c(1, which(diff(data[, 'time']) < 0))[-1],
                          nrow(data)),
                  obs = which(data$observation == 1),
                  t_obs = data[data$observation == 1, 'time'],
                  V_obs = data[data$observation == 1, 'volume'],
                  times = data[, 'time'],
                  obs_per_mouse = matrix(sapply(1:n_mice, function(i) nrow(data[data$mouse == i & data$observation == 1,])), ncol = 1),
                  max_obs = max(sapply(1:n_mice, function(i) nrow(data[data$mouse == i & data$observation == 1,]))),
                  n_params = n_params, 
                  n_tumor_volume_params = n_tumor_volume_params,
                  fraction_index = as.vector(sapply((n_control_mice +1):n_mice, function(i) which(data[data$mouse == i,]$rt_fraction == 1))))

# Model parameters to plot
if (model_id == 'Exponential') {
  params_plot <- c('loggamma_hat', 'logalpha_hat', 'logsigma', 'logV_t0_hat')
} else if (model_id == 'Exponential1' || model_id == 'Exponential2' || model_id == 'Exponential3') {
  params_plot <- c('loggamma_hat', 'logdelta_hat', 'logalpha_hat', 'logsigma', 'logV_t0_hat')
} else if (model_id == 'Exponential4') {
  params_plot <- c('loggamma_hat', 'logdelta_hat', 'logalpha1_hat', 'logalpha2_hat', 'logsigma', 'logV_t0_hat')
} else if (model_id == 'Logistic' || model_id == 'Logistic2' || model_id == 'Gompertz1' || model_id == 'Gompertz2') {
  params_plot <- c('loggamma_hat', 'logK_hat', 'logalpha_hat', 'logsigma', 'logV_t0_hat')
} else if (model_id == 'Logistic3' || model_id == 'Gompertz3') {
  params_plot <- c('loggamma_hat', 'logK_hat', 'logdelta_hat', 'logalpha_hat', 'logsigma', 'logV_t0_hat')
} else if (model_id == 'Logistic4' || model_id == 'Gompertz4') {
  params_plot <- c('loggamma_hat', 'logK_hat', 'logdelta_hat', 'logalpha1_hat', 'logalpha2_hat', 'logsigma', 'logV_t0_hat')
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
    loggamma_hat = rnorm(1, -2, 1),
    logalpha_hat = rnorm(1, -2, 0.1),
    logV_t0_hat = rnorm(1, 215.5, 0.1),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 1),
    logtheta = matrix(rep(c(-2, -2, 215.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Exponential1' || model_id == 'Exponential2' || model_id == 'Exponential3') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2, 1),
    logdelta_hat = rnorm(1, -1, 1),
    logalpha_hat = rnorm(1, -2, 1),
    logV_t0_hat = rnorm(1, 5.5, 0.5),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2, -1, -2, 5.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Exponential4') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2, 1),
    logdelta_hat = rnorm(1, -1, 1),
    logalpha1_hat = rnorm(1, -2, 1),
    logalpha2_hat = rnorm(1, -2, 1),
    logV_t0_hat = rnorm(1, 1.5, 1),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2, -1, -2, -2, 1.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Logistic' || model_id == 'Logistic2') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.5, 1),
    logK_hat = rnorm(1, 2.75, 1),
    logalpha_hat = rnorm(1, -2, 1),
    logV_t0_hat = rnorm(1, 1.5, 1),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2.5, 2.75, -2, 1.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Logistic3') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.5, 1),
    logK_hat = rnorm(1, 2.75, 1),
    logdelta_hat = rnorm(1, -1, 1),
    logalpha_hat = rnorm(1, -2, 1),
    logV_t0_hat = rnorm(1, 1.5, 1),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2.5, 2.75, -1, -2, 1.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Logistic4') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.5, 1),
    logK_hat = rnorm(1, 2.75, 1),
    logdelta_hat = rnorm(1, -1, 1),
    logalpha1_hat = rnorm(1, -2, 1),
    logalpha2_hat = rnorm(1, -2, 1),
    logV_t0_hat = rnorm(1, 5.5, 0.1),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 1),
    logtheta = matrix(rep(c(-2.5, 2.75, -1, -2, -2, 5.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Gompertz1' || model_id == 'Gompertz2') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.5, 1),
    logK_hat = rnorm(1, -9, 1),
    logalpha_hat = rnorm(1, -2, 1),
    logV_t0_hat = rnorm(1, 1.5, 1),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2.5, -9, -2, 1.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Gompertz3') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.5, 1),
    logK_hat = rnorm(1, -9, 1),
    logdelta_hat = rnorm(1, -2, 1),
    logalpha_hat = rnorm(1, -2, 1),
    logV_t0_hat = rnorm(1, 1.5, 1),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2.5, -9, -2, -2, 1.5), ea = n_mice), nrow = n_mice)
  )
} else if (model_id == 'Gompertz4') {
  initialParameterValues = function() list(
    loggamma_hat = rnorm(1, -2.5, 1),
    logK_hat = rnorm(1, -9, 1),
    logdelta_hat = rnorm(1, -2, 1),
    logalpha1_hat = rnorm(1, -2, 1),
    logalpha2_hat = rnorm(1, -2, 1),
    logV_t0_hat = rnorm(1, 1.5, 1),
    omega = exp(rnorm(n_params, log(0.25), 0.5)),
    rho = diag(n_params),
    logsigma = rnorm(1, -1.5, 2),
    logtheta = matrix(rep(c(-2.5, -9, -2, -2, -2, 1.5), ea = n_mice), nrow = n_mice)
  )
} 

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
            chains = 3,
            warmup = 400,
            iter = 1200,
            control = list(adapt_delta = 0.8, #0.8
                           stepsize = 0.1,
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
