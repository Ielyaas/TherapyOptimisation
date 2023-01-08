// Stan code for DIPG7 tumour growth and response to RT. Population-specific 
// rates of RT-dependent transfer to senescence and death.
// Simple model considers the dynamics of proliferating tumour cells (P),
// senescent cells (Q) and dying tumour cells (A).
// dP/dt = gamma*P-(delta1 + delta2)*P, dQ/dt = delta1*P, dA/dt = delta2*P.
// Tumour growth modelled as Exponential growth and population-specific
// decay of senescent cells based on DIPG7 data from 
// Ashley Vardon and data from St Jude hospital

// Ielyaas Cloete (13 November 2022)

functions {
   vector log_lik(vector shared_params, vector ind_params,
                data real[] logV_obs_padded, int[] x_int) {
      int n_params = x_int[2];
      vector[n_params] log_theta_hat = shared_params[1:n_params];
      matrix[n_params,n_params] Omega = to_matrix(shared_params[1+n_params:rows(shared_params) - 2], n_params, n_params);
      real sigma = shared_params[rows(shared_params) - 1];
      vector[n_params] logtheta = ind_params[1:n_params];
      int n = x_int[1];
      vector[n] V_hat_obs = ind_params[n_params+1:n_params+n];
      vector[n] logV_obs = to_vector(logV_obs_padded[1:n]);   
      real lp = multi_normal_lpdf(logtheta | log_theta_hat, Omega);
      real ll = normal_lpdf(logV_obs | log10(V_hat_obs), sigma);
      return [lp + ll]';
   }

  // analytically solved ODEs
  real[,] cell_population_dynamics_num(real[] times,	  // Times at which system dynamics are evaluated
	                                real[] N0,		        // initial conditions
                                  real[] params)   {    // Parameters for system

    real gamma = params[1];                             // Tumour growth rate
    real delta = params[2];                            // Apoptosis rate

    real N[size(times),3];

    for (t in 1:size(times)) {
      
      N[t,1] = N0[1]*exp((gamma-delta)*(times[t] - times[1])); // proliferating compartment dynamics
      N[t,2] = N0[1]/(gamma)*(exp(gamma*(times[t] - times[1])) - 1) + N0[2]; // senescent compartment dynamics
      N[t,3] = N0[1]/(gamma - delta)*(exp((gamma - delta)*(times[t] - times[1])) - 1) + N0[3]; // dying compartment dynamics
    }

    return N;
  }
  
  // Computes the tumor volume dynamics for untreated mice
  real[] control_tumor_volume_dynamics(real[] times, real V_t0,
                                       real gamma, real delta,
                                       real[] x_r, int[] x_i) {
    
    real N_t0[3];                                         // Initial state
    real params[2];                                       // ODE parameters
    matrix[size(times),3] N_hat;                         // Predicted cell population counts
    real V_hat[size(times)];                              // Predicted tumor volume
    
    params[1] = gamma;
    params[2] = delta;

    // initial cell counts
    N_t0[1] = V_t0;
    N_t0[2] = 0;
    N_t0[3] = 0;

    N_hat = to_matrix(cell_population_dynamics_num(times, N_t0, params));
    
    // Sum over cell types to get tumor volume
    for (t in 1:size(times)) {
      V_hat[t] = sum(N_hat[t]);
    }
    return V_hat;
  }
  
  // computes the tumor volume dynamics for mice treated with radiation
  real[] tumor_volume_dynamics(real[] times, real[] rt_times, real[] rt_doses, int[] fraction_index, real V_t0,
                               real gamma, real delta, real alpha1, real alpha2,
                               real[] x_r, int[] x_i) {
    
    real N_t0[3];                                         // Initial state
    real params[2];                                       // ODE parameters
    matrix[size(times), 3] N_hat;                         // Predicted cell population numbers
    real V_hat[size(times)];                              // Predicted tumor volume
    
    params[1] = gamma;
    params[2] = delta;

    // Cell population dynamics before the first fraction of radiation
    N_t0[1] = V_t0;
    N_t0[2] = 0;
    N_t0[3] = 0;

    N_hat[1:fraction_index[1],] = to_matrix(cell_population_dynamics_num(times[1:fraction_index[1]], N_t0, params));
    
    // Cell population dynamics between fractions of radiation
    for (fraction in 1:size(rt_times) - 1) {
      // 'instantaneous' effect of each rt fraction on cell population counts
      N_t0[1] = N_hat[fraction_index[fraction], 1]*exp(-(alpha1 + alpha2)*rt_doses[fraction]);
      N_t0[2] = N_hat[fraction_index[fraction], 2] +
        N_hat[fraction_index[fraction], 1]*(1 - exp(-alpha1*rt_doses[fraction]));
      N_t0[3] = N_hat[fraction_index[fraction], 3] +
        N_hat[fraction_index[fraction], 1]*(1 - exp(-alpha2*rt_doses[fraction]));  
      
      N_hat[fraction_index[fraction] + 1:fraction_index[fraction + 1],] = 
        to_matrix(cell_population_dynamics_num(times[fraction_index[fraction] + 1:fraction_index[fraction + 1]], N_t0, params));
    } 
    
    // Cell population dynamics from final fraction of radiation to final tumor volume measurement
    N_t0[1] = N_hat[fraction_index[size(fraction_index)], 1]*exp(-(alpha1 + alpha2)*rt_doses[size(rt_times)]);
    N_t0[2] = N_hat[fraction_index[size(fraction_index)], 2] +
      N_hat[fraction_index[size(fraction_index)], 1]*exp(-alpha1*rt_doses[size(fraction_index)]);
    N_t0[3] = N_hat[fraction_index[size(fraction_index)], 3] +
      N_hat[fraction_index[size(fraction_index)], 1]*exp(-alpha2*rt_doses[size(fraction_index)]);  
    
    N_hat[fraction_index[size(fraction_index)] + 1:size(times), ] = 
      to_matrix(cell_population_dynamics_num(times[fraction_index[size(fraction_index)] + 1:size(times)], N_t0, params));

    // Sum over cell types to get tumor volume
    for (t in 1:size(times)) {
      V_hat[t] = sum(N_hat[t]);
    }
    return V_hat;
  }
}

data {
  
  int<lower = 1> n_times;                                   // Number of times (observations and fractions)
  int<lower = 1> n_obs;                                     // Number of observations
  int<lower = 1> n_fractions;                               // Number of radiation fractions
  int<lower = 1> n_mice;                                    // Number of mice
  int<lower = 0> n_control_mice;                            // Number of control mice
  int<lower = 0> n_untreated_mice;                          // Number of untreated mice
  real<lower = 0> rt_times[n_fractions];                    // Times of radiation administrations
  real<lower = 0> rt_doses[n_fractions];                    // Radiation doses
  int<lower = 0> schedule_start[n_mice - n_untreated_mice]; // index of first radiation fraction for each mouse
  int<lower = 0> schedule_end[n_mice - n_untreated_mice];   // index of last radiation fraction for each mouse
  int<lower = 0> start[n_mice];                             // index of first observation for each mouse
  int<lower = 0> end[n_mice];                               // index of last observation for each mouse
  int<lower = 0> obs[n_obs];                                // Observation row numbers
  real<lower = 0> t_obs[n_obs];                             // Measurement times
  real<lower = 0> V_obs[n_obs];                             // Measured tumor volumes
  real<lower = 0> times[n_times];                           // Times (observations and fractions)
  int<lower = 0> obs_per_mouse[n_mice, 1];                  // number of observations for each fmouse
  int<lower = 0> max_obs;                                   // max number of observations per mouse
  int<lower = 0> n_params;                                  // number of individual parameters for each mouse
  int<lower = 0> fraction_index[n_fractions];               // index of each radiation fraction           
}

transformed data {
  
  real x_r[0];
  int x_i[0];
  real logV_obs[n_obs];                                   // Log(measured tumor volumes)
  real logV_obs_padded[n_mice, max_obs];                  // padded to be passed into map function

  logV_obs = log10(V_obs);
  
  // add 0's so each row has same number of columns (observations). This rectangular format is necessary for the mapped implementation.
  for (i in 1:n_mice) {
    int pos = 1;
    for (j in 1:n_mice) {
      int n = obs_per_mouse[j][1];
      int stop = pos + n - 1;
      logV_obs_padded[j, 1:n] = to_array_1d(logV_obs[pos:stop]);
      logV_obs_padded[j, (n+1):max_obs] = to_array_1d(rep_array(0, max_obs - n));
      pos += n;
    }
  }
}

parameters {
  
  // population-level parameters on log scale
  real loggamma_hat;                              // Growth rate
  real logdelta_hat;                             // rate of continous transfer from proliferating to dying compartment
  real logalpha1_hat;                             // instantaneous transfer rate at radiation times
  real logalpha2_hat;
  real logV_t0_hat;                               // initial tumor volume
  
  vector[n_params] logtheta[n_mice];              // Individual level parameter values on log scale

  corr_matrix[n_params] rho;                      // Correlation matrix describes correlations between parameters
  vector<lower = 0>[n_params] omega;              // Standard deviations of parameters
  real logsigma;                                  // Residual standard deviation on log scale
  
}

transformed parameters {
  
  vector[n_params] log_theta_hat;                 // Population-level parameters on log scale, as single vector. Define means of distributions of local parameters
  cov_matrix[n_params] Omega;                     // Covariance matrix
  real<lower = 0> sigma;                          // residual standard deviation
  
  // individual-level parameters
  real<lower = 0> gamma[n_mice];                          
  real<lower = 0> delta[n_mice]; 
  real<lower = 0> alpha1[n_mice];
  real<lower = 0> alpha2[n_mice];
  real<lower = 0> V_t0[n_mice]; 

  real<lower = 0> V_hat[n_times];                         // Estimated tumor volume at measurement and treatment times
  real<lower = 0> V_hat_obs[n_obs];                       // Estimated tumor volume at volume measurement times
  vector<lower = 0>[max_obs] V_hat_obs_padded[n_mice];    // Estimated tumor volumes at measurement times padded to be passed to map function

  log_theta_hat[1] = loggamma_hat;                               // Populataion level mean growth rate
  log_theta_hat[2] = logdelta_hat;
  log_theta_hat[3] = logalpha1_hat;
  log_theta_hat[4] = logalpha2_hat;
  log_theta_hat[5] = logV_t0_hat;
  
  sigma = pow(10,logsigma);
  
  Omega = quad_form_diag(rho, omega);                     // diag_matrix(omega) * rho * diag_matrix(omega)
  
  for (i in 1:n_mice) {
    gamma[i] = pow(10,logtheta[i, 1]);
    delta[i] = pow(10, logtheta[i, 2]);
    alpha1[i] = pow(10, logtheta[i, 3]);
    alpha2[i] = pow(10, logtheta[i, 4]);
    V_t0[i] = pow(10, logtheta[i, 5]);
    
    // estimate volumes for untreated mice
    if (i <= n_untreated_mice)
      V_hat[start[i]:end[i]] =
        control_tumor_volume_dynamics(times[start[i]:end[i]], V_t0[i],
                                      gamma[i], delta[i], x_r, x_i);
                                      
    // estimate volumes for mice treated with radiation
    else
      V_hat[start[i]:end[i]] =
        tumor_volume_dynamics(times[start[i]:end[i]],
                              rt_times[schedule_start[i - n_untreated_mice]:schedule_end[i - n_untreated_mice]],
                              rt_doses[schedule_start[i - n_untreated_mice]:schedule_end[i - n_untreated_mice]],
                              fraction_index[schedule_start[i - n_untreated_mice]:schedule_end[i - n_untreated_mice]],
                               V_t0[i],
                              gamma[i], delta[i], alpha1[i], alpha2[i], x_r, x_i);
    
  }
  
  // add 0's so each element has same length vector of observations
  V_hat_obs = V_hat[obs];
  for (i in 1:n_mice) {
    int pos = 1;
    for (j in 1:n_mice) {
      int n = obs_per_mouse[j][1];
      int stop = pos + n - 1;
      V_hat_obs_padded[j, 1:n] = to_vector(V_hat_obs[pos:stop]);
      V_hat_obs_padded[j, (n+1):max_obs] = to_vector(rep_array(0, max_obs - n));
      pos += n;
    }
  }
}

model {
  vector[n_params + max_obs] ind_params[n_mice];
  int x_int[n_mice,3];
  int n_params_rep[n_mice, 1] = rep_array(n_params, n_mice, 1); 

  // Priors
  loggamma_hat ~ normal(-2, 1);
  logdelta_hat ~  normal(-2, 1);
  logalpha1_hat ~ normal(-2, 1);
  logalpha2_hat ~ normal(-2, 1);
  logV_t0_hat ~ normal(1.5,1);

  omega ~ cauchy(0, 1);                                   // Half-cauchy prior (see Bayesian Data Analysis, Gelman et al. p.130)
  rho ~ lkj_corr(1);                                      // Approximately equivalent to uniform distribution on correlations

  logsigma ~ normal(-1.5, 2);

  // format parameters and data to be passed to map function
  for (r in 1:n_mice) {
    ind_params[r] = to_vector(append_row(logtheta[r], V_hat_obs_padded[r]));
    x_int[r] = to_array_1d(append_array(obs_per_mouse[r], n_params_rep[r]));
  }
  
  // compute log likelihood with given parameters
  target += sum(map_rect(log_lik, 
    append_row(append_row(log_theta_hat, append_row(to_vector(Omega), sigma)), n_params), 
    ind_params, 
    logV_obs_padded, 
    x_int
    ));
}

generated quantities {
  
  real<lower = 0> V_pred[n_obs];
  vector[n_obs] log_likelihood;
  
  // Random number generator drawing from normal distribution
  for (ob in 1:n_obs) {
    V_pred[ob] = pow(10, normal_rng(log10(V_hat_obs[ob]), sigma));
    log_likelihood[ob] = normal_rng(log10(V_hat_obs[ob]), sigma);
  }
  
}
