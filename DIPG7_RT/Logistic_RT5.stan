// Stan code for DIPG7 tumour growth and response to RT. 
// Simple model considers the dynamics of proliferating tumour cells (P),
// senescent cells (Q) and dying tumour cells (A).
// dP/dt = gamma*P(1 - P/K) - delta*P, dQ/dt = 0, 
// dA/dt = delta*P.
// Tumour growth modelled as Exponential growth and decay of senescent cells
// based on DIPG7 data from 
// Ashley Vardon and data from St Jude hospital

// Ielyaas Cloete (13 December 2022)

functions {
   

  // analytically solved ODEs
  real[,] cell_population_dynamics_num(real[] times,	  // Times at which system dynamics are evaluated
	                                real[] N0,		        // initial conditions
                                  real[] params)   {    // Parameters for system

    real gamma = params[1];                               // Growth rate
    real k = params[2];                                   // carrying capacity
    real delta = params[3];                              // Death rate

    real N[size(times),3];                      // estimated cell counts

    for (t in 1:size(times)) {
      
      N[t,1] = k*N0[1]*exp(gamma *(times[t] - times[1]))/(k + N0[1]*((exp(gamma*(times[t] - times[1]))) - 1));
      N[t,2] = N0[2]*((k - N0[1])/k)*(1 - exp(-gamma*(times[t] - times[1]))*(N0[1]/(N0[1] - k)))*exp(-gamma*(times[t] - times[1]));
      N[t,3] = delta*N0[3]*k*N0[1]*exp(gamma *(times[t] - times[1]))/(k + N0[1]*((exp(gamma*(times[t] - times[1]))) - 1));
  }
    return N;
  }
  
  // Computes the tumor volume dynamics for untreated mice
  real[] control_tumor_volume_dynamics(real[] times, real V_t0,
                                       real gamma, real k, real delta,
                                       real[] x_r, int[] x_i) {
    
    real N_t0[3];                                         // Initial state
    real params[3];                                       // ODE parameters
    matrix[size(times), 3] N_hat;                         // Predicted cell population numbers
    real V_hat[size(times)];                              // Predicted tumor volume
    
    params[1] = gamma;
    params[2] = k;
    params[3] = delta;

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
  
  // Computes the tumor volume dynamics following radiation given the cell population numbers
  // computed by the cell_population_dynamics() function
  real[] tumor_volume_dynamics(real[] times, real[] tx_times, real[] tx_doses,
                               int n_pre_tx_times, int n_post_tx_times, real t0,  
                               real V_t0, real gamma, real k, real delta, real alpha1, real alpha2,
                               real[] x_r, int[] x_i) {
    
    real N_t0[3];                                         // Initial state
    real params[3];                                       // ODE parameters
    int tx_index[size(tx_times)];
    matrix[size(times), 3] N_hat;                         // Predicted cell population numbers
    real V_hat[size(times)];                              // Predicted tumor volume
    
    params[1] = gamma;
    params[2] = k;
    params[3] = delta;
    
    // Determine indices of treatment administrations in a times array
    for (t in 1:size(times)) {
      for (a in 1:size(tx_times)) {
        if (times[t] == tx_times[a]) {
          tx_index[a] = t;
        }
      }
    }

    // Cell population dynamics before the first fraction of radiation
    N_t0[1] = V_t0;
    N_t0[2] = 0;
    N_t0[3] = 0;

    N_hat[1:tx_index[1],] = to_matrix(cell_population_dynamics_num(times[1:tx_index[1]], N_t0, params));
    
    // Cell population dynamics between fractions of radiation
    for (admin in 1:size(tx_times) - 1) {
      // 'instantaneous' effect of each rt administration on cell population counts
      N_t0[1] = N_hat[tx_index[admin], 1]*exp(-(alpha1 + alpha2)*tx_doses[admin]);
      N_t0[2] = N_hat[tx_index[admin], 2] +
        N_hat[tx_index[admin], 1]*(1 - exp(-alpha1*tx_doses[admin]));
      N_t0[3] = N_hat[tx_index[admin], 3] +
        N_hat[tx_index[admin], 1]*(1 - exp(-alpha2*tx_doses[admin]));  
      
      N_hat[tx_index[admin] + 1:tx_index[admin + 1],] = 
        to_matrix(cell_population_dynamics_num(times[tx_index[admin] + 1:tx_index[admin + 1]], N_t0, params));
    } 
    
    // Cell population dynamics from final fraction of radiation to final tumor volume measurement
    N_t0[1] = N_hat[size(tx_times), 1]*exp(-(alpha1 + alpha2)*tx_doses[size(tx_times)]);
    N_t0[2] = N_hat[size(tx_times), 2] +
      N_hat[size(tx_times), 1]*exp(-alpha1*tx_doses[size(tx_times)]);
    N_t0[3] = N_hat[size(tx_times), 3] +
      N_hat[size(tx_times), 1]*exp(-alpha2*tx_doses[size(tx_times)]);  
    
    N_hat[tx_index[size(tx_index)] + 1:size(times), ] = 
      to_matrix(cell_population_dynamics_num(times[tx_index[size(tx_index)] + 1:size(times)], N_t0, params));

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
  int<lower = 1> n_admin;                                   // Number of treatment administrations
  int<lower = 1> n_mice;                                    // Number of mice
  int<lower = 0> n_control_mice;                            // Number of control mice
  real<lower = 0> tx_times[n_admin];                        // Times of treatment administrations
  real<lower = 0> tx_doses[n_admin];                        // Treatment doses
  int<lower = 0> schedule_start[n_mice];                    // index of first treatment administration for each mouse
  int<lower = 0> schedule_end[n_mice];                      // index of last treatment administration for each mouse
  int<lower = 0> start[n_mice];                             // index of first observation for each mouse
  int<lower = 0> end[n_mice];                               // index of last observation for each mouse
  int<lower = 0> n_pre_tx_times[n_mice];                    // Number of time points before the first treatment
  int<lower = 0> n_post_tx_times[n_mice];                   // Number of time points after the last treatment
  real<lower = 0> t0[n_mice];                               // Initial time
  real<lower = 0> V_t0[n_mice];                             // Tumour volumes at initial time
  int<lower = 0> obs[n_obs];                                // Observation row numbers
  real<lower = 0> t_obs[n_obs];                             // Measurement times
  real<lower = 0> V_obs[n_obs];                             // Measured tumor volumes
  real<lower = 0> times[n_times];                           // Times (observations and fractions)
  
  int<lower = 0, upper = 1> run_estimation;                 // A switch to evaluate the likelihood

}

transformed data {
  
  real x_r[0];
  int x_i[0];
  real logV_obs[n_obs];                                   // Log(measured tumor volumes)

  logV_obs = log10(V_obs);
}

parameters {
  
  // population-level parameters on log scale
  real loggamma_hat;                              // Growth rate
  real logK_hat;                                  // carrying capacity
  real logdelta_hat;                              // Death rate
  real logalpha1_hat;                             // instantaneous transfer rate at radiation times
  real logalpha2_hat;                             // instantaneous transfer rate at radiation times

  vector[5] logtheta[n_mice];                     // Individual level parameter values on log scale
  
  corr_matrix[5] rho;                             // Correlation matrix - within individuals, 
  vector<lower = 0>[5] omega;                     // Standard deviations of individual-level parameters
  real logsigma;                                  // Residual standard deviation on log scale
  
}

transformed parameters {
  
  vector[5] log_theta_hat;                        // Population level parameters on log scale, condensed into vector
  cov_matrix[5] Omega;                            // Covariance matrix
  real<lower = 0> sigma;                          // residual standard deviation
  
  // individual level parameters
  real<lower = 0> gamma[n_mice];                          
  real<lower = 0> K[n_mice];
  real<lower = 0> delta[n_mice];
  real<lower = 0> alpha1[n_mice];
  real<lower = 0> alpha2[n_mice];
  
  real<lower = 0> V_hat[n_times];                         // Estimated tumor volume at measurement and treatment times
  real<lower = 0> V_hat_obs[n_obs];                       // Estimated tumor volume at volume measurement times
  
  log_theta_hat[1] = loggamma_hat;  
  log_theta_hat[2] = logK_hat;
  log_theta_hat[3] = logdelta_hat;
  log_theta_hat[4] = logalpha1_hat;
  log_theta_hat[5] = logalpha2_hat;
  
  sigma = pow(10, logsigma);
  
  Omega = quad_form_diag(rho, omega);                     // diag_matrix(omega) * rho * diag_matrix(omega)
  
  for (i in 1:n_mice) {
    gamma[i] = pow(10, logtheta[i, 1]);
    K[i] = pow(10, logtheta[i, 2]);
    delta[i] = pow(10, logtheta[i, 3]);
    alpha1[i] = pow(10, logtheta[i, 4]);
    alpha2[i] = pow(10, logtheta[i, 5]);
    
    if (i <= n_control_mice)
      V_hat[start[i]:end[i]] =
        control_tumor_volume_dynamics(times[start[i]:end[i]], V_t0[i],
                                      gamma[i], K[i], delta[i], x_r, x_i);
    else
      V_hat[start[i]:end[i]] =
        tumor_volume_dynamics(times[start[i]:end[i]],
                              tx_times[schedule_start[i]:schedule_end[i]],
                              tx_doses[schedule_start[i]:schedule_end[i]],
                              n_pre_tx_times[i], n_post_tx_times[i], t0[i], V_t0[i],
                              gamma[i], K[i], delta[i], alpha1[i], alpha2[i], x_r, x_i);
    
  }

  V_hat_obs = V_hat[obs];
  
}

model {
  
  loggamma_hat ~ normal(-2.5, 1);
  logK_hat ~  normal(2.75, 1);
  logdelta_hat ~ normal(-2, 1);
  logalpha1_hat ~ normal(-2, 1);
  logalpha2_hat ~ normal(-2, 1);

  omega ~ cauchy(0, 1);                                   // Half-cauchy prior (see Bayesian Data Analysis, Gelman et al. p.130)
  rho ~ lkj_corr(1);                                      // Approximately equivalent to uniform distribution on correlations

  logsigma ~ normal(-1.5, 1);


  // compute likelihood with given parameters
  if (run_estimation == 1) {
    logV_obs ~ normal(log10(V_hat_obs), sigma);
  }
}


generated quantities {
  
  real<lower = 0> V_pred[n_obs];
  vector[n_obs] log_likelihood;
  
  // Random number generator drawing from normal distribution
  for (ob in 1:n_obs) {
    V_pred[ob] = pow(10, (normal_rng(log10(V_hat_obs[ob]), sigma)));
    log_likelihood[ob] = normal_rng(log10(V_hat_obs[ob]), sigma);
  }
  
}
