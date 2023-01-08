// Stan code for DIPG7 tumour growth and response to RT and drug. 
// Simple model considers the dynamics of proliferating tumour cells (P),
// senescent cells (Q) and cells undergoing apoptosis (D).
// dP/dt = gamma*log(K/P)*P - (beta_q + beta_a)*P, 
// dQ/dt = beta_q*P - beta_a*Q,
// dD/dt = beta_a*P - delta*D.
// Tumour growth modelled as Gompertz growth and exponential decay of
// senescent and dying cells based on DIPG7 data from 
// Ashley Vardon and data from St Jude hospital

// Ielyaas Cloete (18 November 2022)

functions {
  
  // Computes the solution of the system of ordinary differential equations, given the model parameters,
  // initial conditions and requested solution times
  
  real[] cell_population_dynamics(real t,		              // Time at which derivatives are evaluated
	                                real[] N,		            // System state at which derivatives are evaluated
                                  real[] params,	        // Parameters for system
                                  real[] x_r,	            // Real constants for system that are not fit
                                  int[] x_i) {            // Integer constants for system that are not fit

    real gamma = params[1];                               // Growth rate
    real k = params[2];                                   // Carrying capacity
    real beta_q = params[3];                              // Rate of entry into senescence
    real beta_a = params[4];                              // Rate of entry into apoptosis
    real delta = params[5];                               // Rate of removal of dead cells
    
    real dN_dt[3];

    dN_dt[1] = gamma*log(k/N[1])*N[1] - (beta_q + beta_a)*N[1];       // Proliferating cells
    dN_dt[2] = beta_q*N[1] - beta_a*N[2];                             // Quiescent cells
    dN_dt[3] = beta_a*N[1]*N[2] - delta*N[3];                         // Dying cells
    
    return dN_dt;
  }
  
  // Computes the tumor volume dynamics for untreated mice
  
  real[] control_tumor_volume_dynamics(real[] times, real t0, real V_t0,
                                       real gamma, real k, real beta_q, real beta_a, real delta,
                                       real[] x_r, int[] x_i,
                                       real rel_tol, real abs_tol, int max_num_steps) {
    
    real N_t0[3];                                         // Initial state
    real params[5];                                       // ODE parameters
    real N_t[size(times), 3];
    matrix[size(times), 3] N_hat;                         // Predicted cell population numbers
    real V_hat[size(times)];                              // Predicted tumor volume
    
    params[1] = gamma;
    params[2] = k;
    params[3] = beta_q;
    params[4] = beta_a;
    params[5] = delta;
    
    N_t0[1] = V_t0;
    N_t0[2] = 0;
    N_t0[3] = 0;
    
    N_t = integrate_ode_bdf(cell_population_dynamics,//integrate_ode_rk45(cell_population_dynamics,
                             N_t0, t0,
                             times,
                             params, x_r, x_i,
                             rel_tol, abs_tol, max_num_steps);
    
    for (compartment in 1:size(N_t0)) {
      N_hat[1:size(N_t), compartment] = to_vector(N_t[1:size(N_t), compartment]);
    }
    
    // Sum over cell types to get tumor volume
    for (t in 1:size(times)) {
      V_hat[t] = N_hat[t, 1] + N_hat[t, 2] + N_hat[t, 3];
    }

    return V_hat;
  }

  // Computes the tumor volume dynamics following radiation given the cell population numbers
  // computed by the cell_population_dynamics() function using a recursive approach to handle
  // multiple doses (handles event histories with observation and dosing events)
  
  real[] drug_rt_tumor_volume_dynamics(real[] times, real[] tx_times, real[] tx_doses,
                                       real[] tx_modalities, int n_pre_tx_times, int n_post_tx_times, real t0,
                                       real V_t0, real gamma, real k, real beta_q, real beta_a,
                                       real delta, real alpha, real[] x_r, int[] x_i,
                                       real rel_tol, real abs_tol, int max_num_steps) {
    
    real N_t0[3];                                         // Initial state
    real params[5];                                       // ODE parameters
    int tx_index[size(tx_times)];
    real N_t_1[1, 3];                                     // Cell number for 2 time point
    real N_t_2[2, 3];                                     // Cell number for 3 time points
    real N_t_pre[n_pre_tx_times, 3];
    real N_t_post[n_post_tx_times, 3];
    matrix[size(times), 3] N_hat;                         // Predicted cell population numbers
    real drug;                                            // Response to drug administration
    real RT;                                              // Response to radiation
    real V_hat[size(times)];                              // Predicted tumor volume
    
    params[1] = gamma;
    params[2] = k;
    params[3] = beta_q;
    params[4] = beta_a;
    params[5] = delta;

    // Determine indices of treatment administrations in times array
    for (t in 1:size(times)) {
      for (a in 1:size(tx_times)) {
        if (times[t] == tx_times[a]) {
          tx_index[a] = t;
        }
      }
    }

    // Cell population dynamics before the first treatment administration
    N_t0[1] = V_t0;
    N_t0[2] = 0;
    N_t0[3] = 0;
    
    N_t_pre = integrate_ode_bdf(cell_population_dynamics,//integrate_ode_rk45(cell_population_dynamics,
                                 N_t0, t0,
                                 times[1:tx_index[1]],
                                 params, x_r, x_i,
                                 rel_tol, abs_tol, max_num_steps);
    
    for (compartment in 1:size(N_t0)) {
      N_hat[1:tx_index[1], compartment] = to_vector(N_t_pre[1:size(N_t_pre), compartment]);
    }
    
    // Cell population dynamics between treatment administrations
    for (admin in 1:size(tx_times) - 1) {

      // Drug administration
      if (tx_modalities[admin] == 1) {
        drug = exp(-alpha*tx_doses[admin]);
        
        // Set inital tumor volume and drug concentration (for ODE initial conditions) to final volume following
        // previous administration
        N_t0[1] = N_hat[tx_index[admin], 1]*drug;
        N_t0[2] = N_hat[tx_index[admin], 2]*drug;
        N_t0[3] = N_hat[tx_index[admin], 3]*(1 - drug);
      
      // Radiation administration
      } else if (tx_modalities[admin] == 2) {
        // Calculate DNA damage following radiation
        RT = exp(-alpha*tx_doses[admin]);
      
        // Set inital tumor volume and drug concentration (for ODE initial conditions) to final volume following
        // previous fraction
        N_t0[1] = N_hat[tx_index[admin], 1]*RT;
        N_t0[2] = N_hat[tx_index[admin], 2]*(1 - RT);
        N_t0[3] = N_hat[tx_index[admin], 3];
      }
      
      if (size(times[tx_index[admin] + 1:tx_index[admin + 1]]) == 1) {
        
        N_t_1 = integrate_ode_bdf(cell_population_dynamics,//integrate_ode_rk45(cell_population_dynamics,
                                   N_t0, times[tx_index[admin]],
                                   times[tx_index[admin] + 1:tx_index[admin + 1]],
                                   params, x_r, x_i,
                                   rel_tol, abs_tol, max_num_steps);
        
        for (compartment in 1:size(N_t0)) {
          N_hat[tx_index[admin] + 1:tx_index[admin + 1], compartment] =
            to_vector(N_t_1[1:size(N_t_1), compartment]);
        }
      
      } else if (size(times[tx_index[admin] + 1:tx_index[admin + 1]]) == 2) {
        
        N_t_2 = integrate_ode_bdf(cell_population_dynamics,//integrate_ode_rk45(cell_population_dynamics,
                                   N_t0, times[tx_index[admin]],
                                   times[tx_index[admin] + 1:tx_index[admin + 1]],
                                   params, x_r, x_i,
                                   rel_tol, abs_tol, max_num_steps);
        
        for (compartment in 1:size(N_t0)) {
          N_hat[tx_index[admin] + 1:tx_index[admin + 1], compartment] =
            to_vector(N_t_2[1:size(N_t_2), compartment]);
        }
      }
    }
    // Drug administration
    if (tx_modalities[size(tx_modalities)] == 1) {
      drug = exp(-alpha*tx_doses[size(tx_doses)]);
      
      // Set inital tumor volume and drug concentration (for ODE initial conditions) to final volume following
      // previous administration
      N_t0[1] = N_hat[tx_index[size(tx_index)], 1]*drug;
      N_t0[2] = N_hat[tx_index[size(tx_index)], 2]*drug;
      N_t0[3] = N_hat[tx_index[size(tx_index)], 3]*(1 - drug);
      
    // Radiation administration
    } else if (tx_modalities[size(tx_modalities)] == 2) {
      // Calculate DNA damage following final radiation administration
      RT = exp(-alpha*tx_doses[size(tx_doses)]);
      
      // Cell population dynamics from final fraction of radiation to final tumor volume measurement
      N_t0[1] = N_hat[tx_index[size(tx_index)], 1]*RT;
      N_t0[2] = N_hat[tx_index[size(tx_index)], 2]*(1 - RT);
      N_t0[3] = N_hat[tx_index[size(tx_index)], 3];
    }
    
    N_t_post = integrate_ode_bdf(cell_population_dynamics,//integrate_ode_rk45(cell_population_dynamics,
                                  N_t0, times[tx_index[size(tx_index)]],
                                  times[tx_index[size(tx_index)] + 1:size(times)],
                                  params, x_r, x_i,
                                  rel_tol, abs_tol, max_num_steps);
    
    if (size(N_t_post) == 1) {
      for (compartment in 1:size(N_t0)) {
        N_hat[tx_index[size(tx_index)] + 1, compartment] =
          to_vector(to_matrix(N_t_post))[compartment];
      }
    } else {
      for (compartment in 1:size(N_t0)) {
        N_hat[tx_index[size(tx_index)] + 1:size(times), compartment] =
          to_vector(N_t_post[1:size(N_t_post), compartment]);
      }
    }
    
    // Sum over cell types to get tumor volume
    for (t in 1:size(times)) {
      V_hat[t] = N_hat[t, 1] + N_hat[t, 2] + N_hat[t, 3];
      if (V_hat[t] <= 0.0) {
        V_hat[t] = 1e-10;
      }
    }

    return V_hat;
  }

}

data {
  
  int<lower = 1> n_times;                                 // Number of times (observations and fractions)
  int<lower = 1> n_obs;                                   // Number of observations
  int<lower = 0> n_admin;                                 // Number of treatment administrations
  int<lower = 1> n_mice;                                  // Number of mice
  int<lower = 0> n_control_mice;                          // Number of control/untreated mice
  real<lower = 0> tx_times[n_admin];                      // Times of treatment administrations
  real<lower = 0> tx_doses[n_admin];                      // Treatment doses
  real<lower = 0> tx_modalities[n_admin];                 // Treatment modalities
  int<lower = 0> tx_schedule_start[n_mice];               // Index of first treatment admin for each mouse
  int<lower = 0> tx_schedule_end[n_mice];                 // Index of last treatment admin for each mouse
  int<lower = 0> start[n_mice];                           // Index of first time for each mouse
  int<lower = 0> end[n_mice];                             // Index of last time for each mouse
  int<lower = 0> n_pre_tx_times[n_mice];                  // Number of time points before the first treatment
  int<lower = 0> n_post_tx_times[n_mice];                 // Number of time points after the last treatment
  real<lower = 0> t0[n_mice];                             // Initial time
  real<lower = 0> V_t0[n_mice];                           // Tumor volumes at initial time
  int<lower = 0> obs[n_obs];                              // Observation row numbers
  real<lower = 0> t_obs[n_obs];                           // Measurement times
  real<lower = 0> V_obs[n_obs];                           // Measured tumor volumes
  real<lower = 0> times[n_times];                         // Times (observations and fractions)
  
  int<lower = 0, upper = 1> run_estimation;               // A switch to evaluate the likelihood
  
  real rel_tol;
  real abs_tol;
  int max_num_steps;
  
}

transformed data {

  real x_r[0];
  int x_i[0];
  
  real logV_obs[n_obs];
  
  logV_obs = log(V_obs);                                  // Log(measured tumor volumes)
  
}

parameters {
  
  real<lower = 0> gamma_hat;                              // Growth rate
  real<lower = 0> k_hat;                                  // Carrying capacity
  real<lower = 0> beta_q_hat;                             // Treatment dependent senescence
  real<lower = 0> beta_a_hat;                             // Treatment dependent apoptosis
  real<lower = 0> delta_hat;                              // Death rate
  real<lower = 0> alpha_hat;                              // Dose dependent treatment damage
  
  
  corr_matrix[6] rho;                                     // Correlation matrix
  vector<lower = 0>[6] omega;                             // Standard deviations for inter-individual variability
  real<lower = 0> sigma;                                  // Residual standard deviation
  vector[6] logtheta[n_mice];                             // Individual level parameter values
  
}

transformed parameters {

  vector<lower = 0>[6] theta_hat;                         // Population level means of parameter values
  cov_matrix[6] Omega;                                    // Covariance matrix
  real<lower = 0> gamma[n_mice];
  real<lower = 0> k[n_mice];
  real<lower = 0> beta_q[n_mice];
  real<lower = 0> beta_a[n_mice];
  real<lower = 0> delta[n_mice];
  real<lower = 0> alpha[n_mice];

  real<lower = 0> V_hat[n_times];                         // Output from the ODE solver
  real<lower = 0> V_hat_obs[n_obs];                       // Estimated tumor volume at volume measurement times
  
  // Population level paramters
  theta_hat[1] = gamma_hat;                               // Popultaion level mean growth rate
  theta_hat[2] = k_hat;                                   // Popultaion level mean k
  theta_hat[3] = beta_q_hat;                              // Popultaion level mean beta_quiescent
  theta_hat[4] = beta_a_hat;                              // Popultaion level mean beta_apoptotic
  theta_hat[5] = delta_hat;                               // Popultaion level mean death rate
  theta_hat[6] = alpha_hat;                               // Popultaion level mean alpha
  
  Omega = quad_form_diag(rho, omega);                     // diag_matrix(omega) * rho * diag_matrix(omega)
  
  for (i in 1:n_mice) {

    gamma[i] = exp(logtheta[i, 1]);
    k[i] = exp(logtheta[i, 2]);
    beta_q[i] = exp(logtheta[i, 3]);
    beta_a[i] = exp(logtheta[i, 4]);
    delta[i] = exp(logtheta[i, 5]);
    alpha[i] = exp(logtheta[i, 6]);

    if (i <= n_control_mice) {
      V_hat[start[i]:end[i]] =
        control_tumor_volume_dynamics(times[start[i]:end[i]], t0[i], V_t0[i],
                                      gamma[i], k[i], beta_q[i], beta_a[i], delta[i],
                                      x_r, x_i,
                                      rel_tol, abs_tol, max_num_steps);
    } else {
      V_hat[start[i]:end[i]] =
        drug_rt_tumor_volume_dynamics(times[start[i]:end[i]],
                                      tx_times[tx_schedule_start[i]:tx_schedule_end[i]],
                                      tx_doses[tx_schedule_start[i]:tx_schedule_end[i]],
                                      tx_modalities[tx_schedule_start[i]:tx_schedule_end[i]],
                                      n_pre_tx_times[i], n_post_tx_times[i], t0[i], V_t0[i],
                                      gamma[i], k[i], beta_q[i], beta_a[i], delta[i], alpha[i],
                                      x_r, x_i,
                                      rel_tol, abs_tol, max_num_steps); 
    }
  }
  
  V_hat_obs = V_hat[obs];
}

model {
  
  // Priors
  // Note that priors are also informed by constraints in parameters block
  gamma_hat ~ normal(0.003, 1);
  k_hat ~ normal(562.3, 1);
  beta_q_hat ~ normal(0.01, 1);
  beta_a_hat ~ normal(0.01, 1);
  delta_hat ~ normal(5, 1);
  alpha_hat ~ normal(0.9, 1);
  
  // Half-cauchy prior (see Bayesian Data Analysis, Gelman et al. p.130)
  omega ~ cauchy(0, 1);//2.5                            // Standard deviations
  // Approximately equivalent to uniform distribution on correlations
  rho ~ lkj_corr(2);                                      // Correlation matrix with LKJ prior     
  sigma ~ cauchy(0, 2);//1                                // Residual standard deviation
  
  // Inter-individual variability
  logtheta ~ multi_normal(log(theta_hat), Omega);         // Log transformed individual level parameter values
  
  // Likelihood
  if (run_estimation == 1) {
    logV_obs ~ normal(log(V_hat_obs), sigma);             // Normal error
  }
  
}

generated quantities {
  
  vector[6] logtheta_pred[n_mice];
  real<lower = 0> gamma_pred[n_mice];
  real<lower = 0> k_pred[n_mice];
  real<lower = 0> beta_q_pred[n_mice];
  real<lower = 0> beta_a_pred[n_mice];
  real<lower = 0> delta_pred[n_mice];
  real<lower = 0> alpha_pred[n_mice];
  
  real<lower = 0> V_hat_pred_times[n_times];
  real<lower = 0> V_hat_pred[n_obs];
  real<lower = 0> V_cond[n_obs];
  real<lower = 0> V_pred[n_obs];
  vector[n_obs] log_likelihood;
  
  for (i in 1:n_mice) {
    
    logtheta_pred[i] = multi_normal_rng(log(theta_hat), Omega);
    gamma_pred[i] = exp(logtheta_pred[i, 1]);
    k_pred[i] = exp(logtheta_pred[i, 2]);
    beta_q_pred[i] = exp(logtheta_pred[i, 3]);
    beta_a_pred[i] = exp(logtheta_pred[i, 4]);
    delta_pred[i] = exp(logtheta_pred[i, 5]);
    alpha_pred[i] = exp(logtheta_pred[i, 6]);
    
    if (i <= n_control_mice) {
      V_hat_pred_times[start[i]:end[i]] =
        control_tumor_volume_dynamics(times[start[i]:end[i]], t0[i], V_t0[i],
                                      gamma_pred[i], k_pred[i], beta_q_pred[i], beta_a_pred[i],
                                      delta_pred[i], x_r, x_i,
                                      rel_tol, abs_tol, max_num_steps);
    } else {
      V_hat_pred_times[start[i]:end[i]] =
        drug_rt_tumor_volume_dynamics(times[start[i]:end[i]],
                                      tx_times[tx_schedule_start[i]:tx_schedule_end[i]],
                                      tx_doses[tx_schedule_start[i]:tx_schedule_end[i]],
                                      tx_modalities[tx_schedule_start[i]:tx_schedule_end[i]],
                                      n_pre_tx_times[i], n_post_tx_times[i], t0[i], V_t0[i],
                                      gamma_pred[i], k_pred[i], beta_q_pred[i], beta_a_pred[i],
                                      delta_pred[i], alpha_pred[i], x_r, x_i,
                                      rel_tol, abs_tol, max_num_steps);
    }
  }
  
  V_hat_pred = V_hat_pred_times[obs];
  
  // Random number generator drawing from normal distribution
  for (ob in 1:n_obs) {
    // Individual predictions
    V_cond[ob] = exp(normal_rng(log(V_hat_obs[ob]), sigma));
    // Log likelihood
    log_likelihood[ob] = exp(normal_lpdf(logV_obs[ob] | log(V_hat_obs[ob]), sigma));
    // Population predictions
    V_pred[ob] = exp(normal_rng(log(V_hat_pred[ob]), sigma));
  }
  
}
